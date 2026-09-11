"""
analysis.py — Analysis tab backend (Flask Blueprint).

Wraps the self-targeting NGS pipeline behind a few HTTP endpoints so the whole
thing runs from the browser:

    plate map (.xlsx)  ->  merge per-well FASTQs (cat, per editor+replicate)
                       ->  write config/samples.csv + set library_name
                       ->  snakemake --cores N   (current env, no --use-conda)
                       ->  aggregate.py          ->  <date>_<library>_df.xlsx

Design notes / hard-won lessons baked in:
  * FASTQs are tens of GB -> never uploaded; the browser only sends the small
    plate map. Everything else is addressed by local paths on this machine.
  * merge uses `cat` (gzip streams concatenate) -> avoids the macOS BSD-zcat
    `.Z` failure entirely.
  * snakemake is called directly with --cores N, NO --use-conda (the current
    environment already has cutadapt/biopython/pandas).
  * isolated from the design app: its own _ANALYSIS_JOBS dict and url prefix.

Wire-up in app.py (one line, after `app = Flask(__name__)`):
    import analysis
    app.register_blueprint(analysis.bp)
"""

import os
import glob
import time
import uuid
import threading
import subprocess
import traceback
from pathlib import Path
from collections import defaultdict

import pandas as pd
from flask import Blueprint, request, jsonify, send_file

bp = Blueprint("analysis", __name__)

# isolated from the design app's _JOBS
_ANALYSIS_JOBS = {}


# --------------------------------------------------------------------------- #
# pure helpers (unit-testable without Flask)
# --------------------------------------------------------------------------- #
def _clean_editor(editor):
    """control_A -> control, PEmax_B -> PEmax; otherwise underscores -> dashes.

    Strips a trailing single-letter replicate-group suffix (the '_A'/'_B'/...)
    if present, since that letter is redundant with the replicate number.
    """
    editor = str(editor).strip()
    head, _, tail = editor.rpartition("_")
    if head and len(tail) == 1 and tail.isalpha():
        return head
    return editor.replace("_", "-")


def parse_plate_map(xlsx_path, library_name):
    """Group per-well rows into biological samples by (editor, replicate).

    Expects columns: replicate, editor, r1_file, r2_file (a 'sample' column,
    if present, is ignored - it is per-well, not per biological sample).
    Returns a list of dicts:
        {sample, editor, replicate, r1_files:[...], r2_files:[...]}
    """
    df = pd.read_excel(xlsx_path)
    df.columns = [str(c).strip().lower() for c in df.columns]
    required = {"replicate", "editor", "r1_file", "r2_file"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Plate map is missing columns {sorted(missing)}. "
            f"Found: {list(df.columns)}"
        )
    df = df.dropna(subset=["editor", "replicate", "r1_file", "r2_file"])

    groups = defaultdict(lambda: {"r1": [], "r2": []})
    order = []
    for _, row in df.iterrows():
        editor = str(row["editor"]).strip()
        rep = str(row["replicate"]).strip()
        key = (editor, rep)
        if key not in groups:
            order.append(key)
        g = groups[key]
        r1, r2 = str(row["r1_file"]).strip(), str(row["r2_file"]).strip()
        if r1 not in g["r1"]:
            g["r1"].append(r1)
        if r2 not in g["r2"]:
            g["r2"].append(r2)

    samples = []
    for editor, rep in order:
        base = _clean_editor(editor)
        sample = f"{library_name}-{base}-rep{rep}"
        g = groups[(editor, rep)]
        if len(g["r1"]) != len(g["r2"]):
            raise ValueError(f"{sample}: R1/R2 well count mismatch "
                             f"({len(g['r1'])} vs {len(g['r2'])})")
        samples.append({
            "sample": sample, "editor": base, "replicate": rep,
            "r1_files": g["r1"], "r2_files": g["r2"],
        })

    # sample names must be unique and underscore-free (Snakemake wildcards)
    names = [s["sample"] for s in samples]
    if len(names) != len(set(names)):
        raise ValueError("Derived sample names are not unique - check the plate map.")
    bad = [n for n in names if "_" in n]
    if bad:
        raise ValueError(f"Sample names contain '_' (breaks Snakemake): {bad}")
    return samples


def write_samples_csv(samples, library_name, out_path):
    """7-row config/samples.csv pointing at the *merged* file names."""
    rows = []
    for s in samples:
        rows.append({
            "sample": s["sample"], "library": library_name,
            "replicate": s["replicate"], "editor": s["editor"],
            "r1_file": f"{s['sample']}_R1.fastq.gz",
            "r2_file": f"{s['sample']}_R2.fastq.gz",
        })
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return out_path


def write_combine_csv(samples, raw_dir_name, out_path):
    """Per-well combine-input CSV (one row per well) for combine_techreps_cat.py.

    r1_file/r2_file are prefixed with the raw dir *name* so the combine script,
    run from the pipeline dir, resolves them.
    """
    rows = []
    for s in samples:
        for r1, r2 in zip(s["r1_files"], s["r2_files"]):
            rows.append({
                "sample": s["sample"],
                "r1_file": f"{raw_dir_name}/{r1}",
                "r2_file": f"{raw_dir_name}/{r2}",
            })
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return out_path


def set_config_library_name(config_csv, library_name):
    """Set library_name in the parameter/value config CSV (in place)."""
    df = pd.read_csv(config_csv)
    if "parameter" not in df.columns or "value" not in df.columns:
        raise ValueError("config CSV must have 'parameter' and 'value' columns.")
    mask = df["parameter"] == "library_name"
    if not mask.any():
        raise ValueError("No 'library_name' row in the config CSV.")
    df.loc[mask, "value"] = library_name
    df.to_csv(config_csv, index=False)


def read_config_values(config_csv):
    """Return {parameter: value} from the config CSV (values as strings)."""
    df = pd.read_csv(config_csv)
    if "parameter" not in df.columns or "value" not in df.columns:
        raise ValueError("config CSV must have 'parameter' and 'value' columns.")
    out = {}
    for _, r in df.iterrows():
        p = r["parameter"]
        if pd.isna(p):
            continue
        v = r["value"]
        out[str(p).strip()] = "" if pd.isna(v) else str(v).strip()
    return out


def set_config_params(config_csv, params):
    """Update multiple parameters in the config CSV (in place).

    Only existing rows are updated; unknown keys are ignored (so the UI can't
    inject arbitrary parameters). Empty-string values are skipped, so a blank
    field never wipes a real config value.
    """
    df = pd.read_csv(config_csv)
    if "parameter" not in df.columns or "value" not in df.columns:
        raise ValueError("config CSV must have 'parameter' and 'value' columns.")
    changed = []
    for key, val in params.items():
        if val is None or str(val).strip() == "":
            continue
        mask = df["parameter"] == key
        if mask.any():
            df.loc[mask, "value"] = str(val)
            changed.append(key)
    df.to_csv(config_csv, index=False)
    return changed


def merge_group(raw_dir, out_dir, sample, r1_files, r2_files):
    """cat per-well .gz into one merged pair. Returns (r1_out, r2_out)."""
    raw, out = Path(raw_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for tag, files in (("R1", r1_files), ("R2", r2_files)):
        missing = [f for f in files if not (raw / f).exists()]
        if missing:
            raise FileNotFoundError(
                f"{sample} {tag}: {len(missing)} well file(s) missing under "
                f"{raw_dir}, e.g. {missing[:2]}"
            )
        dst = out / f"{sample}_{tag}.fastq.gz"
        with open(dst, "wb") as fh:
            subprocess.run(["cat"] + [str(raw / f) for f in files],
                           check=True, stdout=fh)
    return out / f"{sample}_R1.fastq.gz", out / f"{sample}_R2.fastq.gz"


# --------------------------------------------------------------------------- #
# background worker
# --------------------------------------------------------------------------- #
def _run_analysis(job_id, params, plate_map_path):
    job = _ANALYSIS_JOBS[job_id]
    try:
        pipeline_dir = Path(params["pipeline_dir"]).expanduser()
        raw_dir = Path(params["raw_dir"]).expanduser()
        library = params["library_name"]
        date = params["date"]
        cores = params["cores"]
        template = params["template"]

        if not pipeline_dir.is_dir():
            raise NotADirectoryError(f"pipeline_dir not found: {pipeline_dir}")
        if not raw_dir.is_dir():
            raise NotADirectoryError(f"raw FASTQ dir not found: {raw_dir}")

        fastq_out = pipeline_dir / "fastqs"
        config_csv = pipeline_dir / "config" / "self_targeting_config.csv"
        samples_csv = pipeline_dir / "config" / "samples.csv"
        analysis_dir = pipeline_dir / "analysis"

        # ---- stage 1: parse plate map ----
        job.update(stage="prepare", message="Parsing plate map...", percent=2)
        samples = parse_plate_map(plate_map_path, library)
        n = len(samples)
        job.update(n_samples=n,
                   message=f"{n} samples: {[s['sample'] for s in samples]}")

        # ---- stage 2: merge FASTQs (cat) ----
        job.update(stage="merge", done=0, total=n,
                   message=f"Merging per-well FASTQs into {n} samples...")
        for i, s in enumerate(samples, 1):
            merge_group(raw_dir, fastq_out, s["sample"], s["r1_files"], s["r2_files"])
            job.update(done=i, percent=2 + int(28 * i / n),
                       message=f"Merged {i}/{n}: {s['sample']}")

        # ---- stage 3: samples.csv + config ----
        job.update(stage="config", message="Writing samples.csv + config...", percent=32)
        write_samples_csv(samples, library, samples_csv)
        # library_name always; plus any config params the UI sent
        cfg_updates = dict(params.get("config", {}))
        cfg_updates["library_name"] = library
        set_config_params(config_csv, cfg_updates)

        # ---- stage 4: snakemake (no --use-conda) ----
        job.update(stage="snakemake", done=0, total=n, percent=34,
                   analysis_dir=str(analysis_dir),
                   message=f"Running snakemake --cores {cores}...")
        proc = subprocess.run(
            ["snakemake", "--cores", str(cores)],
            cwd=str(pipeline_dir), capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                "snakemake failed:\n" + (proc.stderr or proc.stdout)[-3000:])

        # ---- stage 5: aggregate ----
        job.update(stage="aggregate", percent=92, message="Aggregating results...")
        out_xlsx = analysis_dir / "summarized" / f"{date}_{library}_df.xlsx"
        agg = subprocess.run(
            ["python", "aggregate.py",
             "--analysis-dir", str(analysis_dir),
             "--template", template,
             "--samples", str(samples_csv),
             "--library-name", library,
             "--date", date,
             "--out", str(out_xlsx)],
            cwd=str(pipeline_dir), capture_output=True, text=True)
        if agg.returncode != 0:
            raise RuntimeError(
                "aggregate.py failed:\n" + (agg.stderr or agg.stdout)[-3000:])

        job.update(status="done", stage="done", percent=100,
                   message="Done.",
                   result={"ok": True, "n_samples": n,
                           "samples": [s["sample"] for s in samples],
                           "output": str(out_xlsx)})
    except Exception as e:
        job.update(status="error", stage="error",
                   error=f"{type(e).__name__}: {e}",
                   trace=traceback.format_exc())
    finally:
        try:
            os.remove(plate_map_path)
        except OSError:
            pass


# --------------------------------------------------------------------------- #
# routes
# --------------------------------------------------------------------------- #
@bp.route("/api/analysis/config", methods=["GET"])
def api_analysis_config():
    """Read current config values from a pipeline dir, so the UI can prefill
    the B/C-tier fields with the real values instead of blank placeholders."""
    pipeline_dir = request.args.get("pipeline_dir", "").strip()
    if not pipeline_dir:
        return jsonify({"error": "Missing pipeline_dir."}), 400
    cfg = Path(pipeline_dir).expanduser() / "config" / "self_targeting_config.csv"
    if not cfg.is_file():
        return jsonify({"error": f"Config not found: {cfg}"}), 404
    try:
        return jsonify({"ok": True, "values": read_config_values(cfg)})
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 400


@bp.route("/api/analysis/prepare", methods=["POST"])
def api_analysis_prepare():
    """Parse the plate map and write the two sample CSVs ONLY.

    No merging, no snakemake. Fast/synchronous. Lets the user inspect the
    derived sample sheets before committing to a full run.
    Writes:
      <pipeline_dir>/config/samples.csv          (7-row, merged file names)
      <pipeline_dir>/samples_<library>_combine.csv  (per-well, for combine script)
    """
    tmp = None
    try:
        f = request.files.get("plate_map")
        if not f:
            return jsonify({"error": "Upload a plate map .xlsx (field 'plate_map')."}), 400
        pipeline_dir = request.form.get("pipeline_dir", "").strip()
        library = request.form.get("library_name", "").strip()
        raw_dir = request.form.get("raw_dir", "").strip()
        if not pipeline_dir or not library:
            return jsonify({"error": "Missing pipeline_dir or library_name."}), 400

        pdir = Path(pipeline_dir).expanduser()
        if not pdir.is_dir():
            return jsonify({"error": f"pipeline_dir not found: {pdir}"}), 400

        tmp = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           f"_platemap_{uuid.uuid4().hex[:8]}.xlsx")
        f.save(tmp)

        samples = parse_plate_map(tmp, library)

        (pdir / "config").mkdir(exist_ok=True)
        samples_csv = pdir / "config" / "samples.csv"
        combine_csv = pdir / f"samples_{library}_combine.csv"
        write_samples_csv(samples, library, samples_csv)
        # raw dir *name* (basename) for the combine csv prefix; fall back to 'fastqs_raw'
        raw_name = os.path.basename(os.path.normpath(raw_dir)) if raw_dir else "fastqs_raw"
        write_combine_csv(samples, raw_name, combine_csv)

        return jsonify({
            "ok": True,
            "n_samples": len(samples),
            "samples": [{"sample": s["sample"], "editor": s["editor"],
                         "replicate": s["replicate"], "wells": len(s["r1_files"])}
                        for s in samples],
            "samples_csv": str(samples_csv),
            "combine_csv": str(combine_csv),
        })
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 400
    finally:
        if tmp:
            try:
                os.remove(tmp)
            except OSError:
                pass


@bp.route("/api/analysis/download_sheet", methods=["GET"])
def api_analysis_download_sheet():
    """Download a generated sample sheet by absolute path (must be a .csv that
    exists and live under the given pipeline dir, to avoid arbitrary reads)."""
    path = request.args.get("path", "")
    if not path or not path.endswith(".csv") or not os.path.isfile(path):
        return jsonify({"error": "File not found."}), 404
    return send_file(path, as_attachment=True,
                     download_name=os.path.basename(path), mimetype="text/csv")


@bp.route("/api/analysis/run_async", methods=["POST"])
def api_analysis_run_async():
    """Start the full analysis in the background. Only the small plate map is
    uploaded; FASTQs stay on disk and are addressed by path."""
    try:
        f = request.files.get("plate_map")
        if not f:
            return jsonify({"error": "Upload a plate map .xlsx (field 'plate_map')."}), 400
        tmp = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            f"_platemap_{uuid.uuid4().hex[:8]}.xlsx")
        f.save(tmp)

        params = {
            "pipeline_dir": request.form.get("pipeline_dir", "").strip(),
            "raw_dir": request.form.get("raw_dir", "").strip(),
            "library_name": request.form.get("library_name", "").strip(),
            "date": request.form.get("date", "").strip(),
            "template": request.form.get("template", "").strip(),
            "cores": int(request.form.get("cores", 4)),
        }
        for k in ("pipeline_dir", "raw_dir", "library_name", "date", "template"):
            if not params[k]:
                os.remove(tmp)
                return jsonify({"error": f"Missing required field: {k}"}), 400

        # optional config overrides: any form field named cfg_<parameter>
        cfg = {}
        for key in request.form:
            if key.startswith("cfg_"):
                v = request.form.get(key, "").strip()
                if v != "":
                    cfg[key[4:]] = v
        # library_design_file, library_size, date_prefix live in A-tier -> map them in too
        cfg.setdefault("date_prefix", params["date"])
        cfg["library_design_file"] = params["template"]
        params["config"] = cfg
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 400

    job_id = uuid.uuid4().hex[:12]
    _ANALYSIS_JOBS[job_id] = {"status": "running", "stage": "prepare",
                              "percent": 0, "done": 0, "total": 0,
                              "message": "Starting..."}
    threading.Thread(target=_run_analysis, args=(job_id, params, tmp),
                     daemon=True).start()
    return jsonify({"ok": True, "job_id": job_id})


@bp.route("/api/analysis/progress/<job_id>", methods=["GET"])
def api_analysis_progress(job_id):
    job = _ANALYSIS_JOBS.get(job_id)
    if not job:
        return jsonify({"error": "Unknown job id."}), 404

    # live percent during the long snakemake step: count per-sample CSVs on disk
    if job.get("status") == "running" and job.get("stage") == "snakemake":
        adir = job.get("analysis_dir")
        n = max(1, job.get("n_samples", 1))
        if adir and os.path.isdir(adir):
            done = len(glob.glob(os.path.join(adir, "*_analysisdf_focused.csv")))
            job["done"] = done
            job["percent"] = 34 + min(57, int(57 * done / n))  # 34..91 during analysis

    out = {k: job.get(k) for k in ("status", "stage", "percent",
                                   "done", "total", "message")}
    if job.get("status") == "done":
        out["result"] = job.get("result")
    elif job.get("status") == "error":
        out["error"] = job.get("error")
        out["trace"] = job.get("trace")
    return jsonify(out)


@bp.route("/api/analysis/download/<job_id>", methods=["GET"])
def api_analysis_download(job_id):
    job = _ANALYSIS_JOBS.get(job_id)
    if not job or job.get("status") != "done":
        return jsonify({"error": "No finished result for this job."}), 404
    path = job["result"]["output"]
    if not os.path.exists(path):
        return jsonify({"error": f"Output missing on disk: {path}"}), 404
    return send_file(path, as_attachment=True, download_name=os.path.basename(path))
