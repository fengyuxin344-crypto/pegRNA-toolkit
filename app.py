"""
pegRNA Library Designer — Local web app
=======================================
Run:  python app.py
Then open http://localhost:5000 in your browser.

A local tool (single machine, small team). Wraps engine.py in a simple
three-step web UI:
  1. Upload mutation list + gene  -> download PRIDICT batch input
  2. (you run PRIDICT2.0 yourself; the app shows the exact command)
  3. Upload PRIDICT summary        -> download self-targeting library
"""

import io
import json
import os
import glob
import time
import uuid
import tempfile
import threading
import traceback

from flask import Flask, request, jsonify, send_file, Response
import pandas as pd

import engine
try:
    import optiprime_integration as optiprime
except Exception:
    optiprime = None  # OptiPrime is optional; only needed if user selects it

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 32 * 1024 * 1024  # 32 MB uploads

# --- NGS Analysis tab (Blueprint) ---
import analysis
app.register_blueprint(analysis.bp)

# in-memory holding of the last generated files (fine for a local single-user tool)
_STATE = {}
# OptiPrime runs in a completely separate state so the two engine sub-tabs never
# share intermediates (batch.csv, summary, library, QC) — they are independent
# pipelines from the user's point of view.
_STATE_OP = {}


@app.route("/")
def index():
    return Response(INDEX_HTML, mimetype="text/html")


@app.route("/api/lookup_variants", methods=["POST"])
def api_lookup_variants():
    """Turn literature-style lines into unambiguous ClinVar candidates.
    Lines are free-form: 'MAPT, NP_005901.2, P301L', 'TERT promoter, -146C>T',
    'GBA1, NP_000148.2, N409S (previously refer as N370S)'. The USER picks."""
    try:
        text = request.form.get("pasted", "").strip()
        if not text:
            return jsonify({"error": "Paste one variant per line, e.g. 'MAPT, P301L'."}), 400
        lines = [l for l in (ln.strip() for ln in text.splitlines()) if l]
        api_key = request.form.get("ncbi_api_key", "").strip() or None
        results, errors = engine.lookup_variant_table(lines, api_key=api_key)
        return jsonify({"ok": True, "results": results, "errors": errors})
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


@app.route("/api/build_input", methods=["POST"])
def api_build_input():
    """Stage 1: mutation list + gene -> PRIDICT batch CSV.

    Accepts EITHER an uploaded file (.xlsx/.csv) OR pasted text (`pasted` form
    field). Pasted text can be:
      - full CSV with a header row (must include gDNA_mutation, optional gene)
      - simple lines 'GENE, mutation' or 'GENE mutation'
      - simple lines with just a mutation (then the Gene field is used)
    Multi-gene mode kicks in automatically if a gene column/second field exists.
    Mutations may be c. notation OR protein changes (e.g. L483P).
    """
    try:
        mut_df = None
        pasted = request.form.get("pasted", "").strip()
        f = request.files.get("file")

        if pasted:
            mut_df = _parse_pasted_mutations(pasted)
        elif f:
            if f.filename.lower().endswith(".csv"):
                mut_df = pd.read_csv(f)
            else:
                mut_df = _read_mutation_excel(f)
        else:
            return jsonify({"error": "Upload a file or paste a mutation list."}), 400

        if "gDNA_mutation" not in mut_df.columns:
            return jsonify({"error": "Need a 'gDNA_mutation' column (or pasted "
                            "lines like 'GENE, L483P'). "
                            f"Found: {list(mut_df.columns)[:8]}."}), 400

        # multi-gene mode if a 'gene' column is present
        if "gene" in mut_df.columns and mut_df["gene"].notna().any():
            batch_df, report = engine.build_pridict_input_multigene(mut_df)
        else:
            # gene is optional now: genomic (g.), rsID and versioned-transcript
            # inputs carry their own coordinate system and don't need a gene.
            # Rows that DO need a gene (protein / canonical opt-in) without one
            # are reported in `skipped` with a clear message, not silently run.
            gene = request.form.get("gene", "").strip()
            batch_df, report = engine.build_pridict_input(mut_df, gene)

        _STATE["batch_csv"] = batch_df.to_csv(index=False)
        _STATE["source_map"] = report.get("source_map", {})
        return jsonify({"ok": True, "report": report})
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


def _read_mutation_excel(path_or_file):
    """Read a mutation table from an .xlsx that may have MULTIPLE sheets AND a
    header that isn't on the first row (template sheets have a title + notes
    above the real header).

    Strategy: for every sheet, find the row containing 'gDNA_mutation', re-parse
    with it as the header, and keep only sheets that actually have DATA rows
    (a non-empty gDNA_mutation value). Return the first sheet with real data —
    so empty template tabs (header but no rows) are skipped in favour of the tab
    the user actually filled in. Falls back to the first sheet with a header, or
    the first sheet, so the caller still raises a clear error if nothing fits."""
    xl = pd.ExcelFile(path_or_file)
    first_with_header = None
    for name in xl.sheet_names:
        raw = xl.parse(name, header=None)
        for i in range(min(len(raw), 15)):
            rowvals = [str(v).strip() for v in raw.iloc[i].tolist()]
            if "gDNA_mutation" in rowvals:
                df = xl.parse(name, header=i).dropna(how="all")
                # normalise: strip whitespace-only cells to NaN
                df = df.replace(r"^\s*$", pd.NA, regex=True)
                if first_with_header is None:
                    first_with_header = df
                # does this sheet have at least one real mutation value?
                if "gDNA_mutation" in df.columns and df["gDNA_mutation"].notna().any():
                    return df.dropna(subset=["gDNA_mutation"])
                break  # header found on this sheet; move to next sheet
    if first_with_header is not None:
        return first_with_header
    return xl.parse(xl.sheet_names[0])


def _parse_pasted_mutations(text):
    """Turn pasted text into a mutation dataframe.

    Handles three shapes:
      1. CSV with a header row containing 'gDNA_mutation'
      2. lines 'GENE, mutation' / 'GENE mutation' (comma or whitespace separated)
      3. lines with just a mutation
    Returns a dataframe with at least 'gDNA_mutation', optionally 'gene'.
    """
    import io as _io
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    # strip trailing '# ...' comments (e.g. the rsID appended by 'Use selected
    # coordinates') and drop lines that are entirely a comment
    lines = [ln.split("#", 1)[0].strip() for ln in lines]
    lines = [ln for ln in lines if ln]
    if not lines:
        raise ValueError("Pasted text is empty.")

    # case 1: looks like CSV with a header
    header = lines[0].lower()
    if "gdna_mutation" in header or ("gene" in header and "," in lines[0]):
        return pd.read_csv(_io.StringIO(text))

    # case 2/3: parse line by line
    genes, txs, muts = [], [], []
    for ln in lines:
        # split on comma first, else whitespace
        parts = [p.strip() for p in (ln.split(",") if "," in ln else ln.split())]
        parts = [p for p in parts if p]
        if len(parts) >= 3:
            # GENE, TRANSCRIPT, MUTATION
            genes.append(parts[0]); txs.append(parts[1]); muts.append(parts[-1])
        elif len(parts) == 2:
            # GENE, MUTATION
            genes.append(parts[0]); txs.append(None); muts.append(parts[-1])
        else:
            genes.append(None); txs.append(None); muts.append(parts[0])
    df = pd.DataFrame({"gene": genes, "transcript": txs, "gDNA_mutation": muts})
    # drop columns that are entirely empty
    if df["gene"].isna().all():
        df = df.drop(columns=["gene"])
    if df["transcript"].isna().all():
        df = df.drop(columns=["transcript"])
    return df




def _attach_source_variant(lib_df, state=None):
    """Add a `source_variant` column to the library by joining the Stage-1
    source map on sequence_name. This carries the ORIGINAL user input (plus its
    canonical [REF>ALT]) independently of PRIDICT, so Stage-3 QC can confirm each
    construct installs the mutation that was actually requested. No-op if there
    is no source map (e.g. an uploaded summary with unknown provenance).
    `state` selects which pipeline's source map to use (defaults to _STATE)."""
    smap = (state if state is not None else _STATE).get("source_map") or {}
    if not smap or "sequence_name" not in lib_df.columns:
        return
    lib_df["source_variant"] = lib_df["sequence_name"].map(smap).fillna("")


def _bystander_directions():
    """Read the bystander direction choice from the request form and return a
    tuple of directions to generate: () none, ('intro',), ('revert',), or both."""
    v = request.form.get("bystander", "").strip().lower()
    if v in ("intro",):
        return ("intro",)
    if v in ("revert",):
        return ("revert",)
    if v in ("both", "intro_revert", "intro+revert"):
        return ("intro", "revert")
    if v in ("1", "true", "on", "yes"):   # legacy checkbox -> revert (old default)
        return ("revert",)
    return ()   # none / off / empty


def _extract_stage1_inputs():
    """Read the CURRENT request into a plain dict so Stage 1 can run OUTSIDE the
    request context (e.g. in a background thread). Uploaded files are saved to a
    temp dir and referenced by path; nothing here touches `request` afterwards."""
    dna_files = [f for f in request.files.getlist("files")
                 if getattr(f, "filename", "") and f.filename.lower().endswith(".dna")]
    inp = {
        "dna_paths": [],
        "pasted": request.form.get("pasted", "").strip(),
        "gene": request.form.get("gene", "").strip(),
        "directions": _bystander_directions(),
        "window": int(request.form.get("window_codons", 2)),
        "byst_home": request.form.get("pridict_home", "~/PRIDICT2").strip() or "~/PRIDICT2",
        "mut_path": None, "mut_name": None,
    }
    if dna_files:
        tmpdir = tempfile.mkdtemp()
        for f in dna_files:
            p = os.path.join(tmpdir, os.path.basename(f.filename))
            f.save(p)
            inp["dna_paths"].append(p)
    else:
        f = request.files.get("file")
        if f and f.filename:
            tmpdir = tempfile.mkdtemp()
            p = os.path.join(tmpdir, os.path.basename(f.filename))
            f.save(p)
            inp["mut_path"], inp["mut_name"] = p, f.filename
    return inp


def _stage1_batch_core(inp):
    """Build the PRIDICT batch dataframe from an already-extracted inputs dict
    (see _extract_stage1_inputs). No `request` access — safe in a thread."""
    # ---- SnapGene mode (most reliable) ----
    if inp["dna_paths"]:
        paths = inp["dna_paths"]
        batch_df, report = engine.build_pridict_input_from_snapgene(paths)
        directions = inp["directions"]
        if directions:
            window = inp["window"]
            byst_df, byst_report = engine.build_bystander_from_snapgene(
                paths, window_codons=window, directions=directions,
                pridict_home=inp["byst_home"])
            if len(byst_df):
                batch_df = pd.concat([batch_df, byst_df], ignore_index=True)
            report.setdefault("source_map", {}).update(byst_report.get("source_map", {}))
            report["bystander"] = {
                "n_variants": byst_report["n_variants_generated"],
                "per_file": byst_report["per_file"],
                "n_skipped": byst_report["n_skipped"],
                "skipped": byst_report["skipped"],
                "window_codons": window,
                "directions": list(directions),
            }
        return batch_df, report

    # ---- mutation-list mode ----
    if inp["pasted"]:
        mut_df = _parse_pasted_mutations(inp["pasted"])
    elif inp["mut_path"]:
        name = inp["mut_name"] or inp["mut_path"]
        mut_df = (pd.read_csv(inp["mut_path"]) if name.lower().endswith(".csv")
                  else _read_mutation_excel(inp["mut_path"]))
    else:
        raise ValueError("Upload SnapGene .dna files, or paste/upload a mutation list.")

    if "gDNA_mutation" not in mut_df.columns:
        raise ValueError("Need a 'gDNA_mutation' column (or pasted lines like "
                         f"'GENE, L483P'). Found: {list(mut_df.columns)[:8]}.")

    if "gene" in mut_df.columns and mut_df["gene"].notna().any():
        return engine.build_pridict_input_multigene(mut_df)
    if not inp["gene"]:
        raise ValueError("Enter a gene symbol/Ensembl ID, or include the gene in "
                         "each pasted line / a 'gene' column.")
    return engine.build_pridict_input(mut_df, inp["gene"])


def _stage1_batch():
    """Request-based Stage 1 (used by the individual /api/build_input endpoints)."""
    return _stage1_batch_core(_extract_stage1_inputs())



@app.route("/api/build_input_snapgene", methods=["POST"])
def api_build_input_snapgene():
    """Stage 1 (SnapGene mode): annotated .dna files -> PRIDICT batch CSV.

    The most reliable input path — mutations are read straight from the
    annotated features, so no Ensembl / transcript alignment is involved.
    """
    try:
        files = request.files.getlist("files")
        if not files:
            return jsonify({"error": "Upload one or more SnapGene .dna files."}), 400

        import tempfile
        tmpdir = tempfile.mkdtemp()
        paths = []
        for f in files:
            if not f.filename.lower().endswith(".dna"):
                continue
            p = os.path.join(tmpdir, os.path.basename(f.filename))
            f.save(p)
            paths.append(p)
        if not paths:
            return jsonify({"error": "No .dna files found in the upload."}), 400

        directions = _bystander_directions()
        window = int(request.form.get("window_codons", 2))

        batch_df, report = engine.build_pridict_input_from_snapgene(paths)

        if directions:
            byst_home = request.form.get("pridict_home", "~/PRIDICT2").strip() or "~/PRIDICT2"
            byst_df, byst_report = engine.build_bystander_from_snapgene(
                paths, window_codons=window, directions=directions, pridict_home=byst_home)
            # combine standard + bystander inputs
            if len(byst_df):
                batch_df = pd.concat([batch_df, byst_df], ignore_index=True)
            report.setdefault("source_map", {}).update(byst_report.get("source_map", {}))
            report["bystander"] = {
                "n_variants": byst_report["n_variants_generated"],
                "per_file": byst_report["per_file"],
                "n_skipped": byst_report["n_skipped"],
                "skipped": byst_report["skipped"],
                "window_codons": window,
                "directions": list(directions),
            }

        _STATE["batch_csv"] = batch_df.to_csv(index=False)
        _STATE["source_map"] = report.get("source_map", {})
        return jsonify({"ok": True, "report": report})
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


@app.route("/api/download_batch")
def api_download_batch():
    if "batch_csv" not in _STATE:
        return "No batch file generated yet.", 404
    return send_file(io.BytesIO(_STATE["batch_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="pridict_batch_input.csv")


@app.route("/api/run_pridict", methods=["POST"])
def api_run_pridict():
    """Stage 2 automation (Layer 1): run PRIDICT2.0 on the batch.csv from Stage 1.

    BLOCKING — waits for PRIDICT to finish, which can take a long time. The
    browser fetch may time out even though the run continues on disk; when that
    happens the summary still lands in ~/PRIDICT2/predictions/.
    Uses the in-memory batch.csv from Stage 1, or an uploaded file if provided.
    """
    try:
        f = request.files.get("file")
        if f and f.filename:
            batch_text = f.read().decode("utf-8", "replace")
        elif "batch_csv" in _STATE:
            batch_text = _STATE["batch_csv"]
        else:
            return jsonify({"error": "No batch.csv yet — generate it in Step 1 "
                            "first, or upload one here."}), 400

        home = request.form.get("pridict_home", "~/PRIDICT2").strip() or "~/PRIDICT2"
        env = request.form.get("conda_env", "pridict2").strip() or "pridict2"
        cores = int(request.form.get("cores", 3))
        snum = int(request.form.get("summarize_number", 10))
        run_name = request.form.get("run_name", "").strip()

        result = engine.run_pridict(
            batch_text, pridict_home=home, conda_env=env,
            cores=cores, summarize="K562", summarize_number=snum,
            run_name=run_name)

        # stash the summary so Step 3 can build the library without a re-upload
        _STATE["pridict_summary_csv"] = result["summary_text"]
        return jsonify({"ok": True,
                        "n_rows": result["n_rows"],
                        "run_name": result["run_name"],
                        "summary_path": result["summary_path"],
                        "input_path": result["input_path"],
                        "archived_to": result["archived_to"]})
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


@app.route("/api/run_all_optiprime", methods=["POST"])
def api_run_all_optiprime():
    """ONE CLICK for the OptiPrime sub-tab: mutation input -> editseq ->
    OptiPrime scoring -> assembled library -> (Step 3 result). Fully independent
    of the PRIDICT pipeline — everything is stashed in _STATE_OP so the two
    engine tabs never share intermediates.
    """
    if optiprime is None:
        return jsonify({"error": "OptiPrime integration module not available on "
                        "the server."}), 500
    try:
        # ---- Stage 1: mutation input -> batch.csv + per-codon silent options ----
        batch_df, codon_blocks, s1_report = _stage1_optiprime_core(
            _extract_stage1_inputs())
        _STATE_OP["batch_csv"] = batch_df.to_csv(index=False)
        _STATE_OP["source_map"] = s1_report.get("source_map", {})

        # ---- Stage 2: score with OptiPrime -> PRIDICT-style summary ----
        home = request.form.get("optiprime_home", "~/optiprime-src").strip() or "~/optiprime-src"
        env = request.form.get("optiprime_env", "optiprime").strip() or "optiprime"
        summary_df, op_report = optiprime.run_batch_optiprime(
            _STATE_OP["batch_csv"], optiprime_home=home, conda_env=env,
            codon_blocks_map=codon_blocks, settings=_optiprime_settings())
        _STATE_OP["summary_csv"] = summary_df.to_csv(index=False)

        # ---- Stage 3: assemble library (shared engine fn, isolated data) ----
        top_n = int(request.form.get("top_n", 10))
        lib_df, s3_report = engine.build_library(
            summary_df, params={"top_n_per_mutation": top_n})
        _attach_source_variant(lib_df, state=_STATE_OP)
        _STATE_OP["library_csv"] = lib_df.to_csv(index=False)

        return jsonify({
            "ok": True,
            "stage1": s1_report,
            "optiprime": {"n_scored_pegRNAs": op_report["n_scored_pegRNAs"],
                          "n_runs": op_report["n_runs"],
                          "n_input_sequences": op_report["n_input_sequences"],
                          "n_skipped": op_report["n_skipped"],
                          "n_bystander_input": op_report.get("n_bystander_input", 0),
                          "n_bystander_skipped": op_report.get("n_bystander_skipped", 0),
                          "skipped_reason_counts": op_report.get("skipped_reason_counts", {}),
                          "warning": op_report.get("warning", ""),
                          "skipped": op_report["skipped"][:20]},
            "stage3": s3_report,
        })
    except Exception as e:
        print("\n=== OptiPrime run_all_optiprime ERROR ===", flush=True)
        traceback.print_exc()
        print("=== end error ===\n", flush=True)
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


# OptiPrime scores the same pegRNA differently depending on the experimental
# context (cell type, MMR status, editor, scaffold), so these are read from the
# form rather than left at the library's defaults. A mismatch here surfaces as a
# plausible-but-wrong efficiency, never as an error, which is why the values are
# validated against OptiPrime's own vocabularies instead of passed through.
_OP_VOCAB = {
    "cell_type": ["HEK293T", "HeLa", "A549", "HAP1", "K562", "U2OS", "DLD1",
                  "MDA-MB-231", "NIH3T3"],
    "pe_type": ["PE2", "PE4"],
    "cas9_type": ["PE2-Cas9", "PEmax-Cas9", "PE6e-Cas9", "PE6f-Cas9", "PE6g-Cas9"],
    "rt_name": ["PE2-RT", "PE6a-RT", "PE6b-RT", "PE6c-RT", "PE6d-RT"],
    "scaffold": ["SpCas9_OG", "OG_F+E", "BlpI_F+E", "GC_F+E"],
    "motif": ["none", "tevoPreQ1"],
    "cas9_pam": ["SpNGG", "SpNG", "SpNRCH", "SpNRTH", "SpNRRH", "SpG", "SpRY",
                 "SpVRQR", "SpVQR", "SpVRER"],
}


def _optiprime_settings():
    """Read + validate the OptiPrime experimental context from the request."""
    raw = request.form.get("optiprime_settings", "").strip()
    if not raw:
        return {}
    try:
        d = json.loads(raw)
    except Exception:
        return {}
    out = {}
    for k, allowed in _OP_VOCAB.items():
        v = str(d.get(k, "")).strip()
        if v in allowed:
            out[k] = v
    try:
        t = float(d.get("time", 4.0))
        if 0.5 <= t <= 30:
            out["time"] = t
    except (TypeError, ValueError):
        pass
    return out


def _stage1_optiprime_core(inp):
    """Stage 1 for the OptiPrime path.

    Unlike PRIDICT2, OptiPrime does NOT want one input row per silent-bystander
    combination — it wants one run carrying the per-codon options in
    `edit_segments` and combines them itself. So the batch here stays main-edits
    only, and the bystander choice is returned separately as a codon block map.
    """
    if not inp["dna_paths"]:
        # Mutation-list mode has no CDS annotation, so there is no reading frame
        # to derive synonymous codons from. Main edits only, said out loud.
        batch_df, report = _stage1_batch_core({**inp, "directions": ()})
        if inp["directions"]:
            report["bystander_note"] = (
                "Silent bystanders need a reading frame, which only the "
                "SnapGene CDS annotation provides. Main edits only for this run.")
        return batch_df, {}, report

    paths = inp["dna_paths"]
    batch_df, report = engine.build_pridict_input_from_snapgene(paths)
    blocks = {}
    if inp["directions"]:
        blocks, brep = engine.build_codon_blocks_from_snapgene(
            paths, window_codons=inp["window"], directions=inp["directions"])
        report["bystander"] = {
            "mode": "OptiPrime native (edit_segments combinatorics)",
            "window_codons": brep["window_codons"],
            "n_sequences": brep["n_sequences"],
            "n_combinations": brep["n_combinations"],
            "skipped": brep["skipped"],
            "directions": list(inp["directions"]),
        }
    return batch_df, blocks, report


def _run_optiprime_pipeline(job_id, inp, params):
    """Background worker for the OptiPrime one-click, updating _JOBS so
    /api/progress reports a stage-based percentage. OptiPrime itself is a single
    batched call, so its stage shows a slowly-advancing bar rather than an exact
    per-mutation count."""
    job = _JOBS[job_id]
    try:
        job.update(stage="input", percent=8, message="Building input from mutations...")
        batch_df, codon_blocks, s1_report = _stage1_optiprime_core(inp)
        _STATE_OP["batch_csv"] = batch_df.to_csv(index=False)
        _STATE_OP["source_map"] = s1_report.get("source_map", {})
        if len(batch_df) == 0:
            raise ValueError("No sequences were generated from the input.")

        job.update(stage="optiprime", percent=30,
                   message="Searching protospacers & scoring with OptiPrime "
                           "(this can take several minutes)...")
        summary_df, op_report = optiprime.run_batch_optiprime(
            _STATE_OP["batch_csv"], optiprime_home=params["home"],
            conda_env=params["env"], codon_blocks_map=codon_blocks,
            settings=params.get("settings") or {})
        _STATE_OP["summary_csv"] = summary_df.to_csv(index=False)

        job.update(stage="assemble", percent=92,
                   message="Assembling self-targeting library...")
        lib_df, s3_report = engine.build_library(
            summary_df, params={"top_n_per_mutation": params["top_n"]})
        _attach_source_variant(lib_df, state=_STATE_OP)
        _STATE_OP["library_csv"] = lib_df.to_csv(index=False)

        job.update(status="done", stage="done", percent=100, message="Done.",
                   result={"ok": True, "stage1": s1_report,
                           "optiprime": {"n_scored_pegRNAs": op_report["n_scored_pegRNAs"],
                                         "n_runs": op_report["n_runs"],
                                         "n_input_sequences": op_report["n_input_sequences"],
                                         "n_skipped": op_report["n_skipped"],
                                         "n_bystander_input": op_report.get("n_bystander_input", 0),
                                         "n_bystander_skipped": op_report.get("n_bystander_skipped", 0),
                                         "skipped_reason_counts": op_report.get("skipped_reason_counts", {}),
                                         "warning": op_report.get("warning", ""),
                                         "skipped": op_report["skipped"][:20]},
                           "stage3": s3_report})
    except Exception as e:
        print("\n=== OptiPrime async ERROR ===", flush=True)
        traceback.print_exc()
        print("=== end error ===\n", flush=True)
        job.update(status="error", stage="error",
                   error=f"{type(e).__name__}: {e}",
                   trace=traceback.format_exc())


@app.route("/api/run_all_optiprime_async", methods=["POST"])
def api_run_all_optiprime_async():
    """Start the OptiPrime one-click in the background; poll /api/progress."""
    if optiprime is None:
        return jsonify({"error": "OptiPrime integration module not available."}), 500
    try:
        inp = _extract_stage1_inputs()
        params = {
            "home": request.form.get("optiprime_home", "~/optiprime-src").strip() or "~/optiprime-src",
            "env": request.form.get("optiprime_env", "optiprime").strip() or "optiprime",
            "settings": _optiprime_settings(),
            "top_n": int(request.form.get("top_n", 10)),
        }
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 400

    job_id = uuid.uuid4().hex[:12]
    _JOBS[job_id] = {"status": "running", "stage": "input", "percent": 0,
                     "total": 0, "done": 0, "message": "Starting..."}
    t = threading.Thread(target=_run_optiprime_pipeline,
                         args=(job_id, inp, params), daemon=True)
    t.start()
    return jsonify({"ok": True, "job_id": job_id})


@app.route("/api/download_batch_op")
def api_download_batch_op():
    if "batch_csv" not in _STATE_OP:
        return "No OptiPrime batch yet.", 404
    return send_file(io.BytesIO(_STATE_OP["batch_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="batch_optiprime.csv")


@app.route("/api/download_optiprime_summary")
def api_download_optiprime_summary():
    if "summary_csv" not in _STATE_OP:
        return "No OptiPrime summary yet.", 404
    return send_file(io.BytesIO(_STATE_OP["summary_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="optiprime_summary.csv")


@app.route("/api/download_library_op")
def api_download_library_op():
    if "library_csv" not in _STATE_OP:
        return "No OptiPrime library yet.", 404
    return send_file(io.BytesIO(_STATE_OP["library_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="pegRNA_library_optiprime.csv")


@app.route("/api/run_qc_op", methods=["POST"])
def api_run_qc_op():
    if "library_csv" not in _STATE_OP:
        return jsonify({"error": "No OptiPrime library to QC yet."}), 400
    try:
        lib_df = pd.read_csv(io.StringIO(_STATE_OP["library_csv"]))
        clean_df, dropped_df, report = engine.qc_library(lib_df)
        _STATE_OP["qc_csv"] = clean_df.to_csv(index=False)
        _STATE_OP["qc_dropped_csv"] = dropped_df.to_csv(index=False)
        return jsonify({"ok": True, "report": report})
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


@app.route("/api/download_qc_op")
def api_download_qc_op():
    if "qc_csv" not in _STATE_OP:
        return "No OptiPrime QC results yet.", 404
    return send_file(io.BytesIO(_STATE_OP["qc_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="pegRNA_library_optiprime_QC_passed.csv")


@app.route("/api/download_qc_dropped_op")
def api_download_qc_dropped_op():
    if "qc_dropped_csv" not in _STATE_OP:
        return "No OptiPrime QC results yet.", 404
    return send_file(io.BytesIO(_STATE_OP["qc_dropped_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="pegRNA_library_optiprime_QC_dropped.csv")


@app.route("/api/download_pridict_summary")
def api_download_pridict_summary():
    if "pridict_summary_csv" not in _STATE:
        return "No PRIDICT summary yet.", 404
    return send_file(io.BytesIO(_STATE["pridict_summary_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="pridict_summary_K562.csv")


@app.route("/api/run_all", methods=["POST"])
def api_run_all():
    """ONE CLICK: Stage 1 (build input) -> Stage 2 (run PRIDICT) -> Stage 3
    (assemble library), with zero manual download/drop/copy in between.

    Every intermediate is still stashed in _STATE, so the individual
    'Download batch.csv' / 'Download summary' / 'Download library' links keep
    working exactly as before — this just chains the three steps for you.
    """
    try:
        # ---- Stage 1: mutation input -> batch.csv ----
        batch_df, s1_report = _stage1_batch()
        _STATE["batch_csv"] = batch_df.to_csv(index=False)
        _STATE["source_map"] = s1_report.get("source_map", {})

        # ---- Stage 2: run PRIDICT (blocking) ----
        home = request.form.get("pridict_home", "~/PRIDICT2").strip() or "~/PRIDICT2"
        env = request.form.get("conda_env", "pridict2").strip() or "pridict2"
        cores = int(request.form.get("cores", 3))
        snum = int(request.form.get("summarize_number", 10))
        run_name = request.form.get("run_name", "").strip()

        pr = engine.run_pridict(_STATE["batch_csv"], pridict_home=home, conda_env=env,
                                cores=cores, summarize="K562", summarize_number=snum,
                                run_name=run_name)
        _STATE["pridict_summary_csv"] = pr["summary_text"]

        # ---- Stage 3: PRIDICT summary -> self-targeting library ----
        top_n = int(request.form.get("top_n", 10))
        summary_df = pd.read_csv(io.StringIO(pr["summary_text"]))
        lib_df, s3_report = engine.build_library(
            summary_df, params={"top_n_per_mutation": top_n})
        _attach_source_variant(lib_df)
        _STATE["library_csv"] = lib_df.to_csv(index=False)

        return jsonify({
            "ok": True,
            "stage1": s1_report,
            "pridict": {"n_rows": pr["n_rows"], "run_name": pr["run_name"],
                        "archived_to": pr["archived_to"]},
            "stage3": s3_report,
        })
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


# ---------------------------------------------------------------------------
# Background job runner for the one-click pipeline: the HTTP request returns a
# job_id immediately, PRIDICT runs in a thread, and the frontend polls
# /api/progress/<job_id> for a real percentage (so the browser never "times
# out" on the long PRIDICT step). Progress is measured from the number of
# per-sequence prediction files PRIDICT has written into predictions/.
# ---------------------------------------------------------------------------

_JOBS = {}


def _count_pridict_done(pred_dir, since):
    """How many per-sequence predictions PRIDICT has written for this run =
    top-level .csv files in predictions/ (excluding the summary + _archive)
    created at/after `since`."""
    if not pred_dir or not os.path.isdir(pred_dir):
        return 0
    n = 0
    for p in glob.glob(os.path.join(pred_dir, "*.csv")):
        base = os.path.basename(p)
        if "batch_summary" in base:
            continue
        try:
            if os.path.getmtime(p) >= since - 2:
                n += 1
        except OSError:
            pass
    return n


def _run_pipeline(job_id, inp, params):
    """Background worker: Stage 1 -> Stage 2 (PRIDICT) -> Stage 3, updating the
    shared job record so /api/progress can report status + percent."""
    job = _JOBS[job_id]
    try:
        # ---- Stage 1: build input ----
        job.update(stage="input", message="Building PRIDICT input...")
        batch_df, s1_report = _stage1_batch_core(inp)
        _STATE["batch_csv"] = batch_df.to_csv(index=False)
        _STATE["source_map"] = s1_report.get("source_map", {})
        total = len(batch_df)
        if total == 0:
            raise ValueError("No sequences were generated from the input.")

        # ---- Stage 2: run PRIDICT (blocking here, but we're in a thread) ----
        home = params["home"]
        job.update(stage="pridict", total=total, done=0,
                   pred_dir=os.path.join(os.path.expanduser(home), "predictions"),
                   pridict_start=time.time(),
                   message=f"Running PRIDICT on {total} sequences...")
        pr = engine.run_pridict(
            _STATE["batch_csv"], pridict_home=home, conda_env=params["env"],
            cores=params["cores"], summarize="K562",
            summarize_number=params["snum"], run_name=params["run_name"])
        _STATE["pridict_summary_csv"] = pr["summary_text"]

        # ---- Stage 3: assemble library ----
        job.update(stage="assemble", message="Assembling self-targeting library...")
        summary_df = pd.read_csv(io.StringIO(pr["summary_text"]))
        lib_df, s3_report = engine.build_library(
            summary_df, params={"top_n_per_mutation": params["top_n"]})
        _attach_source_variant(lib_df)
        _STATE["library_csv"] = lib_df.to_csv(index=False)

        job.update(status="done", stage="done", percent=100,
                   message="Done.",
                   result={"ok": True, "stage1": s1_report,
                           "pridict": {"n_rows": pr["n_rows"],
                                       "run_name": pr["run_name"],
                                       "archived_to": pr["archived_to"]},
                           "stage3": s3_report})
    except Exception as e:
        job.update(status="error", stage="error",
                   error=f"{type(e).__name__}: {e}",
                   trace=traceback.format_exc())


@app.route("/api/run_all_async", methods=["POST"])
def api_run_all_async():
    """Start the one-click pipeline in the background. Returns a job_id
    immediately; poll /api/progress/<job_id> for status + percent."""
    try:
        inp = _extract_stage1_inputs()
        params = {
            "home": request.form.get("pridict_home", "~/PRIDICT2").strip() or "~/PRIDICT2",
            "env": request.form.get("conda_env", "pridict2").strip() or "pridict2",
            "cores": int(request.form.get("cores", 3)),
            "snum": int(request.form.get("summarize_number", 10)),
            "run_name": request.form.get("run_name", "").strip(),
            "top_n": int(request.form.get("top_n", 10)),
        }
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 400

    job_id = uuid.uuid4().hex[:12]
    _JOBS[job_id] = {"status": "running", "stage": "input", "percent": 0,
                     "total": 0, "done": 0, "message": "Starting..."}
    t = threading.Thread(target=_run_pipeline, args=(job_id, inp, params), daemon=True)
    t.start()
    return jsonify({"ok": True, "job_id": job_id})


@app.route("/api/progress/<job_id>", methods=["GET"])
def api_progress(job_id):
    """Report a running job's status and a real percentage during PRIDICT."""
    job = _JOBS.get(job_id)
    if not job:
        return jsonify({"error": "Unknown job id (it may have been cleared)."}), 404

    # live percentage during the PRIDICT step, from prediction files on disk
    if job.get("status") == "running" and job.get("stage") == "pridict":
        done = _count_pridict_done(job.get("pred_dir"), job.get("pridict_start", 0))
        total = max(1, job.get("total", 1))
        job["done"] = done
        job["percent"] = min(99, int(100 * done / total))  # 100 only when truly done

    out = {k: job.get(k) for k in ("status", "stage", "percent", "total",
                                   "done", "message")}
    if job.get("status") == "done":
        out["result"] = job.get("result")
    elif job.get("status") == "error":
        out["error"] = job.get("error")
        out["trace"] = job.get("trace")
    return jsonify(out)


@app.route("/api/build_library", methods=["POST"])
def api_build_library():
    """Stage 3: PRIDICT summary -> self-targeting library."""
    try:
        f = request.files.get("file")
        if f and f.filename:
            summary_df = pd.read_csv(f)
        elif "pridict_summary_csv" in _STATE:
            # use the summary produced by the Step 2 "Run PRIDICT" button
            summary_df = pd.read_csv(io.StringIO(_STATE["pridict_summary_csv"]))
        else:
            return jsonify({"error": "Upload the PRIDICT K562 summary CSV, or run "
                            "PRIDICT in Step 2 first."}), 400

        top_n = int(request.form.get("top_n", 10))
        params = {"top_n_per_mutation": top_n}
        # library_size kept as an optional overall cap if provided
        if request.form.get("library_size"):
            params["library_size"] = int(request.form.get("library_size"))
        lib_df, report = engine.build_library(summary_df, params=params)
        _attach_source_variant(lib_df)
        _STATE["library_csv"] = lib_df.to_csv(index=False)
        return jsonify({"ok": True, "report": report})
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


@app.route("/api/download_library")
def api_download_library():
    if "library_csv" not in _STATE:
        return "No library generated yet.", 404
    return send_file(io.BytesIO(_STATE["library_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="pegRNA_library.csv")


@app.route("/api/run_qc", methods=["POST"])
def api_run_qc():
    """Stage 3-QC: check the assembled library. Reports pass/fail per construct
    without removing anything. Uses the just-built library if no file is given."""
    try:
        f = request.files.get("file")
        if f:
            lib_df = pd.read_csv(f)
        elif "library_csv" in _STATE:
            lib_df = pd.read_csv(io.StringIO(_STATE["library_csv"]))
        else:
            return jsonify({"error": "Build a library first, or upload a library CSV to QC."}), 400

        clean_df, dropped_df, report = engine.qc_library(lib_df)
        # clean CSV = only constructs that passed every check
        _STATE["qc_csv"] = clean_df.to_csv(index=False)
        # dropped CSV = everything removed, with a QC_drop_reason column
        _STATE["qc_dropped_csv"] = dropped_df.to_csv(index=False)
        return jsonify({"ok": True, "report": report})
    except Exception as e:
        return jsonify({"error": f"{type(e).__name__}: {e}",
                        "trace": traceback.format_exc()}), 500


@app.route("/api/download_qc")
def api_download_qc():
    if "qc_csv" not in _STATE:
        return "No QC results yet.", 404
    return send_file(io.BytesIO(_STATE["qc_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="pegRNA_library_QC_passed.csv")


@app.route("/api/download_qc_dropped")
def api_download_qc_dropped():
    if "qc_dropped_csv" not in _STATE:
        return "No QC results yet.", 404
    return send_file(io.BytesIO(_STATE["qc_dropped_csv"].encode()),
                     mimetype="text/csv", as_attachment=True,
                     download_name="pegRNA_library_QC_dropped.csv")


# Frontend is served from the sibling file for readability
with open(os.path.join(os.path.dirname(__file__), "index.html")) as _fh:
    INDEX_HTML = _fh.read()


def _pick_port(candidates=(5050, 5000, 5001, 8000, 8080)):
    """Pick the first free port. 5050 first because macOS AirPlay grabs 5000."""
    import socket
    for port in candidates:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            s.bind(("127.0.0.1", port))
            s.close()
            return port
        except OSError:
            s.close()
    return 5050  # everything busy? try anyway and let Flask report it


if __name__ == "__main__":
    import threading
    import webbrowser

    port = _pick_port()
    url = f"http://localhost:{port}"
    print(f"\n  pegRNA Library Designer running at  {url}")
    print("  (Leave this window open while you use the tool. Close it or press "
          "Ctrl+C to stop.)\n")

    # open the browser a moment after the server comes up
    threading.Timer(1.5, lambda: webbrowser.open(url)).start()

    app.run(host="127.0.0.1", port=port, debug=False)
