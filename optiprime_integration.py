"""
optiprime_integration.py — run the toolkit's pegRNA design through OptiPrime
instead of PRIDICT2.0, as an independent (Plan A) engine.

The toolkit's Stage-1 output is a set of PRIDICT-style edit strings:

    editseq = "<up><(ref/alt)><down>"        (edit centred, +/-150 bp)

OptiPrime instead wants, per candidate protospacer:

    unedited     : a DNA window with the 20-bp protospacer starting at index 4
                   (PS20_OFFSET), i.e.  NNNN[20-bp protospacer]NGG...
    edit_segments: the EDITED sequence split into segments; single-option
                   segments are fixed, multi-option segments are tried
                   combinatorially (used for silent bystanders).

This module:
  1. find_protospacers()      — enumerate NGG protospacers on both strands whose
                                nick sits upstream of the edit, and lay out the
                                OptiPrime `unedited` window for each.
  2. build_optiprime_runs()   — turn each (editseq, protospacer) into an
                                OptiPrime run dict {name, unedited, edit_segments}.
  3. run_optiprime()          — call DESIGN_PE.py on the JSON (mirrors run_pridict).
  4. optiprime_to_summary()   — translate OptiPrime output back into the exact
                                columns Stage-2 assembly needs, so the rest of the
                                pipeline (assembly, QC, download) is unchanged.

Everything here is deliberately transparent (functions can be printed/tested one
at a time) because the protospacer layout can only be fully validated by running
OptiPrime on your Mac.
"""

import gzip
import json
import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

# OptiPrime lays the protospacer at index 4 of `unedited` (PS20_OFFSET).
PS20_OFFSET = 4
SPACER_LEN = 20
# Official search window: nick may sit up to `searchDist` bp upstream of the edit
# (optiprime-front Utils.js findProtosDirection default searchDist=18).
SEARCH_DIST = 18
# Downstream length kept after the protospacer in the `unedited` window. This is
# start20+71 (a 75-nt window), matching the OptiPrime webserver EXACTLY —
# confirmed byte-for-byte against the webserver's own unedited_seq output. Do
# not change this: the webserver scores RTTs up to ~33 within this same 75-nt
# window, so window length is NOT the lever for score differences.
DOWNSTREAM_LEN = 71
# Per PAM variant, keep at most this many protospacers (nearest-nick first). The
# webserver lets the user pick; in batch mode we auto-take the top few, which is
# what actually reaches the final top-N anyway.
TOP_K_PER_PAMVAR = 3
# Restrict to these PAM variants by default. SpNGG is the standard SpCas9 NGG PAM.
# Set to None to allow every HT-PAMDA variant (SpRY, SpG, ...), like the webserver
# when NGG protospacers are scarce.
DEFAULT_PAM_VARIANTS = ("SpNGG",)

# Stage-2 (self-targeting library) geometry, taken from the reference notebook
# (ALSP_Library_Design_Notebook, cells 22-23) and verified against a PRIDICT2.0
# library: `wide_initial_target` is 99 nt on the protospacer strand with the
# 20 nt genomic protospacer at [9:29], and the spacer is that protospacer with
# its first base SUBSTITUTED by G (never prepended) so it stays exactly 20 nt.
# The notebook hardcodes `spacerlength_with_g = 20` in its 300 bp length budget,
# so a 21 nt spacer silently produces a 301 bp oligo.
WIDE_TARGET_LEN = 99
PROTO_OFFSET_IN_WIDE = 9
SPACER_LEN_WITH_G = 20

_COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# HT-PAMDA table: PAM(4-mer) -> [variant_name, pamda_score]. Copied from the
# OptiPrime webserver (public/HT-PAMDA.json) so protospacer selection matches.
def _load_ht_pamda():
    import json as _json
    here = Path(__file__).resolve().parent
    for cand in (here / "ht_pamda.json", here / "HT-PAMDA.json"):
        if cand.exists():
            return _json.loads(cand.read_text())
    return {}


_HT_PAMDA = _load_ht_pamda()


# ---------------------------------------------------------------------------
# 1. parse a PRIDICT editseq
# ---------------------------------------------------------------------------

# The official PRIDICT2 silent-bystander add-on writes the ALTERED bases in
# LOWERCASE and leaves unchanged bases uppercase, e.g.
#     (CTAGAGCCTCCTG/tTAaAattTatTa)
# An uppercase-only pattern silently rejects every bystander row, so match both
# cases here. Case is not load-bearing for parsing (which base changed is
# recovered by comparing ref to alt), but it IS the marker the rest of the
# toolkit uses to show edits, so `_edit_offsets` preserves that information.
_EDIT_RE = re.compile(r"\(([ACGTacgt-]+)/([ACGTacgt-]+)\)")

# Distance from the start of the OptiPrime window to the nick site
# (4 bp pad + 17 bp into the 20 bp protospacer); see `nick_dist` in
# _find_protos_direction, which measures from window position 21.
NICK_POS_IN_WINDOW = 21


def _edit_offsets(ref: str, alt: str):
    """Positions inside `alt` that actually differ from `ref`.

    For equal-length substitutions this is a positional diff, which reproduces
    exactly what the add-on's lowercasing marks. For indels the whole block is
    treated as changed.
    """
    if len(ref) != len(alt) or "-" in (ref, alt):
        return set(range(len(alt)))
    return {i for i in range(len(alt)) if ref[i].upper() != alt[i].upper()}


def parse_editseq(editseq: str):
    """Split "<up>(ref/alt)<down>" into up/ref/alt/down + edit index.

    Everything is returned UPPERCASE: OptiPrime compares sequences by identity,
    so a stray lowercase base would read as a mismatch. `edit_offsets` keeps the
    information the lowercase carried (which bases are edits).
    """
    m = _EDIT_RE.search(editseq)
    if not m:
        raise ValueError(f"No (ref/alt) edit found in editseq: {editseq[:60]}...")
    ref_raw, alt_raw = m.group(1), m.group(2)
    ref, alt = ref_raw.upper(), alt_raw.upper()
    up = editseq[:m.start()].upper()
    down = editseq[m.end():].upper()
    return {"up": up, "ref": ref, "alt": alt, "down": down,
            "edit_index": len(up),
            "edit_offsets": _edit_offsets(ref_raw, alt_raw)}


def editseq_unedited(parsed) -> str:
    ref = parsed["ref"] if parsed["ref"] != "-" else ""
    return parsed["up"] + ref + parsed["down"]


def editseq_edited(parsed) -> str:
    alt = parsed["alt"] if parsed["alt"] != "-" else ""
    return parsed["up"] + alt + parsed["down"]


def _min_edit(unedited, edited):
    """Port of optiprime-front minEdit: trim the shared prefix/suffix so only the
    changed core remains. Returns (minU, minE, preLen, postLen)."""
    i = 0
    while i < len(unedited) and i < len(edited) and unedited[i] == edited[i]:
        i += 1
    u2, e2 = unedited[i:], edited[i:]
    j = 0
    while j < len(u2) and j < len(e2) and u2[-1 - j] == e2[-1 - j]:
        j += 1
    if j > 0:
        u2, e2 = u2[:-j], e2[:-j]
    return u2, e2, i, j


# ---------------------------------------------------------------------------
# 2. protospacer search  (ported from optiprime-front findProtosDirection)
# ---------------------------------------------------------------------------

def _find_protos_direction(uSeq, eSeq, direction, pam_variants):
    """Port of findProtosDirection. Search one strand for protospacers whose PAM
    (from HT-PAMDA) sits so the nick is within SEARCH_DIST upstream of the edit.
    Returns list of dicts with start20/end20/proto30/unedited/nickDist/pam/... ."""
    minU, minE, preLen, postLen = _min_edit(uSeq, eSeq)
    uDelta = len(uSeq) - len(eSeq) if len(uSeq) > len(eSeq) else 0
    eDelta = len(eSeq) - len(uSeq) if len(eSeq) > len(uSeq) else 0
    preHom = uSeq[:preLen]
    postHom = uSeq[len(uSeq) - postLen:] if postLen else ""
    uLen = len(minU)

    search_start = max(0, preLen - SEARCH_DIST - 21)
    a = min(uLen, 7)
    search = (preHom[search_start:] + minU[:a] + postHom[:max(0, 7 - a)])

    out = []
    # PAM is the 4-mer right after a 24-mer (4 pad + 20 spacer). Match by
    # lookahead so overlapping hits are found, exactly like the JS regex.
    for m in re.finditer(r"(?=[ACGT]{24}[ACGT]{4})", search):
        rel = m.start()
        idx = rel + search_start
        pam = uSeq[idx + 24:idx + 28]
        if len(pam) != 4:
            continue
        entry = _HT_PAMDA.get(pam)
        if entry is None:
            continue
        pam_var, pam_score = entry[0], entry[1]
        if pam_variants and pam_var not in pam_variants:
            continue
        proto30 = uSeq[idx:idx + 30]
        if len(proto30) != 30:
            continue
        nick_dist = len(preHom) - (idx + 21) + 1
        start20 = idx + 4
        end20 = start20 + 20
        # Window = 4 bp pad + 20 bp protospacer + downstream, i.e. 75 nt total
        # (DOWNSTREAM_LEN = 71). Verified byte-for-byte against the OptiPrime
        # webserver, and verified end-to-end: this 75 nt window yields RTT up to
        # 33 and reproduces the webserver's top score. DO NOT LENGTHEN IT — an
        # earlier attempt to extend it to 120 truncated the RTT search instead.
        window = uSeq[start20 - 4:start20 + DOWNSTREAM_LEN + uDelta]
        if window[PS20_OFFSET:PS20_OFFSET + SPACER_LEN] != uSeq[start20:start20 + 20]:
            continue
        out.append({
            "direction": direction,
            "spacer": uSeq[start20:start20 + 20],
            "unedited": window,
            "start20": start20, "end20": end20,
            "proto30": proto30, "pam": pam, "pam_var": pam_var,
            "pam_score": pam_score, "nick_dist": nick_dist,
        })
    return out


def find_protospacers(editseq, pam_variants=DEFAULT_PAM_VARIANTS,
                      top_k_per_pamvar=TOP_K_PER_PAMVAR):
    """Enumerate protospacers for the edit in `editseq`, mirroring the OptiPrime
    webserver: search both strands via HT-PAMDA PAMs, group by PAM variant, and
    keep the nearest-nick `top_k_per_pamvar` per variant (the webserver shows all
    and lets the user pick; batch mode auto-takes the closest few).

    Returns list of dicts with: strand, spacer, unedited, edit_pos_in_window,
    ref, alt, nick_dist, pam, pam_var — laid out so unedited[4:24]==spacer.
    """
    parsed = parse_editseq(editseq)
    uSeq = editseq_unedited(parsed)
    eSeq = editseq_edited(parsed)
    edit_i = parsed["edit_index"]

    fwd = _find_protos_direction(uSeq, eSeq, "+", pam_variants)
    uSeqR, eSeqR = _revcomp(uSeq), _revcomp(eSeq)
    rev = _find_protos_direction(uSeqR, eSeqR, "-", pam_variants)

    entries = fwd + rev

    # group by PAM variant, sort by nick distance (nearest first), keep top-K
    by_var = {}
    for e in entries:
        by_var.setdefault(e["pam_var"], []).append(e)
    kept = []
    for var, lst in by_var.items():
        lst.sort(key=lambda x: x["nick_dist"])
        kept.extend(lst[:top_k_per_pamvar] if top_k_per_pamvar else lst)

    # locate the edited base inside each window, and record ref/alt on that strand
    ref_fwd = parsed["ref"]
    alt_fwd = parsed["alt"]
    results = []
    for e in kept:
        window = e["unedited"]
        if e["direction"] == "+":
            edit_in_window = edit_i - (e["start20"] - 4)
            ref_b, alt_b = ref_fwd, alt_fwd
        else:
            # On the reverse strand the window was cut from uSeqR, so both the
            # edit index and the window start are already in uSeqR coordinates.
            # The edit BLOCK occupies uSeq[edit_i : edit_i + len(ref)], and
            # reverse-complementing puts its START at len - edit_i - len(ref).
            # Using len-1-edit_i (correct only for a 1 bp edit) shifts every
            # multi-base edit by len(ref)-1 — invisible for point mutations,
            # wrong for every bystander block.
            ref_len = len(ref_fwd) if ref_fwd != "-" else 0
            edit_i_rev = len(uSeq) - edit_i - max(ref_len, 1)
            win_start_rev = (e["start20"] - 4)
            edit_in_window = edit_i_rev - win_start_rev
            ref_b, alt_b = _revcomp(ref_fwd), _revcomp(alt_fwd)
        if edit_in_window is None or edit_in_window < 0 or edit_in_window >= len(window):
            continue
        results.append({
            "strand": e["direction"],
            "spacer": e["spacer"],
            "unedited": window,
            # Full genomic context ON THE PROTOSPACER STRAND, plus where the
            # protospacer starts inside it. Stage-2 library assembly needs this:
            # the 75 nt OptiPrime window is a scoring input, not enough sequence
            # to cut a self-targeting target region from.
            "strand_seq": uSeq if e["direction"] == "+" else uSeqR,
            "proto_start": e["start20"],
            "edit_pos_in_window": edit_in_window,
            "ref": ref_b, "alt": alt_b,
            "nick_to_edit": e["nick_dist"],
            "pam": e["pam"], "pam_var": e["pam_var"],
        })
    return results



# ---------------------------------------------------------------------------
# 3. build OptiPrime run dicts (with silent-bystander edit_segments)
# ---------------------------------------------------------------------------

def build_optiprime_runs(sequence_name, editseq, codon_blocks=None,
                         settings=None, bystander_edits=None):
    """Produce OptiPrime run dicts for one editseq (one per protospacer).

    codon_blocks: [(index_in_editseq, [codon_option, ...]), ...] — the per-codon
    synonymous choices for silent bystanders, in editseq coordinates. OptiPrime
    combines them ITSELF and prunes in four rounds (see its README: "The
    edit_segments are silent edit options that should be tried combinatorially"),
    so one run covers every combination. Do NOT pass one pre-enumerated
    combination per call: that is PRIDICT2's interface, not OptiPrime's, and it
    repeats the protospacer search and the full RTT x PBS grid for each one.

    Returns a list of {name, unedited, edit_segments, settings, _meta}.
    """
    if bystander_edits is not None and codon_blocks is None:
        codon_blocks = bystander_edits          # backwards-compatible alias
    runs = []
    protos = find_protospacers(editseq)
    for k, ps in enumerate(protos):
        window = ps["unedited"]
        edit_pos = ps["edit_pos_in_window"]
        alt = ps["alt"] if ps["alt"] != "-" else ""
        win_start = ps["proto_start"] - PS20_OFFSET
        blocks = _blocks_to_window(codon_blocks, ps["strand"],
                                   len(ps["strand_seq"]), win_start, len(window))
        segments = _segments_with_bystanders(window, edit_pos, ps["ref"], alt,
                                             blocks)
        _check_segments(window, segments, edit_pos, ps["ref"], alt,
                        blocks, ps["strand"], sequence_name)
        runs.append({
            "name": f"{sequence_name}__ps{k}_{ps['strand']}",
            "unedited": window,
            "edit_segments": segments,
            "settings": settings or {},
            "_meta": {
                "sequence_name": sequence_name,
                "strand": ps["strand"],
                "spacer": ps["spacer"],
                "edit_pos_in_window": edit_pos,
                "nick_to_edit": ps["nick_to_edit"],
                "strand_seq": ps["strand_seq"],
                "proto_start": ps["proto_start"],
                "n_codon_blocks": len(blocks),
            },
        })
    return runs


def _check_segments(window, segments, edit_pos, ref, alt, bystander_edits,
                    strand, sequence_name):
    """Fail loudly if the segments don't reconstruct the window as expected.

    Two silent-corruption modes have bitten this pipeline before:
      * the trailing fixed segment going missing, which shortens the downstream
        homology and collapses OptiPrime's RTT search to the shortest designs;
      * a ref/alt mix-up, which makes `unedited` already carry the edited base.
    Both are cheap to assert here and expensive to spot in the final library.
    """
    baseline = "".join(seg[0] for seg in segments)
    if len(baseline) != len(window):
        raise ValueError(
            f"{sequence_name} [{strand}]: edit_segments rebuild to "
            f"{len(baseline)} nt but the window is {len(window)} nt. The "
            f"segments must cover the whole window.")
    if window[edit_pos:edit_pos + len(ref)] != ref:
        raise ValueError(
            f"{sequence_name} [{strand}]: window[{edit_pos}] is "
            f"'{window[edit_pos:edit_pos + len(ref)]}' but ref is '{ref}'. "
            f"`unedited` must carry the UNEDITED base — check the ref/alt "
            f"orientation for this strand.")
    # A bystander edit is a BLOCK (e.g. 13 bp) in which only some positions
    # differ, so the whole ref span is legitimate territory for changes.
    edit_span = set(range(edit_pos, edit_pos + max(len(ref), 1)))
    expected_changes = set(edit_span)
    for bi, opts in (bystander_edits or []):
        expected_changes |= set(range(bi, bi + len(opts[0])))
    actual_changes = {i for i in range(len(window)) if window[i] != baseline[i]}
    if not actual_changes <= expected_changes:
        raise ValueError(
            f"{sequence_name} [{strand}]: edited sequence differs from the "
            f"window at {sorted(actual_changes - expected_changes)}, which is "
            f"outside the intended edit positions {sorted(expected_changes)}.")
    if alt and alt != ref and not (actual_changes & edit_span):
        raise ValueError(
            f"{sequence_name} [{strand}]: the edit at {edit_pos} "
            f"({ref}->{alt}) did not make it into edit_segments.")


def _blocks_to_window(codon_blocks, strand, useq_len, win_start, win_len):
    """Map codon option blocks from editseq coordinates into window coordinates.

    On the reverse strand a block occupying editseq[i, i+w) sits at
    uSeqR[len-i-w, len-i), and its options must be reverse-complemented. Blocks
    that do not fit entirely inside the 75 nt window are dropped: an option
    cannot be offered for bases OptiPrime never sees.
    """
    out = []
    for i, opts in (codon_blocks or []):
        w = len(opts[0])
        if any(len(o) != w for o in opts):
            raise ValueError(f"codon block at {i} has options of unequal length")
        if strand == "+":
            wi, wopts = i - win_start, [o.upper() for o in opts]
        else:
            wi = (useq_len - i - w) - win_start
            wopts = [_revcomp(o.upper()) for o in opts]
        if wi < 0 or wi + w > win_len:
            continue
        out.append((wi, wopts))
    out.sort(key=lambda x: x[0])
    return out


def _segments_with_bystanders(window, edit_pos, ref, alt, bystander_edits):
    """Split `window` into OptiPrime edit_segments, applying the main edit and
    offering bystander options combinatorially.

    Positions (edit_pos + bystander indices) are on the SAME `window` coordinate
    system. We walk left-to-right, cutting the sequence at each edit position.
    """
    ref_len = len(ref) if ref != "-" else 0
    main_alt = alt if alt != "-" else ""

    # `bystander_edits` is [(index_in_window, [option, ...]), ...] — the
    # synonymous codon choices OptiPrime combines itself. The block covering the
    # main edit carries that edit in every one of its options, so it replaces the
    # standalone main-edit segment; if no block covers the edit (its codon can
    # run past the end of the 75 nt window) one is added so it is never lost.
    cuts = [(i, len(opts[0]), list(opts)) for i, opts in (bystander_edits or [])]
    covered = any(i <= edit_pos and edit_pos + max(ref_len, 1) <= i + span
                  for i, span, _ in cuts)
    if not covered:
        cuts.append((edit_pos, ref_len, [main_alt]))
    cuts.sort(key=lambda c: c[0])
    for (a, wa, _), (bb, _wb, _o) in zip(cuts, cuts[1:]):
        if a + wa > bb:
            raise ValueError(f"overlapping edit_segments at {a}(+{wa}) and {bb}")

    segments = []
    prev = 0
    for (idx, span, opts) in cuts:
        if idx > prev:
            segments.append([window[prev:idx]])     # fixed lead segment
        segments.append(list(opts))                 # variable (or single) segment
        prev = idx + span
    if prev < len(window):
        segments.append([window[prev:]])            # fixed tail
    return segments


# ---------------------------------------------------------------------------
# 4. run OptiPrime
# ---------------------------------------------------------------------------

def run_optiprime(runs, optiprime_home="~/optiprime-src", conda_env="optiprime",
                  graph_rx="graphs/pe_model.rx", weights_glob="weights/*",
                  out_dir=None, run_name=None):
    """Write the runs to JSON and call DESIGN_PE.py. Mirrors run_pridict.

    Returns the output directory path (contains one subfolder per run).
    """
    home = Path(optiprime_home).expanduser()
    out_dir = Path(out_dir).expanduser() if out_dir else \
        Path(tempfile.mkdtemp(prefix="optiprime_"))
    out_dir.mkdir(parents=True, exist_ok=True)

    # strip _meta before writing (OptiPrime ignores unknown keys, but keep clean)
    payload = [{k: v for k, v in r.items() if not k.startswith("_")} for r in runs]
    job_json = out_dir / "optiprime_input.json"
    with job_json.open("w") as f:
        json.dump(payload, f)

    cmd = (f'cd "{home}" && conda run -n {conda_env} '
           f'python DESIGN_PE.py run '
           f'--run_data "{job_json}" '
           f'--graph_rx {graph_rx} '
           f'--weight_dirs {weights_glob} '
           f'--out_path "{out_dir}"')
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("OptiPrime DESIGN_PE.py failed:\n"
                           + (proc.stderr or proc.stdout)[-3000:])
    return out_dir


# ---------------------------------------------------------------------------
# 5. translate OptiPrime output -> PRIDICT-style summary columns
# ---------------------------------------------------------------------------

# Stage-2 assembly requires exactly these columns:
def _mark_edits_lowercase(rtt, unedited_window):
    """Lowercase the bases of an RT template that differ from the genome.

    PRIDICT2 marks edited bases in RTrevcomp with lowercase, and the toolkit
    relies on that to show where the edits are (a bystander design is one with
    several lowercase bases). OptiPrime emits plain uppercase, so the marker has
    to be reconstructed. The RT template is the reverse complement of the edited
    window starting at the nick, so comparing it against the UNEDITED window at
    the same coordinates recovers exactly the same set of positions.
    """
    rtt = str(rtt).upper().replace("U", "T")
    seg = unedited_window[NICK_POS_IN_WINDOW:NICK_POS_IN_WINDOW + len(rtt)]
    if len(seg) != len(rtt):
        return rtt                      # cannot align; leave as-is
    rc_seg = _revcomp(seg)              # same orientation as the RT template
    return "".join(b.lower() if b != rc_seg[i] else b
                   for i, b in enumerate(rtt))


_REQUIRED = ["Spacer-Sequence", "RTrevcomp", "PBSrevcomp", "PBSlength",
             "RTlength", "wide_initial_target", "Original_Sequence"]


def _read_optiprime_results(run_out_dir: Path) -> pd.DataFrame:
    f = run_out_dir / "full_results.txt.gz"
    with gzip.open(f, "rt") as fh:
        return pd.read_csv(fh, sep="\t")


def optiprime_to_summary(out_dir, runs):
    """Combine every run's full_results into one PRIDICT-style summary DataFrame.

    OptiPrime columns: pegRNA_name, edit_name, RTT_len, PBS_len, spacer, RTT,
    PBS, OptiPrime_score.
    We map/derive the columns Stage-2 needs. RTT/PBS are RNA on OptiPrime's side
    (make_rtt/make_pbs return dna_to_rna(revcomp(...))), i.e. they are already the
    reverse-complement orientation the toolkit calls RTrevcomp/PBSrevcomp — we
    convert RNA (U) back to DNA (T).
    """
    out_dir = Path(out_dir)
    meta_by_name = {r["name"]: r["_meta"] for r in runs}
    frames = []
    skipped = []
    for r in runs:
        rd = out_dir / r["name"]
        if not (rd / "full_results.txt.gz").exists():
            continue
        df = _read_optiprime_results(rd)
        if df.empty:
            continue
        meta = meta_by_name[r["name"]]
        window = r["unedited"]
        edit_pos = meta["edit_pos_in_window"]

        # Stage-2 geometry. The OptiPrime window is only 75 nt, which is enough
        # to score a pegRNA but NOT enough to cut a self-targeting target region
        # from: the notebook needs ~99 nt around the protospacer plus 80-95 nt of
        # genomic sequence downstream of it. Rebuild both from the full editseq,
        # on the protospacer strand, exactly as PRIDICT2.0 reports them.
        strand_seq = meta.get("strand_seq")
        proto_start = meta.get("proto_start")
        if not strand_seq or proto_start is None:
            raise ValueError(
                f"{r['name']}: run is missing strand_seq/proto_start. Re-generate "
                f"the runs with this version of build_optiprime_runs.")
        wide_start = proto_start - PROTO_OFFSET_IN_WIDE
        wide_end = wide_start + WIDE_TARGET_LEN
        if wide_start < 0 or wide_end > len(strand_seq):
            skipped.append(
                f"{r['name']}: only {len(strand_seq)} nt of context on the "
                f"{meta['strand']} strand, not enough for a {WIDE_TARGET_LEN} nt "
                f"wide_initial_target at protospacer position {proto_start}.")
            continue
        wide_initial_target = strand_seq[wide_start:wide_end]
        genomic_proto = strand_seq[proto_start:proto_start + SPACER_LEN_WITH_G]
        # 5' G by SUBSTITUTION, never prepended -> always exactly 20 nt.
        spacer = "G" + genomic_proto[1:]

        # OptiPrime reports the same protospacer (it writes the substituted base
        # as a lowercase g, and in some cases prepends one instead). Cross-check
        # the genomic part so a coordinate slip cannot pass silently.
        op_spacer = str(df["spacer"].iloc[0]).upper()
        if op_spacer[-19:] != genomic_proto[-19:]:
            raise ValueError(
                f"{r['name']}: rebuilt protospacer {genomic_proto} disagrees with "
                f"OptiPrime's {op_spacer}. Check the strand/coordinate mapping.")

        out = pd.DataFrame()
        out["Spacer-Sequence"] = [spacer] * len(df)
        # OptiPrime RTT/PBS are RNA in revcomp orientation -> DNA revcomp
        out["RTrevcomp"] = df["RTT"].apply(
            lambda s: _mark_edits_lowercase(s, window))
        out["PBSrevcomp"] = df["PBS"].str.upper().str.replace("U", "T", regex=False)
        out["RTlength"] = df["RTT_len"].astype(int)
        out["PBSlength"] = df["PBS_len"].astype(int)
        out["OptiPrime_score"] = df["OptiPrime_score"].astype(float)
        # wide_initial_target / Original_Sequence: the toolkit uses these to
        # locate the edit and verify the construct. Use the OptiPrime window
        # (unedited) as Original_Sequence, and a wide slice around the edit as
        # wide_initial_target (mirrors PRIDICT's ~column semantics: the target
        # region the pegRNA acts on). We take the protospacer + downstream.
        # Stage-2 reads these two to recover the downstream genomic sequence and
        # size the target region (notebook cell 22/23), so both must be the FULL
        # context on the protospacer strand — not the 75 nt scoring window.
        out["Original_Sequence"] = strand_seq
        out["wide_initial_target"] = wide_initial_target
        # keep the scoring window for traceability / debugging
        out["OptiPrime_window"] = window
        out["OptiPrime_edit_pos_in_window"] = edit_pos
        # Which silent-codon combination this row is. expand_edit_options letters
        # the options of each VARIABLE segment A, B, C...; the all-A code is the
        # main edit with every codon left at its reference, i.e. no bystander.
        # 'PE' is what OptiPrime emits when no segment varies at all.
        n_var = sum(1 for seg in r["edit_segments"] if len(seg) > 1)
        ref_code = "A" * n_var if n_var else "PE"
        edit_name = (df["edit_name"].astype(str) if "edit_name" in df.columns
                     else pd.Series([ref_code] * len(df), index=df.index))
        out["edit_name"] = edit_name.values
        out["bystander"] = ["no" if e in (ref_code, "PE") else "yes"
                            for e in edit_name]
        out["n_silent_codons"] = [0 if e in (ref_code, "PE")
                                  else sum(1 for c in e if c != "A")
                                  for e in edit_name]
        # carry provenance so the assembler can name/group correctly
        out["sequence_name"] = meta["sequence_name"]
        out["_strand"] = meta["strand"]
        # OptiPrime writes full_results in RTT/PBS enumeration order, NOT by
        # score, so the first row is always the shortest design and usually one
        # of the worst. Anything downstream that takes .head(n) or .iloc[0]
        # would silently pick that. Sort best-first here.
        out = out.sort_values("OptiPrime_score", ascending=False,
                              kind="mergesort").reset_index(drop=True)
        frames.append(out)

    if not frames:
        raise RuntimeError("No OptiPrime results were produced for any run."
                           + ("\n" + "\n".join(skipped) if skipped else ""))
    summary = pd.concat(frames, ignore_index=True)
    # Global best-first ordering, so a downstream "top N per mutation" that
    # relies on row order (rather than sorting by a PRIDICT-specific score
    # column that does not exist here) still gets the best designs.
    summary = summary.sort_values("OptiPrime_score", ascending=False,
                                  kind="mergesort").reset_index(drop=True)

    missing = [c for c in _REQUIRED if c not in summary.columns]
    if missing:
        raise ValueError(f"OptiPrime->summary is missing required columns: {missing}")
    return summary


# ---------------------------------------------------------------------------
# 6. top-level entry: batch.csv (editseq column) -> PRIDICT-style summary
# ---------------------------------------------------------------------------

def run_batch_optiprime(batch_csv_text, optiprime_home="~/optiprime-src",
                        conda_env="optiprime", out_dir=None,
                        codon_blocks_map=None, settings=None,
                        bystander_map=None):
    """Full Stage-2 replacement: take the toolkit's batch.csv (with columns
    `sequence_name` and `editseq`) and produce a PRIDICT-style summary DataFrame
    by scoring every protospacer of every edit with OptiPrime.

    bystander_map: optional {sequence_name: [(idx_in_window, [alts]), ...]} to
    offer synonymous bystanders as edit_segments options. The toolkit's existing
    bystander enumeration supplies these; OptiPrime only scores them.

    Returns (summary_df, report).
    """
    import io as _io
    df = pd.read_csv(_io.StringIO(batch_csv_text))
    if "editseq" not in df.columns or "sequence_name" not in df.columns:
        raise ValueError("batch.csv needs 'sequence_name' and 'editseq' columns; "
                         f"found {list(df.columns)}")

    all_runs = []
    skipped = []
    for _, row in df.iterrows():
        name = str(row["sequence_name"])
        editseq = str(row["editseq"])
        blocks = (codon_blocks_map or bystander_map or {}).get(name)
        try:
            runs = build_optiprime_runs(name, editseq, codon_blocks=blocks,
                                        settings=settings)
        except Exception as e:
            skipped.append((name, f"run-build failed: {e}"))
            continue
        if not runs:
            skipped.append((name, "no valid protospacer found"))
            continue
        all_runs.extend(runs)

    if not all_runs:
        raise RuntimeError("No OptiPrime runs could be built from the batch "
                           "(no protospacers found for any edit).")

    out = run_optiprime(all_runs, optiprime_home=optiprime_home,
                        conda_env=conda_env, out_dir=out_dir)
    summary = optiprime_to_summary(out, all_runs)

    n_bystander_in = int(sum("_byst" in str(n) for n in df["sequence_name"]))
    n_bystander_skipped = int(sum("_byst" in str(n) for n, _ in skipped))
    reason_counts = {}
    for _n, why in skipped:
        key = str(why).split(":")[0].strip()[:60]
        reason_counts[key] = reason_counts.get(key, 0) + 1

    report = {
        "n_input_sequences": len(df),
        "n_runs": len(all_runs),
        "n_scored_pegRNAs": len(summary),
        "n_skipped": len(skipped),
        "skipped": skipped,
        # Surfaced separately so the UI can show a warning instead of a green
        # tick: silently dropping input sequences is how 400 bystander designs
        # went missing without anything in the interface changing.
        "skipped_reason_counts": reason_counts,
        "n_bystander_input": n_bystander_in,
        "n_bystander_skipped": n_bystander_skipped,
        "warning": (
            f"{len(skipped)} of {len(df)} input sequences produced no pegRNA"
            + (f" ({n_bystander_skipped} of them bystander designs)"
               if n_bystander_skipped else "")
            + ". Reasons: "
            + "; ".join(f"{k} x{v}" for k, v in
                        sorted(reason_counts.items(), key=lambda kv: -kv[1])[:3])
        ) if skipped else "",
        "out_dir": str(out),
    }
    if skipped:
        print("\n*** OptiPrime WARNING: " + report["warning"], flush=True)
        for name, why in skipped[:5]:
            print(f"      {name}: {why}", flush=True)
    return summary, report
