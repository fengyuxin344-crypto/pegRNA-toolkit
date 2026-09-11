"""
pegRNA Library Designer — Core Engine
=====================================
Turns a mutation list (for ANY gene) into PRIDICT-ready input sequences,
and later assembles PRIDICT predictions into a self-targeting library.

This module contains ONLY logic (no UI). It is adapted from the original
CSF1R ALSP library-design notebook, generalised to any gene and hardened
against the fragile mutation-parsing in the original.

Two stages:
  Stage 1  build_pridict_input()  — mutation list  -> batch CSV for PRIDICT2.0
  Stage 2  build_library()        — PRIDICT summary -> self-targeting library
"""

import itertools
import ast
import datetime
import glob
import os
import random
import re
import shutil
import subprocess
import time
import requests
import pandas as pd


# ---------------------------------------------------------------------------
# Stage 1a — Fetch gene + build CDS<->genomic coordinate map (any gene)
# ---------------------------------------------------------------------------

ENSEMBL = "https://rest.ensembl.org"


def _get(url, content_type="application/json"):
    r = requests.get(url, headers={"Content-Type": content_type}, timeout=30)
    if not r.ok:
        r.raise_for_status()
    return r.json() if content_type == "application/json" else r.text


def fetch_gene(gene_id):
    """Look up a gene by Ensembl ID (e.g. ENSG00000182578) with transcripts+exons."""
    return _get(f"{ENSEMBL}/lookup/id/{gene_id}?expand=1")


def resolve_gene_id(gene_query):
    """Accept either an Ensembl ID or a gene symbol (e.g. 'CSF1R') and return the ENSG id.

    Symbols are resolved via the human symbol lookup. For non-human species this
    would need a species parameter; kept to human here to match the original work.
    """
    gene_query = gene_query.strip()
    if re.match(r"^ENSG\d+", gene_query, re.IGNORECASE):
        return gene_query
    data = _get(f"{ENSEMBL}/xrefs/symbol/homo_sapiens/{gene_query}?object_type=gene")
    if not data:
        raise ValueError(f"Could not resolve gene symbol '{gene_query}' to an Ensembl gene ID.")
    return data[0]["id"]


def build_coordinate_map(gene_id):
    """Return (gene_info, cds_positions, cds_absolute_positions).

    cds_positions[i] is the position in CDS coordinates (codon-numbered, 1-based
    from the ATG) and cds_absolute_positions[i] is the matching position within
    the fetched genomic sequence. Handles + and - strand genes.
    (Directly generalised from the notebook's main().)
    """
    gene = fetch_gene(gene_id)
    gene_seq = _get(f"{ENSEMBL}/sequence/id/{gene_id}?type=genomic", "text/plain")
    canonical = next(t for t in gene["Transcript"] if t.get("is_canonical") == 1)
    tid = canonical["id"]
    gene_len = gene["end"] - gene["start"] + 1
    strand = gene["strand"]

    # exon relative positions (strand-aware — this is the +/- coordinate flip)
    exon_infos = []
    for exon in canonical["Exon"]:
        if strand == 1:
            rel_start = exon["start"] - gene["start"] + 1
            rel_end = exon["end"] - gene["start"] + 1
        else:
            rel_start = gene_len - (exon["end"] - gene["start"])
            rel_end = gene_len - (exon["start"] - gene["start"])
        exon_infos.append({"exon_id": exon["id"], "rel_start": rel_start, "rel_end": rel_end})
    exon_infos.sort(key=lambda e: e["rel_start"])

    cds_positions, cds_abs = [], []
    exon_start_pos = 1
    start_position = 0
    for idx, ex in enumerate(exon_infos):
        exon_seq = _get(f"{ENSEMBL}/sequence/id/{ex['exon_id']}?type=genomic", "text/plain")
        if idx == 0:
            start_position = exon_seq.find("ATG")
            if start_position == -1:
                raise ValueError("Start codon (ATG) not found in first exon.")
        exon_len = ex["rel_end"] + 1 - ex["rel_start"]
        rel_exon = list(range(exon_start_pos - start_position,
                              exon_start_pos + exon_len - start_position))
        cds_positions += rel_exon
        exon_start_pos += exon_len
        rel_abs = list(range(ex["rel_start"], ex["rel_end"] + 1))
        cds_abs += rel_abs
        if len(rel_abs) != len(rel_exon):
            raise ValueError("Coordinate map length mismatch (relative vs absolute).")

    info = {"gene_id": gene_id, "gene_name": gene["display_name"],
            "strand": strand, "gene_sequence": gene_seq, "transcript_id": tid}
    return info, cds_positions, cds_abs


# ---------------------------------------------------------------------------
# Stage 1b — Parse a mutation string (hardened vs. the original notebook)
# ---------------------------------------------------------------------------

class MutationParseError(Exception):
    pass


def parse_mutation(gdna):
    """Parse an HGVS-style coding mutation string.

    Returns dict(edit_type, location, original_base, mutated_base).
    Raises MutationParseError for anything it cannot handle CONFIDENTLY —
    the caller is expected to report these rather than silently skip them.
    """
    s = str(gdna).strip()

    # substitution: c.1234A>G
    m = re.match(r"^c\.(\d+)([ACGT])>([ACGT])$", s)
    if m:
        return {"edit_type": "substitution", "location": int(m.group(1)),
                "original_base": m.group(2), "mutated_base": m.group(3)}

    # single-base deletion: c.1234delA
    m = re.match(r"^c\.(\d+)del([ACGT])$", s)
    if m:
        return {"edit_type": "deletion", "location": int(m.group(1)),
                "original_base": m.group(2), "mutated_base": None}

    # single-base insertion: c.1234_1235insA  (kept but flagged: notebook didn't handle)
    m = re.match(r"^c\.(\d+)_(\d+)ins([ACGT]+)$", s)
    if m:
        return {"edit_type": "insertion", "location": int(m.group(1)),
                "original_base": None, "mutated_base": m.group(3)}

    raise MutationParseError(
        f"Unsupported mutation format: '{s}'. "
        "Supported: c.NNNX>Y, c.NNNdelX, c.NNN_MMMinsX. "
        "Splice-site, range deletions and complex indels need manual curation."
    )


# ---------------------------------------------------------------------------
# Stage 1b-protein — Convert protein-level (p.) mutations to c. notation
# ---------------------------------------------------------------------------
#
# WHY THIS IS DANGEROUS AND HOW WE GUARD IT
# A protein change like L483P maps to DNA ambiguously (many codons per amino
# acid) AND depends entirely on the reference transcript. This module therefore
# NEVER guesses: it reads the ACTUAL CDS for the gene's Ensembl canonical
# transcript, confirms the reference amino acid matches at that position, and
# reports every uncertainty instead of silently emitting a c. string. A wrong
# c. string is worse than no answer, because everything downstream inherits it.

# standard codon table
_CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

_AA3TO1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
    'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
    'Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
}

# Cross-check anchors: known ClinVar c. notations for well-characterised variants,
# keyed by (gene_symbol_upper, one_letter_protein_change). Used to VERIFY the
# computed result, never to replace the computation. Extend as needed.
_KNOWN_CDNA = {
    ("GBA1", "L483P"): "c.1448T>C",
    ("GBA1", "N409S"): "c.1226A>G",
    ("MAPT", "R406W"): "c.1216C>T",
}


class ProteinConversionError(Exception):
    pass


def _parse_protein_change(p):
    """Parse 'L483P', 'p.Leu483Pro', 'Leu483Pro' -> (orig1, pos, new1)."""
    s = str(p).strip()
    s = re.sub(r"^p\.", "", s)
    # three-letter form
    m = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", s)
    if m and m.group(1) in _AA3TO1 and m.group(3) in _AA3TO1:
        return _AA3TO1[m.group(1)], int(m.group(2)), _AA3TO1[m.group(3)]
    # one-letter form
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", s)
    if m:
        return m.group(1), int(m.group(2)), m.group(3)
    raise ProteinConversionError(f"Could not parse protein change '{p}'.")


def _fetch_cds(gene_id):
    """Fetch the CDS (coding sequence) of the gene's canonical transcript."""
    gene = fetch_gene(gene_id)
    canonical = next(t for t in gene["Transcript"] if t.get("is_canonical") == 1)
    tid = canonical["id"]
    cds = _get(f"{ENSEMBL}/sequence/id/{tid}?type=cds", "text/plain")
    return tid, cds.strip().upper()


def protein_to_cdna(protein_change, gene_query, gene_id=None, cds=None,
                    transcript_id=None):
    """Convert a protein change (e.g. 'L483P') to c. notation for `gene_query`.

    Reads the canonical CDS from Ensembl, verifies the reference amino acid,
    and returns a dict with the computed c. string plus verification metadata.
    NEVER returns a c. string without confirming the original amino acid matches.

    Returns dict(cdna, verified, warnings, ...). Raises ProteinConversionError
    on hard failures (unparseable, reference AA mismatch, ambiguous, etc.).
    """
    orig_aa, pos, new_aa = _parse_protein_change(protein_change)

    if gene_id is None:
        gene_id = resolve_gene_id(gene_query)
    if cds is None:
        transcript_id, cds = _fetch_cds(gene_id)

    warnings = []

    # locate the codon
    codon_start = (pos - 1) * 3
    if codon_start + 3 > len(cds):
        raise ProteinConversionError(
            f"Residue {pos} is beyond the CDS length ({len(cds)//3} aa) for "
            f"{gene_query}. Wrong transcript?")
    codon = cds[codon_start:codon_start + 3]
    translated = _CODON_TABLE.get(codon, "?")

    # HARD CHECK 1: reference amino acid must match
    if translated != orig_aa:
        raise ProteinConversionError(
            f"Reference amino acid mismatch for {gene_query} {protein_change}: "
            f"canonical transcript {transcript_id} has '{translated}' at residue "
            f"{pos}, not '{orig_aa}'. This almost always means the mutation's "
            f"numbering is based on a DIFFERENT transcript/isoform than Ensembl "
            f"canonical. Do NOT use this result — supply the matching transcript.")

    # find single-nucleotide changes within the codon that yield new_aa
    candidates = []
    for i in range(3):
        for base in "ACGT":
            if base == codon[i]:
                continue
            new_codon = codon[:i] + base + codon[i + 1:]
            if _CODON_TABLE.get(new_codon) == new_aa:
                cpos = codon_start + i + 1  # 1-based CDS position
                candidates.append((cpos, codon[i], base))

    if not candidates:
        raise ProteinConversionError(
            f"{gene_query} {protein_change}: no SINGLE-nucleotide change in codon "
            f"'{codon}' produces '{new_aa}'. This substitution needs >1 nt change "
            f"(prime editing can still do it, but it's not a simple 1bp edit — "
            f"handle manually).")

    if len(candidates) > 1:
        warnings.append(
            f"{len(candidates)} different single-nt changes in codon '{codon}' give "
            f"'{new_aa}': {['c.'+str(c[0])+c[1]+'>'+c[2] for c in candidates]}. "
            f"Pick the one matching the actual reported variant.")

    cpos, ref_base, alt_base = candidates[0]
    cdna = f"c.{cpos}{ref_base}>{alt_base}"

    # SOFT CHECK 2: cross-check against known ClinVar anchors
    key = (gene_query.upper(), f"{orig_aa}{pos}{new_aa}")
    known = _KNOWN_CDNA.get(key)
    verified = None
    if known:
        verified = (cdna == known)
        if not verified:
            warnings.append(
                f"Computed {cdna} disagrees with known ClinVar {known} for "
                f"{gene_query} {protein_change}. Using the KNOWN value; investigate "
                f"the transcript.")
            cdna = known  # trust the curated database over computation

    return {
        "protein_change": protein_change, "cdna": cdna, "codon": codon,
        "transcript_id": transcript_id, "n_candidates": len(candidates),
        "cross_checked_against_clinvar": known is not None,
        "matches_clinvar": verified, "warnings": warnings,
    }


def normalize_to_cdna(mutation, gene_query, gene_id=None, cds=None, transcript_id=None):
    """Return a c. mutation string from either a c. string or a p. change.

    If `mutation` already looks like c. notation, return it unchanged. Otherwise
    treat it as a protein change and convert. Returns (cdna_str, info_dict).
    """
    s = str(mutation).strip()
    if s.startswith("c."):
        return s, {"converted": False}
    info = protein_to_cdna(s, gene_query, gene_id=gene_id, cds=cds,
                           transcript_id=transcript_id)
    return info["cdna"], {"converted": True, **info}


# ---------------------------------------------------------------------------
# Stage 1b-variant — Unambiguous variant input classifier + resolvers
# ---------------------------------------------------------------------------
#
# The transcript-ambiguity fix. A variant written as a protein change (L483P)
# or bare c. (c.1448T>C) only means something RELATIVE TO A TRANSCRIPT, and the
# same genomic position numbers differently across transcripts (the MAPT trap:
# 441aa Tau-F vs the 758aa canonical). So we:
#   * recognise "self-coordinate" inputs that pin the coordinate system on their
#     own — genomic HGVS (g.), rsID, and versioned transcript c. (NM_/ENST:c.) —
#     and resolve them straight to a genomic window, no transcript guessing;
#   * REFUSE bare c. / protein changes that arrive WITHOUT a transcript, instead
#     of silently falling back to canonical (which is how you get wrong answers).
# Every resolved row records how it was resolved, for traceability in the report.

class VariantInputError(Exception):
    pass


_NC_CHROM = {f"NC_0000{n:02d}": str(n) for n in range(1, 23)}
_NC_CHROM.update({"NC_000023": "X", "NC_000024": "Y", "NC_012920": "MT"})


def _norm_chrom(c):
    """Normalise a chromosome token to Ensembl style (1..22, X, Y, MT)."""
    c = str(c).strip()
    base = c.split(".")[0].upper()
    if base in _NC_CHROM:
        return _NC_CHROM[base]
    return re.sub(r"^chr", "", c, flags=re.I).upper().replace("MT", "MT")


def _parse_genomic(s):
    """Recognise genomic-coordinate inputs. Returns dict or None.
    Accepts:  NC_000001.11:g.155235252A>G | chr1:g.155235252A>G | 1:g.155235252A>G
              chr1:155235252 A>G | chr1:155235252:A:G | chr1-155235252-A-G
    """
    m = re.match(r"^(?:chr)?([\w.]+):g\.(\d+)([ACGT])>([ACGT])$", s, re.I)
    if m:
        return {"chrom": _norm_chrom(m.group(1)), "pos": int(m.group(2)),
                "ref": m.group(3).upper(), "alt": m.group(4).upper()}
    m = re.match(r"^(?:chr)?([\w.]+)[:\-\s](\d+)[:\-\s]([ACGT])[>:\-\s]?([ACGT])$", s, re.I)
    if m:
        return {"chrom": _norm_chrom(m.group(1)), "pos": int(m.group(2)),
                "ref": m.group(3).upper(), "alt": m.group(4).upper()}
    return None


# ---------------------------------------------------------------------------
# Stage 1b-lookup — ClinVar variant lookup ("literature name" -> coordinates)
# ---------------------------------------------------------------------------
#
# THE BOTTLENECK THIS SOLVES
# Papers write "MAPT P301L". That is a LABEL, not a coordinate: it omits the
# transcript, the codon change, and the genome build. Turning it into something
# a pipeline can use has, until now, been a manual database lookup per variant
# (the reference ALSP notebook did exactly this by hand — see its "CHECKED"
# curated CSV). This module automates the LOOKUP but NOT the DECISION: it
# returns every ClinVar candidate and the caller/user picks. Never auto-select,
# because a protein change can map to several distinct variants (P->L can be
# CCG>CTG or CCG>CTA) and silently choosing one is how libraries go wrong.

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Genes whose community numbering differs from the current HGVS numbering.
# Searching both spellings is what makes "GBA1 N370S" (legacy) and "N409S"
# (HGVS, +39 aa signal peptide) resolve to the same variant.
_LEGACY_ALIASES = {
    "GBA1": {"offset": 39, "note": "GBA1 legacy numbering omits the 39 aa signal peptide"},
    "GBA": {"offset": 39, "note": "GBA1 legacy numbering omits the 39 aa signal peptide"},
}


def _eutils_get(path, params, api_key=None, retries=3):
    """Call E-utilities with polite rate limiting and retry/backoff.
    NCBI allows 3 req/s without a key, 10 with one."""
    params = dict(params)
    params.setdefault("retmode", "json")
    if api_key:
        params["api_key"] = api_key
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(f"{EUTILS}/{path}", params=params, timeout=30)
            if r.status_code == 429:  # rate limited
                time.sleep(1.0 * (attempt + 1))
                continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            last = e
            time.sleep(0.5 * (attempt + 1))
    raise VariantInputError(f"ClinVar/E-utilities request failed: {last}")


def _protein_query_forms(gene, protein_change):
    """Build the ClinVar search terms to try for a protein change, including the
    legacy-numbering alias for genes like GBA1."""
    forms = [f"{gene} {protein_change}"]
    try:
        orig, pos, new = _parse_protein_change(protein_change)
    except ProteinConversionError:
        return forms
    alias = _LEGACY_ALIASES.get(gene.upper())
    if alias:
        # try BOTH the legacy and the HGVS-numbered spelling
        for shifted in (pos + alias["offset"], pos - alias["offset"]):
            if shifted > 0:
                forms.append(f"{gene} {orig}{shifted}{new}")
    return forms


def _protein_from_title(title):
    """Pull the protein change out of a ClinVar title like
    'NM_005910.6(MAPT):c.902C>T (p.Pro301Leu)' -> ('P', 301, 'L').
    Returns None if the title has no usable p. change (e.g. non-coding)."""
    m = re.search(r"\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})\)", str(title))
    if not m:
        return None
    o, pos, n = m.group(1).capitalize(), int(m.group(2)), m.group(3).capitalize()
    if o not in _AA3TO1 or n not in _AA3TO1:
        return None
    return (_AA3TO1[o], pos, _AA3TO1[n])


def _transcript_from_title(title):
    """Pull the transcript accession out of a ClinVar title."""
    m = re.match(r"\s*((?:NM_|NR_|ENST)\d+(?:\.\d+)?)", str(title))
    return m.group(1) if m else ""


def _cdot_from_title(title):
    """Pull the c. change out of a ClinVar title, e.g. 'c.-146C>T' or 'c.902C>T'."""
    m = re.search(r"(c\.-?\*?\d+[ACGT]>[ACGT])", str(title))
    return m.group(1) if m else ""


# a non-coding / promoter change written the way papers write it: -146C>T
_NONCODING_RE = re.compile(r"^(?:c\.)?(-\d+)([ACGT])>([ACGT])$", re.I)


def _normalise_gene(g):
    """'TERT promoter' -> 'TERT'. Strips descriptive suffixes so what goes to
    ClinVar is the actual gene symbol."""
    g = str(g).strip()
    g = re.sub(r"\s*[-_]?\s*(promoter|gene|locus|region)\s*$", "", g, flags=re.I)
    return g.strip()


def parse_lookup_line(line):
    """Parse ONE lookup line as it is really written in a lab mutation list:

        TERT promoter, -146C>T
        GBA1, NP_000148.2, L483P, (previously refer as L444P)
        GBA1, NP_000148.2, N409S (previously refer as N370S)
        MAPT, NP_005901.2, P301L
        ALDH1L2,NP_001029345.2, I141N
        MAPT, P301L
        GBA1, chr1:g.155235843T>C            (already unambiguous -> passthrough)

    Returns dict: gene, kind ('protein'|'noncoding'|'passthrough'), change,
    cdot (for non-coding), accession, aliases (mined from the parentheses),
    warnings.
    """
    raw = str(line).strip()
    if not raw:
        raise VariantInputError("empty line")

    # 1) lift out parenthetical notes; mine them for alias spellings (L444P etc.)
    notes = re.findall(r"\(([^)]*)\)", raw)
    core = re.sub(r"\([^)]*\)", " ", raw)
    aliases = []
    for n in notes:
        for m in re.finditer(r"\b([A-Z])(\d+)([A-Z])\b", n):
            aliases.append(m.group(0))

    parts = [p.strip() for p in core.split(",") if p.strip()]
    if not parts:
        raise VariantInputError(f"'{raw}': nothing left after removing notes")

    gene = _normalise_gene(parts[0])
    if not gene:
        raise VariantInputError(f"'{raw}': no gene symbol")

    accession, change, warnings = None, None, []
    for p in parts[1:]:
        if re.match(r"^(NM_|NR_|XM_|XR_|ENST|NP_|XP_|ENSP)", p, re.I):
            accession = p
        elif change is None:
            change = p
    if change is None:
        raise VariantInputError(
            f"'{raw}': no variant found. Write it as 'GENE, CHANGE' — e.g. "
            "'MAPT, P301L' or 'TERT promoter, -146C>T'.")

    # NP_/XP_ is a PROTEIN accession, not a transcript: it cannot carry a c.
    # position, so it is not used for the query. Kept only as a cross-check.
    if accession and re.match(r"^(NP_|XP_|ENSP)", accession, re.I):
        warnings.append(
            f"'{accession}' is a protein accession, not a transcript — not used "
            "for the query. Check the transcript shown on the ClinVar hit.")

    m = _NONCODING_RE.match(change)
    if m:
        return {"gene": gene, "kind": "noncoding", "change": change,
                "cdot": f"c.{m.group(1)}{m.group(2).upper()}>{m.group(3).upper()}",
                "accession": accession, "aliases": aliases,
                "warnings": warnings, "raw": raw}

    try:
        _parse_protein_change(change)
        return {"gene": gene, "kind": "protein", "change": change, "cdot": "",
                "accession": accession, "aliases": aliases,
                "warnings": warnings, "raw": raw}
    except ProteinConversionError:
        pass

    if classify_variant(change)["kind"] in ("genomic", "rsid", "transcript_c"):
        return {"gene": gene, "kind": "passthrough", "change": change, "cdot": "",
                "accession": accession, "aliases": aliases,
                "warnings": warnings, "raw": raw}

    raise VariantInputError(
        f"'{raw}': can't tell what '{change}' is. Use a protein change (L483P), "
        "a promoter change (-146C>T), a genomic coordinate, or an rsID.")


def _from_spdi(spdi):
    """Parse ClinVar's canonical_spdi, e.g. 'NC_000017.11:46010388:C:T'.
    SPDI positions are 0-BASED, so we add 1 to get the HGVS/1-based position.
    This is the reliable location field — variation_loc often omits ref/alt."""
    parts = str(spdi or "").split(":")
    if len(parts) != 4:
        return None
    acc, pos, ref, alt = parts
    chrom = _NC_CHROM.get(acc.split(".")[0].upper())
    if not chrom:
        return None
    ref, alt = ref.upper(), alt.upper()
    if not (len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT"):
        return None  # not a clean SNV — the genomic path only takes SNVs
    try:
        return chrom, int(pos) + 1, ref, alt
    except ValueError:
        return None


def _summary_to_candidates(docsum, gene=None, resolve=True):
    """Pull candidates out of one ClinVar esummary record.

    Getting the genomic coordinate is done with a FALLBACK CHAIN, because the
    esummary schema is inconsistent about where the location lives:
      1. canonical_spdi          (0-based; carries ref/alt)
      2. variation_loc GRCh38    (only when it actually has ref/alt)
      3. rsID  -> Ensembl variant_recoder
      4. title's NM_:c. -> Ensembl variant_recoder
    Steps 3-4 reuse the recoder path that is already proven to work in this tool.
    Returns (candidates, drops) — drops explains anything that couldn't resolve.
    """
    out, drops = [], []
    title = docsum.get("title", "")
    sig = ((docsum.get("germline_classification") or {}).get("description")
           or (docsum.get("clinical_significance") or {}).get("description") or "")

    for v in (docsum.get("variation_set") or []):
        vtitle = v.get("variation_name") or title
        rsid = ""
        for x in (v.get("variation_xrefs") or []):
            if str(x.get("db_source", "")).lower() == "dbsnp":
                rsid = "rs" + str(x.get("db_id"))
                break

        locs, how = [], ""
        # 1) canonical SPDI
        spdi = _from_spdi(v.get("canonical_spdi"))
        if spdi:
            locs, how = [spdi], "spdi"
        # 2) GRCh38 variation_loc
        if not locs:
            for loc in (v.get("variation_loc") or []):
                if str(loc.get("assembly_name", "")).upper() != "GRCH38":
                    continue
                chrom, start = loc.get("chr"), loc.get("start")
                ref = (loc.get("ref") or "").upper()
                alt = (loc.get("alt") or "").upper()
                if (chrom and start and len(ref) == 1 and len(alt) == 1
                        and ref in "ACGT" and alt in "ACGT"):
                    locs, how = [(str(chrom), int(start), ref, alt)], "variation_loc"
                    break
        # 3) rsID via Ensembl recoder     4) title's NM_:c. via Ensembl recoder
        if not locs and resolve:
            for src, q in (("rsid", rsid),
                           ("title_cdot", (f"{_transcript_from_title(vtitle)}:"
                                           f"{_cdot_from_title(vtitle)}")
                            if _transcript_from_title(vtitle) and _cdot_from_title(vtitle)
                            else "")):
                if not q or q.endswith(":"):
                    continue
                try:
                    locs, how = [_recoder_genomic(q)], src
                    break
                except Exception as e:
                    drops.append(f"{vtitle}: {src} '{q}' failed ({type(e).__name__})")

        if not locs:
            drops.append(f"{vtitle}: no usable GRCh38 SNV location "
                         f"(spdi={v.get('canonical_spdi')!r}, rsid={rsid!r})")
            continue

        for chrom, pos, ref, alt in locs:
            out.append({
                "gene": gene,
                "genomic": f"chr{chrom}:g.{pos}{ref}>{alt}",
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "rsid": rsid,
                "clinvar_id": docsum.get("uid", ""),
                "significance": sig,
                "title": vtitle,
                "resolved_by": how,
                "protein": _protein_from_title(vtitle),
                "protein_field": str(docsum.get("protein_change") or "").strip(),
                "cdot": _cdot_from_title(vtitle),
                "transcript": _transcript_from_title(vtitle),
            })
    return out, drops


def _protein_from_field(s):
    """Parse ClinVar's own one-letter protein_change field, e.g. 'P301L'.
    That field can list several changes; take the first parseable one."""
    for tok in re.split(r"[,\s;]+", str(s or "")):
        try:
            return _parse_protein_change(tok)
        except ProteinConversionError:
            continue
    return None


def _protein_matches(cand_protein, want, gene):
    """How well does a candidate's protein change match what was asked for?

    Returns:
      'exact'    same amino acids AND same position (or the known legacy offset)
      'aa_only'  same amino-acid change, DIFFERENT position — this is normal:
                 ClinVar titles a record on ITS chosen transcript, which may
                 number differently from the one the literature uses (MAPT is
                 the classic case: 441 aa Tau-F vs the longer MANE isoform).
                 KEPT, but flagged so the user verifies.
      False      different amino acids (P301L vs P301Q) — objectively not the
                 requested variant, so it is dropped.
      None       can't tell — kept and flagged unverified.
    """
    if cand_protein is None or want is None:
        return None
    c_o, c_pos, c_n = cand_protein
    w_o, w_pos, w_n = want
    if (c_o, c_n) != (w_o, w_n):
        return False                      # different amino-acid change: drop
    if c_pos == w_pos:
        return "exact"
    alias = _LEGACY_ALIASES.get(str(gene).upper())
    if alias and abs(c_pos - w_pos) == alias["offset"]:
        return "exact"                    # legacy vs HGVS numbering, same variant
    return "aa_only"                      # transcript numbering differs — verify


def _query_forms(spec):
    """Every ClinVar search term worth trying for one parsed lookup line.
    Several spellings are tried because ClinVar's free-text search is picky:
    a bare 'MAPT P301L' can miss where 'MAPT[gene] AND P301L' hits."""
    gene, kind, change = spec["gene"], spec["kind"], spec["change"]
    forms = []

    if kind == "noncoding":
        cd = spec["cdot"]                       # c.-146C>T
        forms += [f"{gene}[gene] AND {cd}", f"{gene} {cd}", f"{gene}[gene] AND {change}"]
        return forms

    # protein change: try the plain, the field-qualified, and the 3-letter form
    spellings = [change] + list(spec.get("aliases") or [])
    try:
        o, pos, n = _parse_protein_change(change)
        alias = _LEGACY_ALIASES.get(gene.upper())
        if alias:                                # GBA1 legacy numbering (+/- 39 aa)
            for shifted in (pos + alias["offset"], pos - alias["offset"]):
                if shifted > 0:
                    spellings.append(f"{o}{shifted}{n}")
        one2three = {v: k for k, v in _AA3TO1.items()}
        if o in one2three and n in one2three:
            spellings.append(f"p.{one2three[o]}{pos}{one2three[n]}")
    except ProteinConversionError:
        pass

    seen = set()
    for s in spellings:
        if s in seen:
            continue
        seen.add(s)
        forms += [f"{gene}[gene] AND {s}", f"{gene} {s}"]
    return forms


def lookup_variants(spec, api_key=None, retmax=20, strict=True):
    """Look one parsed lookup line (from parse_lookup_line) up in ClinVar.

    Returns (candidates, diag). Each candidate carries an unambiguous GRCh38
    coordinate plus rsID, significance, and the ClinVar title (which shows the
    transcript and c. that ClinVar used, so you can verify).

    When `strict`, candidates that objectively do NOT match are dropped: if you
    asked for P301L, a P301Q record at the same position is a different variant.
    Only genuinely equivalent candidates survive — those are real ambiguity and
    the USER picks. `diag` reports what each query term returned, so a zero
    result can be debugged instead of just saying "not found".
    """
    want = None
    if spec["kind"] == "protein":
        try:
            want = _parse_protein_change(spec["change"])
        except ProteinConversionError:
            pass
    want_cdot = spec.get("cdot", "").replace(" ", "").upper()

    seen, candidates, diag = set(), [], []
    for term in _query_forms(spec):
        es = _eutils_get("esearch.fcgi",
                         {"db": "clinvar", "term": term, "retmax": retmax},
                         api_key=api_key)
        ids = ((es.get("esearchresult") or {}).get("idlist") or [])
        step = {"term": term, "n_ids": len(ids), "n_docs": 0,
                "n_raw": 0, "n_kept": 0, "drops": [], "filtered": 0}
        if not ids:
            diag.append(step)
            time.sleep(0.34)
            continue
        su = _eutils_get("esummary.fcgi",
                         {"db": "clinvar", "id": ",".join(ids)}, api_key=api_key)
        result = su.get("result") or {}
        uids = result.get("uids") or [k for k in result.keys() if k != "uids"]
        step["n_docs"] = len(uids)
        for uid in uids:
            cands, drops = _summary_to_candidates(result.get(uid, {}), gene=spec["gene"])
            step["n_raw"] += len(cands)
            step["drops"] += drops[:2]
            for cand in cands:
                key = (cand["chrom"], cand["pos"], cand["ref"], cand["alt"])
                if key in seen:
                    continue
                if spec["kind"] == "noncoding":
                    got = (cand.get("cdot") or "").replace(" ", "").upper()
                    match = (got == want_cdot) if got else None
                else:
                    # check BOTH protein sources: the title (transcript-dependent
                    # numbering) and ClinVar's own protein_change field. Take the
                    # best verdict — only drop when we are sure the amino-acid
                    # change itself is different.
                    verdicts = []
                    for src in (cand.get("protein"),
                                _protein_from_field(cand.get("protein_field"))):
                        verdicts.append(_protein_matches(src, want, spec["gene"]))
                    if "exact" in verdicts:
                        match = "exact"
                    elif "aa_only" in verdicts:
                        match = "aa_only"
                    elif any(v is None for v in verdicts):
                        match = None
                    else:
                        match = False
                if strict and match is False:
                    step["filtered"] += 1
                    continue
                cand["protein_match"] = match
                cand["query"] = term
                seen.add(key)
                candidates.append(cand)
                step["n_kept"] += 1
        diag.append(step)
        time.sleep(0.34)          # stay under NCBI's 3 req/s no-key limit
        if candidates:
            break                 # a term hit — no need to try more spellings

    _rank = {"exact": 0, "aa_only": 1, None: 2}
    candidates.sort(key=lambda c: (_rank.get(c.get("protein_match"), 3),
                                   "athogenic" not in (c.get("significance") or "")))
    return candidates, diag


def _ensembl_fallback(spec):
    """When ClinVar has NO record (novel / rare variants like ALDH1L2 I141N),
    resolve the protein change against the gene's Ensembl transcript CDS instead.
    protein_to_cdna verifies the REFERENCE amino acid, so a wrong transcript
    raises rather than silently producing a wrong coordinate.
    Returns a one-element candidate list, or [] if it can't be done."""
    if spec["kind"] != "protein":
        return []
    gid = resolve_gene_id(spec["gene"])
    tid, cds = _fetch_cds(gid)
    conv = protein_to_cdna(spec["change"], spec["gene"], gene_id=gid,
                           cds=cds, transcript_id=tid)
    chrom, pos, ref, alt = _recoder_genomic(f"{tid}:{conv['cdna']}")
    warns = list(conv.get("warnings") or [])
    warns.append(f"Not in ClinVar — resolved against Ensembl transcript {tid}. "
                 "Reference amino acid was verified, but confirm this is the "
                 "transcript your literature numbers on.")
    return [{
        "gene": spec["gene"],
        "genomic": f"chr{chrom}:g.{pos}{ref}>{alt}",
        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
        "rsid": "", "clinvar_id": "", "significance": "not in ClinVar",
        "title": f"{tid}:{conv['cdna']} ({spec['change']}) — via Ensembl CDS",
        "resolved_by": "ensembl_cds",
        "protein": None, "protein_field": "", "cdot": conv["cdna"],
        "transcript": tid,
        "protein_match": "exact",
        "notes": warns,
    }]


def lookup_variant_table(lines, api_key=None):
    """Look up many free-form lines. Returns (results, errors).
    Each result: {gene, query, kind, candidates, warnings, diag}."""
    results, errors = [], []
    for line in lines:
        line = str(line).strip()
        if not line:
            continue
        try:
            spec = parse_lookup_line(line)
        except VariantInputError as e:
            errors.append((line, str(e)))
            continue

        if spec["kind"] == "passthrough":
            # already unambiguous (genomic / rsID / accession:c.) — no lookup needed
            results.append({
                "gene": spec["gene"], "query": spec["change"], "kind": "passthrough",
                "warnings": spec["warnings"], "diag": [],
                "candidates": [{"genomic": spec["change"], "title":
                                "already unambiguous — used as given",
                                "rsid": "", "significance": "", "protein_match": True}],
            })
            continue

        try:
            cands, diag = lookup_variants(spec, api_key=api_key)
        except VariantInputError as e:
            errors.append((line, f"ClinVar request failed: {e}"))
            continue

        if not cands:
            # ClinVar has nothing usable — try resolving against the Ensembl CDS
            fb_err = ""
            try:
                cands = _ensembl_fallback(spec)
                if cands:
                    spec["warnings"] = list(spec["warnings"]) + list(cands[0].get("notes") or [])
            except Exception as e:
                fb_err = f"{type(e).__name__}: {e}"

        if not cands:
            bits = []
            for d in diag[:4]:
                bits.append(f"{d['term']}: {d['n_ids']} id / {d.get('n_docs',0)} doc / "
                            f"{d.get('n_raw',0)} parsed / {d.get('filtered',0)} filtered out")
            why = [x for d in diag for x in d.get("drops", [])][:2]
            msg = "no usable ClinVar candidate. " + " | ".join(bits)
            if why:
                msg += " || " + " ; ".join(why)
            if fb_err:
                msg += " || Ensembl fallback also failed: " + fb_err
            msg += " -- supply a genomic coordinate or rsID directly instead."
            errors.append((line, msg))
            continue

        results.append({"gene": spec["gene"], "query": spec["change"],
                        "kind": spec["kind"], "candidates": cands,
                        "warnings": spec["warnings"], "diag": diag})
    return results, errors


def classify_variant(raw, transcript=None):
    """Classify one mutation input. Returns a dict with at least `kind` in
    {genomic, rsid, transcript_c, bare_c, protein, unknown}, the parsed pieces,
    the transcript if one was supplied, and `needs_transcript` (True when the
    input is transcript-dependent but no transcript was given)."""
    s = str(raw).strip()
    tx = str(transcript).strip() if transcript not in (None, "") else None
    if tx and tx.lower() in ("nan", "none"):
        tx = None

    # inline accession:c.  (NM_000157.4:c.1448T>C  /  ENST00000368373.8:c.1226A>G)
    m = re.match(r"^\s*((?:NM_|NR_|XM_|XR_|ENST)\d+(?:\.\d+)?)\s*:\s*(c\..+)$", s, re.I)
    if m:
        return {"kind": "transcript_c", "transcript": m.group(1),
                "cdot": m.group(2).strip(), "raw": s, "needs_transcript": False}

    # rsID
    if re.match(r"^rs\d+$", s, re.I):
        return {"kind": "rsid", "rsid": s.lower(), "raw": s, "needs_transcript": False}

    # genomic coordinate
    g = _parse_genomic(s)
    if g:
        return {"kind": "genomic", "raw": s, "needs_transcript": False, **g}

    # bare c.  (transcript-dependent)
    if s.lower().startswith("c."):
        return {"kind": "bare_c", "cdot": s, "transcript": tx, "raw": s,
                "needs_transcript": tx is None}

    # protein change (transcript-dependent)
    try:
        _parse_protein_change(s)
        return {"kind": "protein", "protein": s, "transcript": tx, "raw": s,
                "needs_transcript": tx is None}
    except ProteinConversionError:
        pass

    return {"kind": "unknown", "raw": s, "needs_transcript": False}


_MISSING_TX_MSG = (
    "no transcript specified. c./protein numbering is transcript-dependent, and "
    "guessing the canonical transcript is exactly how positions end up wrong "
    "(e.g. MAPT is numbered on NM_005910/441aa, not the canonical 758aa isoform). "
    "Fix by giving one of: a genomic coordinate (e.g. chr1:g.155235252A>G), an "
    "rsID, a versioned transcript (e.g. NM_000157.4:c.1448T>C), or add a "
    "'transcript' column. To deliberately use the canonical transcript, put "
    "'canonical' in the transcript field."
)


def _recoder_genomic(hgvs):
    """Ask Ensembl variant_recoder for the genomic (g.) HGVS of any variant HGVS
    or rsID. Returns (chrom, pos, ref, alt). Network call (Ensembl)."""
    data = _get(f"{ENSEMBL}/variant_recoder/human/{hgvs}?")
    # response: list of dicts keyed by allele; find an hgvsg entry
    hgvsg = None
    for entry in (data or []):
        for _, v in entry.items():
            if isinstance(v, dict) and v.get("hgvsg"):
                hgvsg = v["hgvsg"][0]
                break
        if hgvsg:
            break
    if not hgvsg:
        raise VariantInputError(f"Ensembl variant_recoder returned no genomic "
                                f"mapping for '{hgvs}'.")
    g = _parse_genomic(hgvsg)
    if not g:
        raise VariantInputError(f"Could not parse recoder genomic HGVS '{hgvsg}'.")
    return g["chrom"], g["pos"], g["ref"], g["alt"]


def _window_from_genomic(chrom, pos, ref, alt, context=150):
    """Fetch the +/- context window around a genomic SNV and build the
    introduce/revert PRIDICT strings. Verifies the reference base. Network call."""
    start, end = pos - context, pos + context
    seq = _get(f"{ENSEMBL}/sequence/region/human/{chrom}:{start}..{end}:1",
               "text/plain").strip().upper()
    if len(seq) != (2 * context + 1):
        raise VariantInputError(
            f"Fetched {len(seq)} bp for {chrom}:{start}..{end}, expected "
            f"{2*context+1}. Check the genome build / coordinate.")
    up, refbase, down = seq[:context], seq[context], seq[context + 1:]
    if refbase != ref.upper():
        raise VariantInputError(
            f"Reference base at {chrom}:{pos} is '{refbase}', but the variant "
            f"says '{ref}'. Wrong genome build or coordinate?")
    introduce = f"{up}({ref}/{alt}){down}"
    revert = f"{up}({alt}/{ref}){down}"
    return introduce, revert





def make_pridict_inputs(mut, gene_seq, cds_positions, cds_abs, context=150):
    """Given a parsed mutation, build the introduce/revert PRIDICT input strings.

    Format: <up>(<orig>/<edit>)<down>  with `context` bp either side.
    Returns (introduce_str, revert_str). Raises if the reference base doesn't
    match the genome (a strong signal the mutation or coordinate map is off).
    """
    loc = mut["location"]
    if loc not in cds_positions:
        raise MutationParseError(f"CDS position {loc} not found in coordinate map.")
    idx = cds_positions.index(loc)
    abs_pos = cds_abs[idx]

    start = abs_pos - context - 1
    end = abs_pos + context
    ctx = gene_seq[start:end]
    up = ctx[:context]
    down = ctx[-context:]
    ref = ctx[context:-context] if context else ctx

    if mut["edit_type"] == "substitution":
        if ref != mut["original_base"]:
            raise MutationParseError(
                f"Reference base '{ref}' at position {loc} does not match expected "
                f"'{mut['original_base']}'. Check the mutation or transcript.")
        introduce = f"{up}({mut['original_base']}/{mut['mutated_base']}){down}"
        revert = f"{up}({mut['mutated_base']}/{mut['original_base']}){down}"
    elif mut["edit_type"] == "deletion":
        introduce = f"{up}({ref}/-){down}"
        revert = f"{up}(-/{ref}){down}"
    elif mut["edit_type"] == "insertion":
        introduce = f"{up}(+{mut['mutated_base']}){down}"
        revert = f"{up}({mut['mutated_base']}/-){down}"  # best-effort
    else:
        raise MutationParseError(f"Unhandled edit type: {mut['edit_type']}")

    return introduce, revert


def _resolve_variant_row(cls, gene_query, context, canonical_bits):
    """Turn one classified variant into (introduce, revert, trace). Raises
    VariantInputError / ProteinConversionError / MutationParseError on anything
    it can't do confidently — the caller reports these, never silently drops."""
    kind = cls["kind"]

    if kind == "unknown":
        raise VariantInputError(
            f"unrecognised variant format '{cls['raw']}'. Use a genomic coordinate "
            "(chr1:g.123A>G), an rsID, a versioned transcript (NM_/ENST:c.…), or a "
            "protein change with a 'transcript' column.")

    # --- self-coordinate inputs: resolve straight to a genomic window ----------
    if kind == "genomic":
        intro, rev = _window_from_genomic(cls["chrom"], cls["pos"], cls["ref"],
                                          cls["alt"], context)
        return intro, rev, {"coordinate_system": "genomic", "transcript_used": None,
                            "genomic": f"chr{cls['chrom']}:g.{cls['pos']}{cls['ref']}>{cls['alt']}",
                            "warnings": []}
    if kind == "rsid":
        chrom, pos, ref, alt = _recoder_genomic(cls["rsid"])
        intro, rev = _window_from_genomic(chrom, pos, ref, alt, context)
        return intro, rev, {"coordinate_system": "rsID -> genomic", "transcript_used": None,
                            "rsid": cls["rsid"],
                            "genomic": f"chr{chrom}:g.{pos}{ref}>{alt}", "warnings": []}
    if kind == "transcript_c":
        chrom, pos, ref, alt = _recoder_genomic(f"{cls['transcript']}:{cls['cdot']}")
        intro, rev = _window_from_genomic(chrom, pos, ref, alt, context)
        return intro, rev, {"coordinate_system": "transcript c. -> genomic",
                            "transcript_used": cls["transcript"],
                            "genomic": f"chr{chrom}:g.{pos}{ref}>{alt}", "warnings": []}

    # --- transcript-dependent inputs: a transcript is REQUIRED -----------------
    if cls.get("needs_transcript"):
        raise VariantInputError(_MISSING_TX_MSG)

    tx = cls["transcript"]

    # explicit opt-in to the canonical transcript (never the silent default)
    if tx.lower() in ("canonical", "mane"):
        gid, info, cpos, cabs, tid, cds = canonical_bits(gene_query)
        if kind == "protein":
            cdna, conv = normalize_to_cdna(cls["protein"], gene_query, gene_id=gid,
                                           cds=cds, transcript_id=tid)
            warns = conv.get("warnings", [])
        else:
            cdna, warns = cls["cdot"], []
        mut = parse_mutation(cdna)
        intro, rev = make_pridict_inputs(mut, info["gene_sequence"], cpos, cabs, context)
        return intro, rev, {"coordinate_system": "canonical CDS (explicit opt-in)",
                            "transcript_used": tid, "cdna": cdna, "warnings": warns}

    # bare c. with a user-given accession -> recoder maps c. to genomic directly
    if kind == "bare_c":
        chrom, pos, ref, alt = _recoder_genomic(f"{tx}:{cls['cdot']}")
        intro, rev = _window_from_genomic(chrom, pos, ref, alt, context)
        return intro, rev, {"coordinate_system": "transcript c. -> genomic",
                            "transcript_used": tx,
                            "genomic": f"chr{chrom}:g.{pos}{ref}>{alt}", "warnings": []}

    # protein change with a specific ENST -> convert against THAT transcript's CDS
    if re.match(r"^ENST", tx, re.I):
        cds = _get(f"{ENSEMBL}/sequence/id/{tx}?type=cds", "text/plain").strip().upper()
        conv = protein_to_cdna(cls["protein"], gene_query, cds=cds, transcript_id=tx)
        chrom, pos, ref, alt = _recoder_genomic(f"{tx}:{conv['cdna']}")
        intro, rev = _window_from_genomic(chrom, pos, ref, alt, context)
        return intro, rev, {"coordinate_system": "protein -> transcript c. -> genomic",
                            "transcript_used": tx, "cdna": conv["cdna"],
                            "warnings": conv.get("warnings", [])}

    # protein + RefSeq NM (or other non-ENST): Ensembl won't serve its CDS
    raise VariantInputError(
        f"protein change '{cls['protein']}' with a RefSeq transcript '{tx}': supply "
        f"the c. change instead (e.g. {tx}:c.####X>Y), or the matching ENST "
        f"accession, or a genomic/rsID input. Converting protein->c. needs the "
        f"transcript's CDS, which Ensembl doesn't serve for RefSeq NM_ ids.")


def _source_string(raw, trace):
    """Build the source_variant string carried to Stage-3 QC: the user's original
    input plus the canonical single-base change in brackets, e.g.
    'chr17:g.46010389C>T [C>T]'. The bracketed [REF>ALT] is what QC parses to
    confirm the construct installs the mutation the user actually provided."""
    ref = alt = None
    for key in ("genomic", "cdna"):
        m = re.search(r"([ACGT])>([ACGT])\s*$", str(trace.get(key, "")), re.I)
        if m:
            ref, alt = m.group(1).upper(), m.group(2).upper()
            break
    return f"{raw} [{ref}>{alt}]" if ref else str(raw)


def build_pridict_input(mutation_df, gene_query, context=150, mutation_col="gDNA_mutation",
                        check_col="Check", transcript_col="transcript"):
    """Stage 1 entry point for a SINGLE gene.

    Classifies each input and routes it so the coordinate system is never
    guessed. Genomic (g.), rsID and versioned-transcript (NM_/ENST:c.) inputs
    resolve straight to a genomic window. Bare c. or protein changes REQUIRE a
    transcript (inline accession or a 'transcript' column); without one they are
    reported in `skipped`, never silently pushed through canonical. Every
    resolved row is recorded in report['resolved'] for traceability.
    """
    _canon = {}

    def _canonical_bits(gq):
        if gq not in _canon:
            gid = resolve_gene_id(gq)
            info, cpos, cabs = build_coordinate_map(gid)
            tid, cds = _fetch_cds(gid)
            _canon[gq] = (gid, info, cpos, cabs, tid, cds)
        return _canon[gq]

    rows, skipped, warnings, resolved, source_map = [], [], [], [], {}
    for i, r in mutation_df.iterrows():
        raw = r.get(mutation_col)
        if check_col in mutation_df.columns and str(r.get(check_col)).strip().lower() == "x":
            skipped.append((str(raw), "flagged 'x' for manual curation"))
            continue
        tx = r.get(transcript_col) if transcript_col in mutation_df.columns else None
        cls = classify_variant(raw, transcript=tx)
        try:
            intro, revert, trace = _resolve_variant_row(cls, gene_query, context, _canonical_bits)
        except (VariantInputError, ProteinConversionError, MutationParseError) as e:
            skipped.append((str(raw), str(e)))
            continue
        if trace.get("warnings"):
            warnings.append((str(raw), "; ".join(trace["warnings"])))
        tag = re.sub(r"[^A-Za-z0-9]", "", str(raw))
        # source_variant: the ORIGINAL user input plus its canonical base change,
        # carried independently of PRIDICT so Stage-3 QC can confirm the construct
        # installs the mutation the user actually asked for.
        src = _source_string(str(raw), trace)
        for nm in (f"{tag}_intro_{i}", f"{tag}_revert_{i}"):
            source_map[nm] = src
        rows.append({"sequence_name": f"{tag}_intro_{i}", "editseq": intro})
        rows.append({"sequence_name": f"{tag}_revert_{i}", "editseq": revert})
        resolved.append({"input": str(raw), "kind": cls["kind"], **trace})

    batch_df = pd.DataFrame(rows)
    report = {
        "gene_name": gene_query, "gene_query": gene_query, "strand": None,
        "n_input_mutations": len(mutation_df),
        "n_sequences_generated": len(batch_df),
        "n_skipped": len(skipped), "skipped": skipped,
        "warnings": warnings,
        "resolved": resolved,
        "source_map": source_map,
    }
    return batch_df, report



def build_pridict_input_multigene(mutation_df, context=150,
                                  mutation_col="gDNA_mutation", gene_col="gene",
                                  check_col="Check"):
    """Stage 1 entry point for MULTIPLE genes in one table.

    The table must have a `gene_col` (default 'gene') naming the gene for each
    row, plus `mutation_col`. Rows are grouped by gene, each gene fetched from
    Ensembl once, then processed. Returns (batch_df, report) with a per-gene
    breakdown. Mutations can be c. or protein notation.
    """
    if gene_col not in mutation_df.columns:
        raise ValueError(
            f"Multi-gene mode needs a '{gene_col}' column naming the gene per row. "
            f"Found columns: {list(mutation_df.columns)}. "
            f"For a single gene, use single-gene mode instead.")

    all_rows, per_gene, all_skipped, all_warnings, all_resolved = [], [], [], [], []
    source_map = {}
    for gene_query, sub in mutation_df.groupby(gene_col, dropna=False):
        # rows with no gene are fine as long as they are self-coordinate
        # (genomic / rsID / accession:c.); pass an empty gene through.
        gq = "" if (gene_query is None or (isinstance(gene_query, float) and pd.isna(gene_query))) else str(gene_query)
        try:
            batch_df, rep = build_pridict_input(
                sub.reset_index(drop=True), gq, context=context,
                mutation_col=mutation_col, check_col=check_col)
        except Exception as e:
            per_gene.append({"gene": gq or "(no gene)", "error": f"{type(e).__name__}: {e}",
                             "n_sequences": 0})
            all_skipped += [(f"[{gq or 'no gene'}] {m}", "gene-level failure: " + str(e))
                            for m in sub[mutation_col].astype(str)]
            continue
        # prefix sequence names with gene to keep them unique across genes
        prefix = (gq + "_") if gq else ""
        if len(batch_df):
            batch_df = batch_df.copy()
            batch_df["sequence_name"] = prefix + batch_df["sequence_name"]
            all_rows.append(batch_df)
        # re-key the per-gene source_map with the SAME prefix so it still matches
        for nm, sv in rep.get("source_map", {}).items():
            source_map[prefix + nm] = sv
        per_gene.append({"gene": rep["gene_name"] or "(no gene)", "strand": rep.get("strand"),
                         "n_sequences": rep["n_sequences_generated"],
                         "n_skipped": rep["n_skipped"]})
        all_skipped += [(f"[{gq or 'no gene'}] {m}", why) for m, why in rep["skipped"]]
        all_warnings += [(f"[{gq or 'no gene'}] {m}", w) for m, w in rep.get("warnings", [])]
        all_resolved += [{"gene": gq, **rv} for rv in rep.get("resolved", [])]

    combined = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame(
        columns=["sequence_name", "editseq"])
    report = {
        "mode": "multi-gene",
        "n_genes": mutation_df[gene_col].nunique(),
        "per_gene": per_gene,
        "n_sequences_generated": len(combined),
        "source_map": source_map,
        "n_skipped": len(all_skipped), "skipped": all_skipped,
        "warnings": all_warnings,
        "resolved": all_resolved,
    }
    return combined, report


# ---------------------------------------------------------------------------
# Stage 1-snapgene — Read mutations directly from annotated SnapGene .dna files
# ---------------------------------------------------------------------------
#
# This is the MOST RELIABLE input path: the user has already annotated each
# mutation as a feature on the sequence in SnapGene, so we don't need Ensembl,
# transcript alignment, isoform disambiguation, or CDS coordinate systems.
# We read the feature's location + the sequence around it directly. Promoter
# mutations (e.g. TERT -146C>T) work naturally because they're just positions
# on the provided sequence.
#
# The base change is taken from the feature LABEL. We parse several label
# conventions; if a label carries no base change, we DO NOT guess — we report
# it as needing curation (per the user's choice to trust only label info).

def _parse_label_base_change(label, ref_seq):
    """Extract the base change from a SnapGene feature label.

    Recognises conventions like:
      '-146C/T'            -> ('C','T', None)      single base, offset unknown
      'P301L_C>T'          -> ('C','T', None)
      'N370S_AAC>AGC'      -> ('A','G', 1)         codon change: offset 1 within feature
      'L444P'              -> None                 no base info in label
    Returns (orig, alt, offset) where offset is the 0-based position of the
    changed base WITHIN the feature/codon (or None if unknown), or None if the
    label carries no usable base change. The offset removes the ambiguity when a
    codon contains repeated bases (e.g. finding the right 'A' in 'AAC').
    """
    lab = str(label)

    # codon change like AAC>AGC or GTG>ATG -> use the differing POSITION
    m = re.search(r"([ACGT]{2,})\s*>\s*([ACGT]{2,})", lab, re.I)
    if m:
        a, b = m.group(1).upper(), m.group(2).upper()
        if len(a) == len(b):
            diffs = [(i, a[i], b[i]) for i in range(len(a)) if a[i] != b[i]]
            if len(diffs) == 1:
                i, o, n = diffs[0]
                return o, n, i  # offset within codon is unambiguous

    # single-base change:  X/Y  or  X>Y
    m = re.search(r"(?<![A-Za-z])([ACGT])\s*[/>]\s*([ACGT])(?![A-Za-z])", lab, re.I)
    if m:
        return m.group(1).upper(), m.group(2).upper(), None

    return None  # no base change in label — do not guess


def read_snapgene_mutations(path, context=150):
    """Parse an annotated SnapGene .dna file into PRIDICT inputs.

    Finds misc_feature annotations whose label looks like a point mutation,
    reads the reference base and (from the label) the alternate base, and builds
    introduce/revert PRIDICT input strings from the surrounding sequence.

    Returns (rows, skipped) where rows are dicts {sequence_name, editseq} and
    skipped is a list of (label, reason). Features whose label carries no base
    change are skipped and reported (never guessed).
    """
    try:
        from Bio import SeqIO
    except ImportError:
        raise RuntimeError("Biopython is required to read SnapGene files "
                           "(pip install biopython).")

    rec = SeqIO.read(path, "snapgene")
    seq = str(rec.seq).upper()
    rows, skipped = [], []

    # a label is a candidate mutation if it looks like -NNNx/y or PnnnX or has X>Y
    def _looks_like_mutation(lab):
        return bool(
            re.search(r"-?\d+\s*[ACGT]\s*[/>]\s*[ACGT]", lab, re.I) or
            re.search(r"[ACGT]{2,}\s*>\s*[ACGT]{2,}", lab, re.I) or
            re.match(r"^[A-Z]\d+[A-Z]", lab)  # protein change like L444P/R406W
        )

    seen = set()
    for feat in rec.features:
        if feat.type != "misc_feature":
            continue
        label = feat.qualifiers.get("label", [""])[0]
        if not label or not _looks_like_mutation(label):
            continue
        start, end = int(feat.location.start), int(feat.location.end)
        key = (label, start, end)
        if key in seen:
            continue
        seen.add(key)

        ref_here = seq[start:end]
        change = _parse_label_base_change(label, ref_here)
        if change is None:
            skipped.append((label,
                "Label has no base change (only a protein change). Re-annotate in "
                "SnapGene as e.g. 'L444P_T>C', or supply the base change."))
            continue
        orig_base, alt_base, offset = change

        # locate the exact edited base.
        if len(ref_here) == 1:
            # feature pinpoints the single edited base already — trust it,
            # regardless of whether the label was single-base or codon-level.
            edit_pos = start
        elif offset is not None:
            # codon-level label told us WHICH position changed — unambiguous
            if offset >= len(ref_here):
                skipped.append((label,
                    f"Codon offset {offset} is outside the {len(ref_here)}-bp "
                    "feature. Check that the feature spans the whole codon, or "
                    "annotate just the single edited base."))
                continue
            edit_pos = start + offset
        else:
            # single-base label but multi-bp feature: search, accept if unambiguous
            hits = [start + i for i, b in enumerate(ref_here) if b == orig_base]
            if len(hits) != 1:
                skipped.append((label,
                    f"Base '{orig_base}' appears {len(hits)} times in feature "
                    f"'{ref_here}'. Use a codon label like 'AAC>AGC' so the tool "
                    "knows which position changed, or annotate the single base."))
                continue
            edit_pos = hits[0]

        # sanity: reference base under the feature must match the label's orig
        if seq[edit_pos] != orig_base:
            skipped.append((label,
                f"Reference base '{seq[edit_pos]}' at the annotated position does "
                f"not match label's '{orig_base}'. Check the annotation."))
            continue

        up = seq[edit_pos - context:edit_pos]
        down = seq[edit_pos + 1:edit_pos + 1 + context]
        if len(up) < context or len(down) < context:
            skipped.append((label,
                f"Not enough flanking sequence ({context} bp) around the mutation "
                "in this file. Use a longer SnapGene region."))
            continue

        introduce = f"{up}({orig_base}/{alt_base}){down}"
        revert = f"{up}({alt_base}/{orig_base}){down}"
        tag = re.sub(r"[^A-Za-z0-9]", "", label)[:30]
        src = f"{label} [{orig_base}>{alt_base}]"
        # `_edit_pos/_orig/_alt` are carried for callers that need to place the
        # edit back on the SnapGene coordinate (e.g. the codon-block builder).
        # build_pridict_input_from_snapgene only reads sequence_name/editseq.
        extra = {"_source": src, "_edit_pos": edit_pos,
                 "_orig": orig_base, "_alt": alt_base}
        rows.append({"sequence_name": f"{tag}_intro", "editseq": introduce, **extra})
        rows.append({"sequence_name": f"{tag}_revert", "editseq": revert, **extra})

    return rows, skipped


def build_pridict_input_from_snapgene(paths, context=150):
    """Stage 1 entry point for one or more SnapGene .dna files.

    `paths` is a list of file paths. Each file's annotated mutations become
    PRIDICT inputs. Returns (batch_df, report) with a per-file breakdown and a
    transparent list of anything skipped.
    """
    all_rows, per_file, all_skipped = [], [], []
    for path in paths:
        name = path.split("/")[-1]
        try:
            rows, skipped = read_snapgene_mutations(path, context=context)
        except Exception as e:
            per_file.append({"file": name, "error": f"{type(e).__name__}: {e}",
                             "n_sequences": 0})
            continue
        # prefix names with file stem to keep unique across files
        stem = re.sub(r"\.dna$", "", name)
        stem = re.sub(r"[^A-Za-z0-9]", "", stem)[:20]
        for r in rows:
            r["sequence_name"] = f"{stem}_{r['sequence_name']}"
        all_rows += rows
        all_skipped += [(f"[{name}] {lab}", why) for lab, why in skipped]
        per_file.append({"file": name, "n_sequences": len(rows),
                         "n_skipped": len(skipped)})

    source_map = {r["sequence_name"]: r.get("_source", "") for r in all_rows}
    batch_df = pd.DataFrame(
        [{"sequence_name": r["sequence_name"], "editseq": r["editseq"]} for r in all_rows]
    ) if all_rows else pd.DataFrame(columns=["sequence_name", "editseq"])
    report = {
        "mode": "snapgene",
        "n_files": len(paths),
        "per_file": per_file,
        "n_sequences_generated": len(batch_df),
        "n_skipped": len(all_skipped), "skipped": all_skipped,
        "source_map": source_map,
    }
    return batch_df, report


# ---------------------------------------------------------------------------
# Stage 1-bystander — Generate silent bystander variants (SnapGene-based)
# ---------------------------------------------------------------------------
#
# A silent bystander is an EXTRA synonymous change introduced alongside the real
# edit — it changes DNA bases but NOT the protein, and can raise prime-editing
# efficiency (e.g. by evading mismatch repair). Only makes sense for coding
# substitutions where we know the reading frame.
#
# Reliability strategy (no internet, no transcript alignment):
#   - The reading frame comes from the SnapGene annotation itself: either a 3-bp
#     codon feature, or a single-base feature whose codon boundary is recovered
#     by trying the 3 frames and keeping the one that translates to the mutation's
#     original amino acid (unambiguous when combined with the target AA).
#   - Every generated bystander is VERIFIED to be synonymous (translates to the
#     same protein window) before being emitted. Non-synonymous combinations are
#     discarded, never shipped.
#   - Promoter / non-coding mutations are refused (no codons).

# reverse codon table: amino acid -> list of synonymous codons
_EDITSEQ_RE = re.compile(r"\(([ACGTacgt-]+)/([ACGTacgt-]+)\)")
_SYN_CODONS = {}
for _c, _a in _CODON_TABLE.items():
    _SYN_CODONS.setdefault(_a, []).append(_c)


def _bystander_rc(s):
    return str(s).translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def _cds_aa_window(fullseq, cds_start, cds_end, strand, lo, hi):
    """Translate the CDS codons overlapping the plus-strand window [lo, hi) in the
    CORRECT reading frame and strand (from the CDS annotation). Returns a dict
    {codon_plus_start: amino_acid}. This is the INDEPENDENT truth used to decide
    whether bystander edits are silent — it does not rely on how the edits were
    enumerated, so it catches wrong-frame AND wrong-strand mistakes."""
    aa = {}
    if cds_start is None or cds_end is None:
        return aa
    if strand == -1:
        # codons are anchored at cds_end; codon i occupies plus [cds_end-3(i+1), cds_end-3i)
        i_lo = max(0, (cds_end - hi) // 3)
        i_hi = (cds_end - lo) // 3 + 1
        for i in range(i_lo, i_hi):
            ps = cds_end - 3 * (i + 1)
            if ps < cds_start or ps + 3 > cds_end:
                continue
            aa[ps] = _CODON_TABLE.get(_bystander_rc(fullseq[ps:ps + 3]), "?")
    else:
        first = cds_start + ((max(lo, cds_start) - cds_start) // 3) * 3
        for ps in range(first, hi, 3):
            if ps < cds_start or ps + 3 > cds_end:
                continue
            aa[ps] = _CODON_TABLE.get(fullseq[ps:ps + 3], "?")
    return aa


def _bystanders_silent(orig_seq, byst_seq, cds_start, cds_end, strand, lo, hi):
    """True iff `byst_seq` (original amino acids + candidate bystander edits)
    translates to the SAME protein as `orig_seq` across the window — i.e. every
    bystander change is synonymous in the real CDS frame/strand. If no CDS info
    is available, returns None (can't verify)."""
    if cds_start is None:
        return None
    a = _cds_aa_window(orig_seq, cds_start, cds_end, strand, lo, hi)
    b = _cds_aa_window(byst_seq, cds_start, cds_end, strand, lo, hi)
    return a == b


class BystanderError(Exception):
    pass





def _revcomp(s):
    return str(s).translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


# ---------------------------------------------------------------------------
# Stage 1-bystander (OptiPrime native) — codon OPTION BLOCKS
# ---------------------------------------------------------------------------
# PRIDICT2 wants every bystander combination as its own input row (the add-on
# enumerates them). OptiPrime wants the opposite: ONE run per protospacer with
# the per-codon synonymous options handed over as `edit_segments`, which it then
# combines itself and prunes in four rounds. Feeding it the enumerated rows is
# both wrong for its interface and hundreds of times more expensive.
#
# Difference from the add-on, deliberate: a codon is only offered if it lies
# ENTIRELY inside the annotated CDS. The add-on walks the reading frame linearly
# along the input sequence, so where the +-window crosses a splice junction it
# emits "synonymous" variants for intronic bases. For GBA1 N370S the CDS starts
# at the edited codon, so its -1 "codon" is the intron's last 3 bases, TAG, and
# the offered TAA/TGA both destroy the 3' splice acceptor AG. Those changes pass
# a translate-and-compare silence check precisely because they touch no coding
# base, which is what makes them dangerous rather than harmless.

def _syn_codons_for(codon):
    """Synonymous codons for `codon`, the codon itself first."""
    aa = _CODON_TABLE.get(codon)
    if aa is None or aa == "*":
        return [codon]
    opts = sorted(_SYN_CODONS.get(aa, [codon]))
    return [codon] + [c for c in opts if c != codon]


def build_codon_blocks_from_snapgene(paths, window_codons=2, context=150,
                                     directions=("intro", "revert")):
    """Codon option blocks for OptiPrime's `edit_segments`, per input sequence.

    Returns ({sequence_name: [(index_in_editseq, [option, ...]), ...]}, report).
    Indices and options are in the SAME orientation as the editseq that
    `build_pridict_input_from_snapgene` produces for that name, so they can be
    dropped straight into the OptiPrime run builder.
    """
    from Bio import SeqIO            # imported lazily, as elsewhere in this module
    blocks_by_name, per_file, skipped = {}, [], []
    for path in paths:
        fname = path.split("/")[-1]
        stem = re.sub(r"\.dna$", "", fname)
        stem = re.sub(r"[^A-Za-z0-9]", "", stem)[:20]
        try:
            rec = SeqIO.read(path, "snapgene")
        except Exception as e:
            per_file.append({"file": fname, "error": str(e)})
            continue
        seq = str(rec.seq).upper()
        cds_list = [(int(f.location.start), int(f.location.end),
                     (f.location.strand if f.location.strand in (1, -1) else 1))
                    for f in rec.features if f.type == "CDS"]

        rows, _sk = read_snapgene_mutations(path, context=context)
        n_here = 0
        for r in rows:
            name = f"{stem}_{r['sequence_name']}"
            direction = "intro" if name.endswith("_intro") else "revert"
            if direction not in directions:
                continue
            ep = r["_edit_pos"]
            cds_hit = next(((s, e, st) for s, e, st in cds_list if s <= ep < e), None)
            if cds_hit is None:
                skipped.append((name, "No CDS covers this edit; silent bystanders "
                                      "need a reading frame. Non-coding variant?"))
                continue
            cds_s, cds_e, strand = cds_hit
            if strand == -1:
                k = (cds_e - ep - 1) // 3
                main_cs = cds_e - 3 * (k + 1)
            else:
                main_cs = cds_s + ((ep - cds_s) // 3) * 3

            # the base this editseq ENDS UP with at the edit position
            edited_base = r["_alt"] if direction == "intro" else r["_orig"]
            edited_seq = seq[:ep] + edited_base + seq[ep + 1:]

            blocks, n_out = [], 0
            for j in range(-window_codons, window_codons + 1):
                cs = main_cs + 3 * j
                if cs < cds_s or cs + 3 > cds_e:
                    n_out += 1
                    continue                      # outside the CDS -> never offer
                plus = edited_seq[cs:cs + 3]
                coding = _bystander_rc(plus) if strand == -1 else plus
                opts_coding = _syn_codons_for(coding)
                if len(opts_coding) < 2 and j != 0:
                    continue                      # nothing to vary here
                opts_plus = ([_bystander_rc(c) for c in opts_coding]
                             if strand == -1 else opts_coding)
                blocks.append((cs - ep + context, opts_plus))
            if n_out:
                skipped.append((name, f"{n_out} codon(s) in the +-{window_codons} "
                                      f"window lie outside the CDS and were not "
                                      f"offered (splice-site protection)."))
            if blocks:
                blocks_by_name[name] = blocks
                n_here += 1
        per_file.append({"file": fname, "n_sequences": n_here})

    n_combos = {}
    for name, blocks in blocks_by_name.items():
        n = 1
        for _i, opts in blocks:
            n *= len(opts)
        n_combos[name] = n
    return blocks_by_name, {"per_file": per_file, "skipped": skipped,
                            "window_codons": window_codons,
                            "n_sequences": len(blocks_by_name),
                            "n_combinations": n_combos}


# ---------------------------------------------------------------------------
# Stage 1-bystander — Silent bystander generation via the OFFICIAL PRIDICT2
# add-on (addons/silentbystander).  We do NOT re-implement the enumeration;
# we load the add-on's own functions from the user's PRIDICT2 install at run
# time and call them, so the silent-bystander logic is byte-for-byte the
# authors' code (and stays in sync with `git pull`).  Our job is only to feed
# it correctly: build the intro/revert PRIDICT input in CODING orientation and
# in-frame (compute ORF_start from the SnapGene CDS annotation; reverse-
# complement for minus-strand genes).  This is the reference-notebook workflow.
# ---------------------------------------------------------------------------

_ADDON_CACHE = {}


def _load_silentbystander_addon(pridict_home="~/PRIDICT2"):
    """Load the official silent-bystander functions from the user's PRIDICT2
    install (addons/silentbystander/notebook_silent_bystander_input.ipynb).
    Only cells that contain purely imports / function / class definitions are
    executed, so the example/driver cells with hard-coded paths never run.
    Returns a namespace dict with `bystander_creation_for_pridict` etc."""
    import ast
    home = os.path.expanduser(pridict_home)
    nb_path = os.path.join(home, "addons", "silentbystander",
                           "notebook_silent_bystander_input.ipynb")
    key = os.path.abspath(nb_path)
    if key in _ADDON_CACHE:
        return _ADDON_CACHE[key]
    if not os.path.exists(nb_path):
        raise BystanderError(
            "Official silent-bystander add-on not found at " + nb_path + ". "
            "It ships with PRIDICT2 (addons/silentbystander/). Check that "
            "'PRIDICT folder' points to your PRIDICT2 install, or `git pull` it.")
    import json
    nb = json.load(open(nb_path))
    ns = {}
    for c in nb.get("cells", []):
        if c.get("cell_type") != "code":
            continue
        src = "".join(c.get("source", []))
        if not src.strip():
            continue
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        if tree.body and all(isinstance(n, (ast.Import, ast.ImportFrom,
                             ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef))
                             for n in tree.body):
            exec(compile(tree, nb_path, "exec"), ns)  # noqa: S102 (official code)
    if "bystander_creation_for_pridict" not in ns:
        raise BystanderError(
            "Loaded the add-on notebook but could not find "
            "`bystander_creation_for_pridict`. The add-on format may have changed.")
    _ADDON_CACHE[key] = ns
    return ns


def _changes_outside_cds(editseq, genomic_seq, cds_start, cds_end):
    """Genomic positions this bystander changes that are NOT inside the CDS.

    The add-on walks the reading frame linearly along its input, so where the
    +-window crosses a splice junction it offers "synonymous" variants for
    intronic bases. Those pass a translate-and-compare silence check precisely
    because they touch no coding base — which is what makes them dangerous
    rather than harmless. For GBA1 N370S the CDS starts at the edited codon, so
    its -1 "codon" is the intron's last three bases and the offered variants
    change the 3' splice acceptor AG itself.

    The add-on re-emits the editseq with its own flank lengths (299 nt where the
    main-edit row has 301, starting 2 bp further along), so the mapping back to
    genomic coordinates is found by ALIGNING rather than assumed from `context`.
    Anchoring on the first 40 bp is safe: the edit sits ~150 bp in, so the
    anchor is identical in both the reference and the (possibly mutant) input.
    """
    m = _EDITSEQ_RE.search(str(editseq))
    if not m:
        return []
    ref, alt = m.group(1), m.group(2)
    if len(ref) != len(alt) or "-" in (ref, alt):
        return []                      # indel: no positional mapping, leave it
    es = str(editseq)
    un = (es[:m.start()] + ref + es[m.end():]).upper()

    j = genomic_seq.find(un[:40])
    flip = False
    if j < 0:
        j = genomic_seq.find(_bystander_rc(un)[:40])
        flip = True
    if j < 0:
        return []                      # cannot place it; do not guess

    out = []
    for k in range(len(ref)):
        if ref[k].upper() == alt[k].upper():
            continue
        i = m.start() + k
        g = (j + len(un) - 1 - i) if flip else (j + i)
        if not (cds_start <= g < cds_end):
            out.append(g)
    return sorted(out)


def build_bystander_from_snapgene(paths, window_codons=2, context=150,
                                  max_variants_per_mut=400, total_edit_limit=40,
                                  directions=("revert",), pridict_home="~/PRIDICT2",
                                  max_edit_length=10):
    """Generate silent-bystander PRIDICT inputs from annotated SnapGene files,
    using the OFFICIAL PRIDICT2 silent-bystander add-on for the enumeration.

    For each coding substitution we build the intro/revert PRIDICT input in the
    CDS's coding orientation and in reading frame (ORF_start from the CDS), then
    hand it to the add-on's `bystander_creation_for_pridict` with silent='yes',
    change_edit_bases='no'. Promoter / non-coding mutations (no CDS) are skipped.
    `window_codons` maps to the add-on's `silent_surrounding_AA_nr`.
    """
    from Bio import SeqIO

    def _aa_from_label(label):
        m = re.match(r"^([A-Z])(\d+)([A-Z])", str(label))
        return (m.group(1), m.group(3)) if m else (None, None)

    addon = _load_silentbystander_addon(pridict_home)
    bystander_creation = addon["bystander_creation_for_pridict"]

    all_rows, per_file, skipped = [], [], []
    byst_source_map = {}
    n_outside = 0
    for path in paths:
        name = path.split("/")[-1]
        stem = re.sub(r"\.dna$", "", name)
        stem = re.sub(r"[^A-Za-z0-9]", "", stem)[:20]
        try:
            rec = SeqIO.read(path, "snapgene")
        except Exception as e:
            per_file.append({"file": name, "error": str(e)})
            continue
        seq = str(rec.seq).upper()
        cds_list = [(int(f.location.start), int(f.location.end),
                     (f.location.strand if f.location.strand in (1, -1) else 1))
                    for f in rec.features if f.type == "CDS"]

        def _cds_for(pos):
            for cs, ce, st in cds_list:
                if cs <= pos < ce:
                    return cs, ce, st
            return None

        n_here = 0
        for feat in rec.features:
            if feat.type != "misc_feature":
                continue
            label = feat.qualifiers.get("label", [""])[0]
            orig_aa, target_aa = _aa_from_label(label)
            if orig_aa is None:
                if re.search(r"-?\d+\s*[ACGT]\s*[/>]\s*[ACGT]", label, re.I) and \
                   not re.match(r"^[A-Z]\d+[A-Z]", label):
                    skipped.append((f"[{name}] {label}",
                        "Non-coding / promoter mutation — no codons, bystander not applicable."))
                continue

            s, e = int(feat.location.start), int(feat.location.end)
            change = _parse_label_base_change(label, seq[s:e])
            if change is None:
                skipped.append((f"[{name}] {label}",
                    "Label has no base change; cannot place the edit for bystander."))
                continue
            orig_base, alt_base, offset = change

            try:
                if (e - s) == 1:
                    ep = s
                elif (e - s) == 3:
                    ep = s + (offset if offset is not None else
                              next((i for i in range(3) if seq[s + i] == orig_base), 0))
                else:
                    ep = s + offset if offset is not None else s

                cds_hit = _cds_for(ep)
                if cds_hit is None:
                    raise BystanderError(
                        "No CDS feature covers this edit; cannot establish the reading "
                        "frame required for silent bystanders. Add a CDS annotation.")
                cds_s, cds_e, strand = cds_hit

                # reading frame / edit offset from the CDS
                if strand == -1:
                    k = (cds_e - ep - 1) // 3
                    codon_start = cds_e - 3 * (k + 1)
                    edit_off_plus = ep - codon_start
                    edit_off_coding = 2 - edit_off_plus
                else:
                    codon_start = cds_s + ((ep - cds_s) // 3) * 3
                    edit_off_plus = ep - codon_start
                    edit_off_coding = edit_off_plus
                if not (0 <= edit_off_coding < 3):
                    raise BystanderError("Edit falls outside a clean CDS codon.")
                orf_start = (3 - edit_off_coding) % 3

                up = seq[ep - context:ep]
                down = seq[ep + 1:ep + 1 + context]
                if len(up) < context or len(down) < context:
                    raise BystanderError(
                        f"Need >= {context} bp flanking; use a longer SnapGene region.")

                # coding-orientation PRIDICT inputs for each direction
                if strand == -1:
                    c_up, c_down = _revcomp(down), _revcomp(up)
                    o_cod, a_cod = _revcomp(orig_base), _revcomp(alt_base)
                else:
                    c_up, c_down = up, down
                    o_cod, a_cod = orig_base, alt_base
                inputs = {
                    "intro": f"{c_up}({o_cod}/{a_cod}){c_down}",
                    "revert": f"{c_up}({a_cod}/{o_cod}){c_down}",
                }

                tag = re.sub(r"[^A-Za-z0-9]", "", label)[:24]
                byst_src = f"{label} [{orig_base}>{alt_base}]"
                made = 0
                for direction in directions:
                    if direction not in inputs:
                        continue
                    outdf = bystander_creation(
                        inputs[direction], window_codons, orf_start,
                        f"{stem}_{tag}_{direction}", 94, total_edit_limit,
                        max_edit_length, silent="yes", change_edit_bases="no")
                    for vi, row in enumerate(outdf.itertuples(index=False)):
                        nchg = int(getattr(row, "total_nr_of_base_changes"))
                        outside = _changes_outside_cds(row.editseq, seq,
                                                       cds_s, cds_e)
                        if outside:
                            n_outside += 1
                            continue
                        nm = f"{stem}_{tag}_{direction}_byst{vi}_{nchg}nt"
                        all_rows.append({"sequence_name": nm, "editseq": row.editseq})
                        byst_source_map[nm] = byst_src
                        made += 1
                        if made >= max_variants_per_mut * len(directions):
                            break
                n_here += made
            except BystanderError as ex:
                skipped.append((f"[{name}] {label}", str(ex)))
                continue
            except Exception as ex:
                skipped.append((f"[{name}] {label}",
                                f"add-on error: {type(ex).__name__}: {ex}"))
                continue

        per_file.append({"file": name, "n_bystander_variants": n_here})

    batch_df = pd.DataFrame(all_rows) if all_rows else pd.DataFrame(
        columns=["sequence_name", "editseq"])
    report = {
        "mode": "bystander",
        "n_files": len(paths),
        "per_file": per_file,
        "n_variants_generated": len(batch_df),
        "n_skipped": len(skipped), "skipped": skipped,
        "window_codons": window_codons,
        "source_map": byst_source_map,
        "engine": "official PRIDICT2 silentbystander add-on",
        # Variants the add-on produced that change bases outside the annotated
        # CDS. They translate identically (they touch no coding base), so the
        # silence check cannot see them, but they can hit a splice site.
        "n_dropped_outside_cds": n_outside,
    }
    return batch_df, report


# ---------------------------------------------------------------------------
# Stage 2 — Assemble self-targeting library from PRIDICT summary
# ---------------------------------------------------------------------------

# Default fixed elements — these reproduce the Schwank-lab CSF1R design exactly.
# All are exposed as parameters so a different vector / cloning system can be used.
DEFAULTS = {
    "fiveprimeoverhang": "GTGGAAAGGACGAAACACC",
    "bsmbi_spacer": "GTTTAGAGACGGTAGCTGTCGTCTCTGTGC",
    "tevopreQ1": "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA",
    "polyT": "TTTTTTT",
    "threeprimeoverhang": "GTGACTCCTATGACGCTTCT",
    "spacerlength_with_g": 20,
    "barcodelength": 6,
    "selftargetmaxlength": 300,
    "library_size": 1000,
    "top_n_per_mutation": 10,
}


def _occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start)
        if start == -1:
            return count
        count += 1
        start += 1


def _revcomp_dna(seq):
    return str(seq).translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def build_library(summary_df, params=None, seed=42):
    """Stage 2 entry point.

    summary_df: a PRIDICT2.0 K562 summary. Must contain the pegRNA component
    columns (Spacer-Sequence, RTrevcomp, PBSrevcomp, PBSlength, RTlength),
    plus wide_initial_target and Original_Sequence (used to recover the
    genomic downstream sequence and cut the target). params: overrides for
    DEFAULTS. Returns (library_df, report).

    This assembles the full self-targeting construct, exactly reproducing the
    notebook logic (Cell 22/23): recover downstream sequence, size the target,
    build the construct, then run QC (length, poly-base, stray BsmBI/primer
    sites). Rows that fail hard requirements are dropped and reported, not
    silently kept.
    """
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    random.seed(seed)

    df = summary_df.copy().reset_index(drop=True)

    required = ["Spacer-Sequence", "RTrevcomp", "PBSrevcomp", "PBSlength",
                "RTlength", "wide_initial_target", "Original_Sequence"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError("Summary file is missing columns needed to build "
                         f"constructs: {missing}. Is this a PRIDICT2.0 K562 summary?")

    static_len = (len(p["fiveprimeoverhang"]) + len(p["bsmbi_spacer"]) +
                  len(p["tevopreQ1"]) + len(p["polyT"]) + len(p["threeprimeoverhang"]) +
                  p["spacerlength_with_g"])

    # ---- classify intro/revert + mutation label from sequence_name ----
    # bystander sequences are named like  <mother>_byst<idx>_<n>nt  (and may also
    # carry an _intro/_revert). We strip the bystander suffix so every bystander
    # variant is grouped back UNDER its mother mutation, and flag it as bystander.
    def _is_bystander(name):
        return bool(re.search(r"_byst\d+_\d+nt", str(name)))

    def _direction(name):
        n = str(name).lower()
        if "_revert" in n:
            return "revert"
        if "_intro" in n:
            return "intro"
        return "unknown"

    def _mutation(name):
        s = str(name)
        s = re.sub(r"_byst\d+_\d+nt.*$", "", s)   # strip bystander suffix first
        s = re.sub(r"_(intro|revert).*$", "", s)  # then strip direction
        return s

    if "sequence_name" in df.columns:
        df["intro_or_revert"] = df["sequence_name"].apply(_direction)
        df["mutation"] = df["sequence_name"].apply(_mutation)
        # PRIDICT2 gives each bystander combination its own input sequence, so
        # the name carries the flag. OptiPrime enumerates the combinations
        # INSIDE one run and labels them in an `edit_name` column, so the name is
        # identical for all of them — trust the column when the summary has one.
        if "bystander" in df.columns:
            df["bystander"] = df["bystander"].fillna("no")
        else:
            df["bystander"] = df["sequence_name"].apply(
                lambda n: "yes" if _is_bystander(n) else "no")
    else:
        df["intro_or_revert"] = "unknown"
        df["mutation"] = df.index.astype(str)
        df["bystander"] = "no"

    # ---- pre-select top-N per mutation BEFORE the expensive assembly / barcode
    # steps, so ONLY the winners are ever assigned a barcode. This is what keeps
    # the barcode pool from exploding: a few hundred constructs need a barcode,
    # not every one of the (possibly tens of thousands of) bystander candidates.
    # A buffer above top_n is kept here so the downstream / QC drops below don't
    # leave fewer than top_n; the set is trimmed to EXACTLY top_n at the very end.
    score_col = next((c for c in df.columns if "K562" in c and "Score" in c), None)
    _top_n = p.get("top_n_per_mutation")

    def _select_top(d, n):
        if not n:
            return d
        if score_col and score_col in d.columns:
            d = d.sort_values(score_col, ascending=False)
        gcols = ["mutation"]
        if "intro_or_revert" in d.columns and d["intro_or_revert"].nunique() > 1:
            gcols.append("intro_or_revert")
        return d.groupby(gcols, group_keys=False).head(int(n)).reset_index(drop=True)

    if _top_n:
        df = _select_top(df, int(_top_n) + 20)   # over-select with a buffer

    # ---- recover downstream genomic sequence (notebook Cell 22) ----
    notes = []
    df["extension_length"] = df["PBSlength"] + df["RTlength"]
    downstream = []
    for _, row in df.iterrows():
        wit = str(row["wide_initial_target"])
        orig = str(row["Original_Sequence"])
        end19 = wit[-19:]
        pos = orig.find(end19)
        if pos == -1 or orig.find(end19, pos + 1) != -1:
            downstream.append(None)  # not found or ambiguous
        else:
            downstream.append(orig[pos + 19:])
    df["downstream_sequence"] = downstream

    n_before = len(df)
    df = df[df["downstream_sequence"].notna()].reset_index(drop=True)
    if len(df) < n_before:
        notes.append(f"{n_before - len(df)} rows dropped: could not locate a unique "
                     "downstream sequence.")

    df["extended_wide_initial_target"] = df["wide_initial_target"].astype(str) + df["downstream_sequence"].astype(str)

    # ---- size the target region (notebook Cell 23) ----
    df["targetlength"] = (p["selftargetmaxlength"] - len(p["fiveprimeoverhang"]) -
                          p["spacerlength_with_g"] - len(p["bsmbi_spacer"]) -
                          df["extension_length"] - len(p["tevopreQ1"]) -
                          len(p["polyT"]) - p["barcodelength"] - len(p["threeprimeoverhang"]))

    # drop rows whose target region would be too short to be usable
    n_before = len(df)
    df = df[df["targetlength"] >= 20].reset_index(drop=True)
    if len(df) < n_before:
        notes.append(f"{n_before - len(df)} rows dropped: extension too long to fit "
                     f"a target within {p['selftargetmaxlength']} bp.")

    df["targetseq"] = df.apply(
        lambda x: x["extended_wide_initial_target"][5:5 + int(x["targetlength"])], axis=1)
    df["targetseq_length"] = df["targetseq"].str.len()

    # The notebook's `availabletargetminusminimumtarget` check, made enforcing:
    # if extended_wide_initial_target is shorter than 5 + targetlength, the slice
    # above silently returns a SHORT target and the construct comes out under
    # selftargetmaxlength instead of failing. That is how a predictor supplying
    # too little genomic context produces a quietly-wrong library.
    n_before = len(df)
    short = df["targetseq_length"] < df["targetlength"]
    if short.any():
        deficit = int((df.loc[short, "targetlength"] - df.loc[short, "targetseq_length"]).max())
        notes.append(
            f"{int(short.sum())} rows dropped: not enough genomic context to fill the "
            f"target region (short by up to {deficit} bp). The predictor summary must "
            f"supply a full-length Original_Sequence and a {99}-nt wide_initial_target.")
        df = df[~short].reset_index(drop=True)

    df["targetseq_revcomp"] = df["targetseq"].apply(_revcomp_dna)

    if df.empty:
        raise ValueError(
            "No constructs left to assemble. Reasons: " + "; ".join(notes)
            if notes else "No constructs left to assemble.")

    # ---- assign unique barcodes ----
    chars = ["A", "T", "G", "C"]
    barcodes = ["".join(c) for c in itertools.product(chars, repeat=p["barcodelength"])]
    random.shuffle(barcodes)
    if len(df) > len(barcodes):
        raise ValueError(f"Need {len(df)} barcodes but only {len(barcodes)} of "
                         f"length {p['barcodelength']} exist. Increase barcodelength.")
    df["barcode"] = barcodes[:len(df)]

    # ---- assemble the construct (notebook Cell 23) ----
    df["selftargetconstruct"] = df.apply(
        lambda x: (p["fiveprimeoverhang"] + x["Spacer-Sequence"] + p["bsmbi_spacer"] +
                   x["RTrevcomp"] + x["PBSrevcomp"] + p["tevopreQ1"] + p["polyT"] +
                   x["targetseq_revcomp"] + x["barcode"] + p["threeprimeoverhang"]), axis=1)
    df["selftargetconstruct_length"] = df["selftargetconstruct"].str.len()

    # ---- QC (notebook Cell 24) ----
    df["polybases"] = df["selftargetconstruct"].apply(
        lambda x: sum(_occurrences(x, b * 16) for b in "ACGT"))
    df["BsmBI_count"] = df["selftargetconstruct"].apply(
        lambda x: _occurrences(x, "CGTCTC") + _occurrences(x, "GAGACG"))
    df["fwprimcount"] = df["selftargetconstruct"].apply(
        lambda x: _occurrences(x, "AACACCG") + _occurrences(x, "CGGTGTT"))
    df["revprimcount"] = df["selftargetconstruct"].apply(
        lambda x: _occurrences(x, "GTGACTCC") + _occurrences(x, "GGAGTCAC"))

    n_before = len(df)
    df = df[(df["polybases"] == 0) &
            (df["selftargetconstruct_length"] <= p["selftargetmaxlength"])].reset_index(drop=True)
    qc_dropped = n_before - len(df)
    if qc_dropped:
        notes.append(f"{qc_dropped} rows dropped in QC (poly-base run \u226516 or over length).")

    # note stray cloning/primer sites but don't auto-drop (informational)
    n_stray_bsmbi = int((df["BsmBI_count"] > 2).sum())

    # ---- final trim: keep EXACTLY top-N per mother mutation (bystanders
    # compete within it). The heavy over-selection already happened before
    # barcoding; this just trims the survivors down to the requested top_n.
    if _top_n:
        df = _select_top(df, int(_top_n))
    else:
        # legacy fallback: global cap by library_size when no top_n is set
        if score_col:
            df = df.sort_values(score_col, ascending=False).reset_index(drop=True)
        if len(df) > p["library_size"]:
            df = df.iloc[:p["library_size"]].reset_index(drop=True)

    df = df.sort_values(["mutation", "intro_or_revert"]).reset_index(drop=True)

    # move key columns to the front
    front = ["mutation", "intro_or_revert", "bystander", "barcode",
             "selftargetconstruct", "selftargetconstruct_length"]
    front = [c for c in front if c in df.columns]
    df = df[front + [c for c in df.columns if c not in front]]

    report = {
        "n_pegRNAs_in": len(summary_df),
        "n_constructs_built": len(df),
        "n_mutations": int(df["mutation"].nunique()) if len(df) else 0,
        "top_n_per_mutation": _top_n,
        "n_bystander_constructs": int((df["bystander"] == "yes").sum()) if len(df) else 0,
        "construct_length": (int(df["selftargetconstruct_length"].iloc[0])
                             if len(df) else None),
        "all_300bp": bool((df["selftargetconstruct_length"] == p["selftargetmaxlength"]).all())
                     if len(df) else False,
        "n_with_extra_bsmbi": n_stray_bsmbi,
        "notes": notes,
    }
    return df, report


# ---------------------------------------------------------------------------
# Stage 3-QC — Quality control of assembled self-targeting constructs
# ---------------------------------------------------------------------------
#
# Verifies every construct against the rules that make a self-targeting prime-
# editing library actually work. Reports pass/fail per check per construct;
# NEVER auto-drops anything — the researcher decides. Checks are grounded in the
# biology of the construct (see each check's docstring).

def _revcomp_qc(s):
    return str(s).translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def qc_library(library_df, params=None):
    """Run QC on an assembled library. Returns (annotated_df, report).

    annotated_df is library_df plus per-check boolean columns and a 'QC_pass'
    summary column. report summarises counts and lists failing constructs with
    reasons. Nothing is removed.
    """
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    df = library_df.copy().reset_index(drop=True)

    if "selftargetconstruct" not in df.columns:
        raise ValueError("This file has no 'selftargetconstruct' column — QC runs "
                         "on an assembled library (Step 3 output).")

    checks = {}   # check_name -> boolean Series (True = pass)
    con = df["selftargetconstruct"].astype(str)

    # 1) length is within the synthesis limit (oligo pool max, e.g. Twist 300 bp).
    #    <= is allowed: shorter constructs synthesise fine; only exceeding the
    #    max is a hard failure.
    checks["length_max"] = con.str.len() <= p["selftargetmaxlength"]

    # 2) required fixed elements present (5' primer, BsmBI/scaffold spacer,
    #    tevopreQ1, polyT, 3' primer) and in the right order
    def _elements_ok(c):
        order = [p["fiveprimeoverhang"], p["bsmbi_spacer"], p["tevopreQ1"],
                 p["polyT"], p["threeprimeoverhang"]]
        pos = -1
        for el in order:
            i = c.find(el, pos + 1)
            if i == -1 or i <= pos:
                return False
            pos = i
        return True
    checks["elements_present_ordered"] = con.apply(_elements_ok)

    # 3) spacer consistency: the spacer PRIDICT designed for this row must appear
    #    in the construct right after the 5' primer (catches mis-joins / shifts)
    if "Spacer-Sequence" in df.columns:
        # 3a) the spacer must be exactly `spacerlength_with_g` nt. The 300 bp
        #     length budget in build_library uses that constant while the
        #     construct is assembled from the ACTUAL spacer, so a 21 nt spacer
        #     (a leading G prepended instead of substituted) yields a 301 bp
        #     oligo that exceeds the synthesis limit. Length is not implied by
        #     the "spacer appears in the construct" check below.
        checks["spacer_length"] = df["Spacer-Sequence"].astype(str).str.len() == \
            p["spacerlength_with_g"]

        def _spacer_ok(row):
            c = str(row["selftargetconstruct"])
            sp = str(row["Spacer-Sequence"])
            after5 = c.find(p["fiveprimeoverhang"])
            if after5 == -1:
                return False
            expected_at = after5 + len(p["fiveprimeoverhang"])
            # spacer sits immediately after the 5' primer (allowing the leading G)
            return c[expected_at:expected_at + len(sp)] == sp or \
                   c[expected_at + 1:expected_at + 1 + len(sp)] == sp
        checks["spacer_matches_pegRNA"] = df.apply(_spacer_ok, axis=1)

    # 4) self-targeting integrity: the target sequence embedded (revcomp) must be
    #    the reverse complement of this row's targetseq — i.e. the pegRNA's own
    #    target is what's recorded next to it. This is the heart of self-targeting.
    if "targetseq" in df.columns:
        def _selftarget_ok(row):
            c = str(row["selftargetconstruct"])
            trc = _revcomp_qc(str(row["targetseq"]))
            return trc in c
        checks["selftarget_integrity"] = df.apply(_selftarget_ok, axis=1)

    # 5) barcode uniqueness (duplicates make reads unassignable)
    if "barcode" in df.columns:
        dup = df["barcode"].duplicated(keep=False)
        checks["barcode_unique"] = ~dup

    # 6) no stray BsmBI/Esp3I sites beyond the two intended cloning sites
    checks["no_extra_bsmbi"] = con.apply(
        lambda c: (_occurrences(c, "CGTCTC") + _occurrences(c, "GAGACG")) <= 2)

    # 7) no stray primer-binding sites inside (would disturb PCR amplification)
    checks["no_extra_primer_sites"] = con.apply(
        lambda c: (_occurrences(c, "AACACCG") + _occurrences(c, "CGGTGTT") +
                   _occurrences(c, "GTGACTCC") + _occurrences(c, "GGAGTCAC")) <= 2)

    # 8) polyT terminator present and intact
    checks["polyT_intact"] = con.str.contains(p["polyT"], regex=False)

    # 9) no long homopolymer runs (synthesis/sequencing errors)
    checks["no_long_homopolymer"] = con.apply(
        lambda c: all(b * 16 not in c for b in "ACGT"))

    # 10) spacer / PBS / RT seamlessly reconstruct the target, and
    # 11) the RT template installs the SPECIFIC edit the mother mutation asked
    #     for — not merely "something differs in the RT window". The base change
    #     the RT writes at the edit position is mapped back to the mutation the
    #     user provided (parsed from the mutation label), matching either strand
    #     and either direction (intro installs it, revert installs the reverse).
    _seamless_cols = ["PBSlocation", "RT_mutated_location",
                      "protospacerlocation_only_initial", "PBSrevcomp",
                      "RTrevcomp", "wide_initial_target", "wide_mutated_target",
                      "deepeditposition"]
    notes = []
    if all(c in df.columns for c in _seamless_cols):
        _RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        _CB = {"A": "T", "T": "A", "G": "C", "C": "G"}

        def _rc(s):
            return str(s).translate(_RC)[::-1]

        def _loc(x):
            try:
                v = ast.literal_eval(x) if isinstance(x, str) else x
                return [int(v[0]), int(v[1])]
            except Exception:
                return None

        def _edit_from_label(name):
            """Pull the base change out of a mutation label. Handles both the
            codon form (I141N->'ATCAAC' = ATC>AAC) and the SNV form
            (TERT_124CT = C>T). Returns (from_base, to_base) for a single-base
            change, or None if it can't be read unambiguously."""
            m = re.search(r"([ACGTacgt]+)$", str(name))
            if not m:
                return None
            run = m.group(1).upper()
            if len(run) < 2 or len(run) % 2 != 0:
                return None
            k = len(run) // 2
            a, b = run[:k], run[k:]
            diffs = [(a[i], b[i]) for i in range(k) if a[i] != b[i]]
            return diffs[0] if len(diffs) == 1 else None

        def _seam_edit(row):
            """(seamless_ok, rt_edit_ok, installed_str)."""
            pl = _loc(row["protospacerlocation_only_initial"])
            pbl = _loc(row["PBSlocation"])
            rl = _loc(row["RT_mutated_location"])
            if pl is None or pbl is None or rl is None:
                return (False, False, "")
            wi = str(row["wide_initial_target"]).upper()
            wm = str(row["wide_mutated_target"]).upper()
            pbs = _rc(row["PBSrevcomp"]).upper()
            rt = _rc(row["RTrevcomp"]).upper()
            nick = pbl[1]
            rt_matches_window = (rt == wm[rl[0]:rl[1]])
            seamless = (pl[0] <= nick <= pl[1]
                        and pbl[1] == rl[0]
                        and pbs == wi[pbl[0]:pbl[1]]
                        and rt_matches_window)

            # --- closed-loop RT-carries-the-edit ---
            # Every base change the RT installs is the set of positions where the
            # initial and mutated targets differ INSIDE the RT window (rt_matches
            # _window already guarantees the RT reproduces all of them). The
            # mother mutation — parsed from the label — must be one of them. Extra
            # changes are allowed: those are the silent bystander edits.
            hi = min(rl[1], len(wi), len(wm))
            diffs = [(i, wi[i], wm[i]) for i in range(rl[0], hi) if wi[i] != wm[i]]
            label = _edit_from_label(row.get("mutation", ""))
            if label:
                lf, lt = label
                allowed = {(lf, lt), (lt, lf),
                           (_CB.get(lf, lf), _CB.get(lt, lt)),
                           (_CB.get(lt, lt), _CB.get(lf, lf))}
                mother = [(i, f, t) for (i, f, t) in diffs if (f, t) in allowed]
                mother_present = len(mother) > 0
                if mother:
                    i, f, t = mother[0]
                    installed = f"{f}>{t}@{i}"
                    extra = len(diffs) - len(mother)
                    if extra > 0:
                        installed += f" (+{extra} silent)"
                else:
                    installed = ("no " + f"{lf}>{lt}" + " in RT window; changes: "
                                 + ",".join(f"{f}>{t}@{i}" for i, f, t in diffs[:3]))
            else:
                # label unreadable — fall back to "at least one real edit present"
                mother_present = len(diffs) > 0
                installed = ",".join(f"{f}>{t}@{i}" for i, f, t in diffs[:3])

            rt_ok = bool(rt_matches_window and mother_present)
            return (bool(seamless), rt_ok, installed)

        triples = df.apply(_seam_edit, axis=1)
        checks["spacer_pbs_rt_seamless"] = triples.apply(lambda t: t[0])
        checks["rt_encodes_edit"] = triples.apply(lambda t: t[1])
        df["rt_edit_installed"] = triples.apply(lambda t: t[2])
    else:
        missing = [c for c in _seamless_cols if c not in df.columns]
        notes.append("Seamless / RT-edit checks skipped — library is missing "
                     "columns: " + ", ".join(missing) + ". Re-run Step 3 build "
                     "to include them.")

    # assemble annotated dataframe (FULL — every construct, with QC_ columns)
    for name, series in checks.items():
        df["QC_" + name] = series.values
    check_cols = ["QC_" + n for n in checks]
    df["QC_pass"] = df[check_cols].all(axis=1)

    # build report from the FULL set, so the GUI still shows EVERY construct that
    # failed any check — including ones that will be dropped from the CSV below.
    per_check = {name: int((~checks[name]).sum()) for name in checks}  # n failing
    failing = df[~df["QC_pass"]]
    fail_details = []
    for _, row in failing.head(500).iterrows():
        failed = [n.replace("QC_", "") for n in check_cols if not row[n]]
        label = row.get("mutation", "")
        name = row.get("sequence_name", label)
        fail_details.append((str(name), failed))

    # By default now: EVERY construct that fails ANY check is dropped from the
    # clean CSV. (Previously only "hard" checks like extra BsmBI were dropped.)
    # Pass params={"drop_checks": [...]} to drop only specific checks instead.
    if params is None:
        params = {}
    if "drop_checks" in params:
        drop_checks = [c for c in params["drop_checks"] if c in checks]
        if drop_checks:
            drop_mask = ~df[["QC_" + c for c in drop_checks]].all(axis=1)
        else:
            drop_mask = pd.Series(False, index=df.index)
    else:
        # default: drop anything that isn't a full QC_pass
        drop_checks = list(checks.keys())
        drop_mask = ~df["QC_pass"]

    csv_df = df[~drop_mask].reset_index(drop=True)
    n_dropped = int(drop_mask.sum())

    # Build a separate table of everything dropped, WITH a human-readable reason
    # (which checks it failed), so it can be downloaded for the record.
    dropped_df = df[drop_mask].copy().reset_index(drop=True)
    if len(dropped_df):
        reasons = []
        for _, row in dropped_df.iterrows():
            failed = [n.replace("QC_", "") for n in check_cols if not row[n]]
            reasons.append("; ".join(failed))
        # put the reason column first for readability
        dropped_df.insert(0, "QC_drop_reason", reasons)
    else:
        dropped_df["QC_drop_reason"] = pd.Series(dtype=str)

    report = {
        "n_total": len(df),
        "n_pass": int(df["QC_pass"].sum()),
        "n_fail": int((~df["QC_pass"]).sum()),
        "checks_run": list(checks.keys()),
        "failures_per_check": per_check,
        "fail_details": fail_details,
        "notes": notes,
        # CSV was filtered: constructs failing these checks were removed from it
        "drop_checks": list(drop_checks),
        "n_dropped_from_csv": n_dropped,
        "n_in_csv": len(csv_df),
    }
    return csv_df, dropped_df, report

# ---------------------------------------------------------------------------
# Stage 2-automation — Layer 1: run PRIDICT2.0 from the app (blocking)
# ---------------------------------------------------------------------------
#
# "Layer 1" automation: copy the batch CSV into PRIDICT's input/ folder, invoke
# PRIDICT in its own conda env, then read back the summary it writes. This runs
# SYNCHRONOUSLY (blocks until PRIDICT finishes) — see the note in app.py about
# browser timeouts on long runs. It does NOT replace running PRIDICT by hand;
# it just wires the same command up to a button.

def _find_conda():
    """Locate a conda executable. Prefer PATH; fall back to common install dirs
    so the app works even when launched from a shell where conda isn't on PATH."""
    exe = shutil.which("conda")
    if exe:
        return exe
    for cand in (
        "~/miniforge3/bin/conda", "~/miniconda3/bin/conda",
        "~/anaconda3/bin/conda", "/opt/miniconda3/bin/conda",
        "/opt/homebrew/Caskroom/miniforge/base/bin/conda",
    ):
        p = os.path.expanduser(cand)
        if os.path.exists(p):
            return p
    return "conda"  # last resort — will raise a clear error if missing


def _sanitize_run_name(name):
    """Make a filesystem-safe run name; blank -> timestamped default."""
    name = (name or "").strip()
    if not name:
        name = "batch_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    # keep it tame: letters, digits, dash, underscore, dot
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    if name.lower().endswith(".csv"):
        name = name[:-4]
    return name


def _archive_predictions(pred_dir, run_name):
    """Move everything currently in predictions/ into predictions/_archive/<stamp>/
    so the next run starts clean and no manual clearing is needed. Nothing is
    deleted. Returns the archive path (or None if there was nothing to move)."""
    if not os.path.isdir(pred_dir):
        os.makedirs(pred_dir, exist_ok=True)
        return None
    archive_root = os.path.join(pred_dir, "_archive")
    entries = [e for e in os.listdir(pred_dir) if e != "_archive"]
    if not entries:
        return None
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    dest = os.path.join(archive_root, f"{run_name}__{stamp}")
    os.makedirs(dest, exist_ok=True)
    for e in entries:
        shutil.move(os.path.join(pred_dir, e), os.path.join(dest, e))
    return dest


def run_pridict(batch_csv_text, pridict_home="~/PRIDICT2", conda_env="pridict2",
                cores=3, summarize="K562", summarize_number=10, timeout=None,
                script="pridict2_pegRNA_design.py", run_name=None, archive=True):
    """Run PRIDICT2.0 on a batch CSV and return its summary.

    Parameters
    ----------
    batch_csv_text : str   the batch.csv content (as produced by Stage 1)
    pridict_home   : str   path to the PRIDICT2 install (contains the script + input/)
    conda_env      : str   conda env name PRIDICT lives in
    cores          : int   parallel CPU cores
    summarize      : str   cell line used to RANK candidates (not a filter)
    summarize_number : int candidates kept per target in the summary
    timeout        : int|None  seconds before giving up (None = wait indefinitely)
    run_name       : str|None  name for this run's input file (blank -> timestamp).
                     The input is written as input/<run_name>.csv, so PRIDICT's
                     summary is named after it and past runs stay identifiable.
    archive        : bool  before running, sweep the existing predictions/ into
                     predictions/_archive/<run_name>__<timestamp>/ so the folder
                     is clean and you never have to empty it by hand.

    Returns
    -------
    dict with keys: summary_text, summary_path, n_rows, run_name, input_path,
                    archived_to, stdout, stderr
    Raises RuntimeError with a readable message on any failure.
    """
    home = os.path.expanduser(pridict_home)
    if not os.path.isdir(home):
        raise RuntimeError(
            f"PRIDICT folder not found at {home}. "
            "Set the correct path (default ~/PRIDICT2).")
    script_path = os.path.join(home, script)
    if not os.path.exists(script_path):
        raise RuntimeError(
            f"Can't find {script} inside {home}. Is this the PRIDICT2 folder?")

    input_dir = os.path.join(home, "input")
    pred_dir = os.path.join(home, "predictions")
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(pred_dir, exist_ok=True)

    run_name = _sanitize_run_name(run_name)

    # 0) archive the previous run's output, so this one starts on a clean folder
    archived_to = _archive_predictions(pred_dir, run_name) if archive else None

    # 1) drop this run's input into PRIDICT's input/ folder, named by run_name
    input_fname = f"{run_name}.csv"
    input_path = os.path.join(input_dir, input_fname)
    with open(input_path, "w") as fh:
        fh.write(batch_csv_text)

    # remember which summaries already exist, so we can spot the NEW one
    pattern = os.path.join(pred_dir, f"*_summary_{summarize}_batch_summary.csv")
    before = set(glob.glob(pattern))

    # 2) run PRIDICT in its conda env, from the PRIDICT folder
    conda = _find_conda()
    cmd = [conda, "run", "-n", conda_env,
           "python", script, "batch",
           "--input-fname", input_fname,
           "--cores", str(cores),
           "--summarize", summarize,
           "--summarize_number", str(summarize_number)]
    try:
        proc = subprocess.run(cmd, cwd=home, capture_output=True,
                              text=True, timeout=timeout)
    except FileNotFoundError:
        raise RuntimeError(
            "Could not launch conda. Start the app from a shell where conda "
            "works (one where `conda activate` succeeds), or edit _find_conda().")
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"PRIDICT still running after {timeout}s and was stopped. It may "
            "just be slow — check ~/PRIDICT2/predictions manually, or raise the "
            "timeout / switch to background (Layer 2) execution.")

    if proc.returncode != 0:
        raise RuntimeError(
            f"PRIDICT exited with code {proc.returncode}.\n"
            f"--- stderr (tail) ---\n{(proc.stderr or '')[-1500:]}")

    # 3) find the summary PRIDICT just wrote. Prefer one named after this run,
    #    then a brand-new file, then newest overall (folder is clean anyway).
    after = set(glob.glob(pattern))
    named = sorted(
        [p for p in after if os.path.basename(p).startswith(run_name)],
        key=os.path.getmtime)
    new = sorted(after - before, key=os.path.getmtime)
    hits = named or new or sorted(after, key=os.path.getmtime)
    if not hits:
        raise RuntimeError(
            "PRIDICT finished but no summary matched "
            f"'*_summary_{summarize}_batch_summary.csv' in {pred_dir}. "
            "If the summary is empty, check that the pridict2 env has "
            "pandas 2.x (pandas 3.x produces silent empty output).")
    summary_path = hits[-1]
    summary_df = pd.read_csv(summary_path)
    return {
        "summary_text": summary_df.to_csv(index=False),
        "summary_path": summary_path,
        "n_rows": int(len(summary_df)),
        "run_name": run_name,
        "input_path": input_path,
        "archived_to": archived_to,
        "stdout": proc.stdout or "",
        "stderr": proc.stderr or "",
    }
