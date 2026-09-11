# Modifications to the upstream NGS pipeline

The Snakemake pipeline is Nicolas Mathis's work. This page records what was
changed and added on top of it, so the differences from his original are
explicit and reversible.

## Modified files

### `combine_individual/combine_techreps.py`

Merges the raw FASTQ files.

**What changed**

- Replaced `zcat {files} | gzip` with `cat {files}`
- Removed the `stderr=DEVNULL` that was swallowing error messages

**Why**

macOS ships BSD `zcat`, which does not handle these `.gz` files the way the
pipeline expects — the workaround would be to decompress, concatenate and
recompress, which is slow. `cat` concatenates gzip streams directly, with no
dependency on `zcat` at all. Suppressing stderr meant real failures looked
like silent no-ops, so that was removed too.

## New files

### `aggregate.py` (pipeline working directory)

Based on Nicolas's aggregation notebook
(`notebook_editing_analysis_focused_ALSP-HSPCpilot.ipynb`), rewritten as a
command-line script so the cells do not have to be opened and run by hand.

**What it does**: merges the per-sample CSVs, subtracts control background,
averages across replicates, and writes the final results file.

**Bug fixed along the way**: the key column in the analysis CSVs is named
`Name`, not `uniquename`. The notebook assumed the latter.

### `analysis.py` (toolkit)

Added 2026-08-22. Integrates the whole analysis workflow — merge, sample
sheets and config, snakemake, aggregate — into the toolkit web interface, so
the pipeline can be driven from Tab 2 rather than from the command line.
