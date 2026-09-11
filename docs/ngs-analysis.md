# NGS analysis (Tab 2)

The toolkit drives the whole pipeline, but the pipeline's working directory
has to be prepared first.

## Before each run

In `~/Documents/pegRNA/ngs_analysis/`:

1. Drop the raw sequencing `.gz` files into `fastqs_raw/`.
2. Copy your library construction CSV into the working directory. An example
   file is included in the folder so you can see the expected format.
3. Clear the leftovers from the previous run, or Snakemake reports
   "Nothing to be done":

```bash
cd ~/Documents/pegRNA/ngs_analysis
rm -rf analysis/* temp_output/* fastqs/*
```

## Running

In Tab 2, **Run full pipeline** does everything in one click: builds the input
sheets from your plate map, merges FASTQs, runs Snakemake, and aggregates into
the control-corrected results file.

Set these for each new run:

| Setting | What it is |
|---------|-----------|
| `library_design_file` | Your library construction CSV |
| `library_size` | The template's real row count |
| `date_prefix` | Date stamp for the output filename |
| `library_name` | Used to build sample names — must stay consistent |
| Pipeline directory | Path to your `ngs_analysis` folder |

It takes hours. Keep the tab open.

## Plate map requirements

The plate map must contain these four columns. Names are case-insensitive,
but do not put spaces inside a column name.

| Column | Meaning | Example |
|--------|---------|---------|
| `replicate` | Biological replicate number | `1`, `2` |
| `editor` | Editor or condition name | `control_A`, `PEmax_B` |
| `r1_file` | R1 FASTQ filename, exact | `HL_IIB_1_A01_..._R1_001.fastq.gz` |
| `r2_file` | R2 FASTQ filename, exact | `HL_IIB_1_A01_..._R2_001.fastq.gz` |

If any is missing, the interface stops with a "missing columns" error.

A `sample` column is treated as a per-well label, not the biological sample
name, so its contents do not matter — spaces and underscores are fine.

### How wells become samples

All rows sharing the same `(editor, replicate)` pair are one biological
sample, and their R1/R2 FASTQs are merged. The name is built as:

```
<library_name>-<editor>-rep<replicate>
```

With `library_name = ALSP-pilot`, every well with `editor = control_A` and
`replicate = 1` becomes `ALSP-pilot-control-rep1`.

## Results

The final output is `analysis/summarized/<name>.csv`. The column
`PEmax_averageedited` is each construct's true editing efficiency.

!!! warning "Save downstream outputs as CSV, not XLSX"
    Excel silently reformats sequence-like values, which corrupts them.

## Running the pipeline by hand

If you would rather not use Tab 2, **Generate sample sheets** writes the two
CSVs (`samples_<lib>_combine.csv` and `config/samples.csv`). Download them
into the working directory, then:

```bash
conda activate snake_libanalysis
cd ~/Documents/pegRNA/ngs_analysis
```

**1. Dry run.** Checks files, columns and that the DAG builds. Fast, and
catches most errors before you commit hours to a full run.

```bash
./run.sh --dry-run
```

**2. Full run.** Adjust `--cores` to your Mac.

```bash
caffeinate -i snakemake --cores 6
```

**3. Aggregate the per-sample CSVs.**

```bash
python aggregate.py \
  --analysis-dir analysis \
  --template YOUR_LIBRARY_CONSTRUCTION.csv \
  --samples config/samples.csv \
  --library-name YOUR-LIBRARY-NAME \
  --date 20260808 \
  --out analysis/summarized/20260808_YOUR-LIBRARY-NAME_df.csv
```

!!! warning "`library_name` must match"
    It has to be identical between the pipeline run and the aggregate step,
    or file matching fails.
