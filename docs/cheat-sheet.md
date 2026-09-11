# Cheat sheet

## Start the toolkit

```bash
conda activate snake_libanalysis
cd ~/Documents/pegRNA/pegRNA-toolkit
caffeinate -dimsu python app.py
# then open http://localhost:5050
```

## Update the toolkit

```bash
cd ~/Documents/pegRNA/pegRNA-toolkit
git pull
```

## Reset before a new NGS run

```bash
cd ~/Documents/pegRNA/ngs_analysis
rm -rf analysis/* temp_output/* fastqs/*
```

## Run PRIDICT manually

```bash
conda activate pridict2
cd ~/Documents/pegRNA/PRIDICT2
python pridict2_pegRNA_design.py batch \
  --input-fname batch.csv \
  --cores 3 \
  --summarize K562 \
  --summarize_number 10
```

## Run the pipeline manually

```bash
conda activate snake_libanalysis
cd ~/Documents/pegRNA/ngs_analysis
./run.sh --dry-run
caffeinate -i snakemake --cores 6
```

## Which environment for what

| Task | Environment | Folder |
|------|-------------|--------|
| Web app | `snake_libanalysis` | `pegRNA-toolkit/` |
| NGS pipeline | `snake_libanalysis` | `ngs_analysis/` |
| PRIDICT | `pridict2` | `PRIDICT2/` |
| OptiPrime | `optiprime` | `pegRNA-toolkit/optiprime-src/` |

## Terminal basics

| | |
|---|---|
| `cd` | change directory |
| **Control + C** | stop the running command |
| **↑** | recall the previous command |
| **Cmd + K** | clear the screen |
