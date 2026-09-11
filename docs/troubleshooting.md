# Troubleshooting

Problems that have actually come up, and what fixed them.

!!! tip "Keep this page growing"
    Every question a colleague asks is a missing entry. Add it, and answer
    with a link next time.

## Quick table

| Symptom | Fix |
|---------|-----|
| Analysis tab cannot find snakemake | You started the app from the wrong environment. Stop it, `conda activate snake_libanalysis`, restart. |
| PRIDICT result is empty | `conda activate pridict2` then `pip install "pandas>=2,<3"` |
| Pipeline `.gz` / `zcat` error | `conda activate snake_libanalysis` then `conda install -c conda-forge gzip` |
| Pipeline says "Nothing to be done" | Clear `analysis/`, `temp_output/`, `fastqs/` before the run |
| Mac sleeps during a long run | Start with `caffeinate` |
| Not sure which environment is active | The prompt shows `(pridict2)` or `(snake_libanalysis)`. Switch with `conda activate <name>`. |
| OptiPrime imports fail | A dependency was upgraded. Re-pin — see below. |
| `git pull` refuses to run | You have local edits to a tracked file. `git stash`, pull, `git stash pop`. |

## Installation

### `conda: command not found`, or the prompt has no `(base)`

You did not reopen the Terminal after `conda init zsh`. Close the window
completely, open a new one, and run `conda --version` again.

### `brew: command not found` after installing Homebrew

The two PATH lines Homebrew printed at the end did not take effect. Run them,
open a new Terminal, and check again.

### OptiPrime imports fail, or scores look wrong

Something upgraded a pinned dependency. Check the three critical versions:

```bash
conda activate optiprime
python -c "import jax; print('jax', jax.__version__)"
python -c "import numpy; print('numpy', numpy.__version__)"
python -c "import sklearn; print('sklearn', sklearn.__version__)"
```

You must see `jax 0.4.30`, `numpy 1.26.4`, `sklearn 1.0.2`. If NumPy shows
2.x:

```bash
pip install numpy==1.26.4
```

Do not "upgrade" anything in this environment. Newer versions of NumPy,
scikit-learn, LightGBM and JAX all break OptiPrime.

### The toolkit cannot find an engine

The engine paths are set in the interface, not hard-coded. Check that the
path you entered matches where you actually cloned the engine.

## Running

### macOS `zcat` fails on `.gz` files

macOS ships BSD `zcat`, which does not handle these files the way the
pipeline expects.

```bash
conda activate snake_libanalysis
conda install -c conda-forge gzip
```

### Snakemake reports "Nothing to be done"

Output from a previous run is still in place, so Snakemake thinks the work is
already finished.

```bash
cd ~/Documents/pegRNA/ngs_analysis
rm -rf analysis/* temp_output/* fastqs/*
```

## Results look wrong

### Fewer designs than expected

A run that returns fewer rows than you thought looks like success. Check the
counts in the QC report against what you expected, every time.

<!-- TODO: list the likely causes in order — restriction-site conflicts in the
     assembly architecture deserve their own entry, since a target that simply
     cannot be built should be identified quickly rather than retried. -->

### Efficiencies differ from an engine's own web server

Check that you are comparing the same design. Prediction output is not always
written in score order, so taking the first row rather than the best-scoring
row will disagree with what the web server reports for the same target.

<!-- TODO: state whether this is handled internally now, and how a user can
     verify one design against the web server themselves. -->

## Still stuck

Open an issue, and include: the input file, the options you set, the error
text, and your macOS version.

<https://github.com/fengyuxin344-crypto/pegRNA-toolkit/issues>
