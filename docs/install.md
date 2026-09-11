# Installation

One-time setup. Budget about an hour the first time — much of it is downloads
that run on their own.

Everything below is typed into the **Terminal** app.

| Part | What |
|------|------|
| 1 | Install the basics: Command Line Tools, Homebrew, Miniforge |
| 2 | Lay out your folders and get the toolkit |
| 3a | Install PRIDICT2.0 (environment: `pridict2`) |
| 3b | Install OptiPrime — optional (environment: `optiprime`) |
| 4 | Set up the NGS pipeline and toolkit (environment: `snake_libanalysis`) |

Then see [Starting the toolkit](starting.md).

---

## Part 1 — Install the basics

### 1a. Apple Command Line Tools

This is what gives you `git`.

```bash
xcode-select --install
```

A window pops up. Click **Install** and wait. If it says already installed,
move on.

### 1b. Homebrew

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

It asks for your password and takes a few minutes. At the end it prints two
lines for you to run, which add Homebrew to your PATH. **Run the lines it
prints** rather than copying from here — the path differs between Apple
Silicon and Intel Macs. On Apple Silicon they are:

```bash
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)"
```

Check:

```bash
brew --version
```

A version number (e.g. `Homebrew 4.x.x`) means it worked. If you get
`command not found`, the two PATH lines did not take effect — run them, then
open a new Terminal window and check again.

### 1c. Miniforge (conda for Mac)

```bash
brew install --cask miniforge
conda init zsh
```

**Close the Terminal window and open a new one.** This is required — it
reloads your shell settings. The new prompt should start with `(base)`.

```bash
conda --version
```

A version number plus `(base)` at the start of your prompt means conda is
ready. If you see `command not found: conda` or no `(base)`, you almost
certainly did not reopen the Terminal. Close it completely, open a new one,
and check again.

---

## Part 2 — Folders and the toolkit

One top folder holds everything.

```bash
mkdir -p ~/Documents/pegRNA
cd ~/Documents/pegRNA
git clone https://github.com/fengyuxin344-crypto/pegRNA-toolkit.git
```

You now have `~/Documents/pegRNA/pegRNA-toolkit/`. The remaining parts add the
engines and the pipeline beside it.

Target layout when you are done:

```
~/Documents/pegRNA/
├── pegRNA-toolkit/         <- the toolkit app (this repo)
│   ├── app.py
│   ├── engine.py
│   ├── analysis.py
│   ├── optiprime_integration.py
│   ├── index.html
│   └── optiprime-src/      <- Part 3b
├── PRIDICT2/               <- Part 3a
└── ngs_analysis/           <- Part 4, Nicolas's pipeline working directory
    ├── env.yaml            <- conda recipe, builds snake_libanalysis
    ├── Snakefile
    ├── run.sh              <- convenience runner (./run.sh --dry-run)
    ├── aggregate.py        <- per-sample CSVs -> final results file
    ├── combine_techreps_cat.py
    ├── config/
    │   ├── self_targeting_config.csv
    │   └── samples.csv
    ├── scripts/
    ├── workflow/
    │   ├── rules/
    │   └── envs/
    ├── combine_individual/
    ├── fastqs_raw/         <- you drop raw sequencing .gz files here
    ├── fastqs/             <- merged, created by the run
    ├── analysis/           <- per-sample results + summarized/
    ├── temp_output/        <- cutadapt intermediates, safe to delete
    └── logs/
```

---

## Part 3a — PRIDICT2.0 (environment: `pridict2`)

```bash
cd ~/Documents/pegRNA
git clone https://github.com/uzh-dqbm-cmi/PRIDICT2.git
cd PRIDICT2

conda env create -f pridict2_repo.yml
conda activate pridict2

pip install torch==2.0.1 --index-url https://download.pytorch.org/whl/cpu
```

Apple Silicon (M1/M2/M3/M4) also needs TensorFlow:

```bash
conda install -c conda-forge tensorflow
```

Test, still in `pridict2` and inside `~/Documents/pegRNA/PRIDICT2`:

```bash
python -c "import torch; print('torch ok', torch.__version__)"
```

`torch ok 2.0.1` means PRIDICT is ready.

PRIDICT needs pandas 2.x. If a prediction ever comes back empty:

```bash
conda activate pridict2
pip install "pandas>=2,<3"
```

---

## Part 3b — OptiPrime (optional, environment: `optiprime`)

OptiPrime is an alternative prediction engine. It uses JAX rather than
PyTorch, so it gets its own environment, completely separate from `pridict2`.
Skip this part if you only want PRIDICT.

**Read this first.** OptiPrime's dependencies are fragile: newer versions of
NumPy, scikit-learn, LightGBM and JAX all break it. The commands below pin the
exact versions that work together. Install them **in this order** and do not
upgrade anything afterwards. Python must be 3.10 — scikit-learn 1.0.2, which
OptiPrime's `rs3` dependency needs, does not support 3.11.

**1. Get the code**

```bash
cd ~/Documents/pegRNA/pegRNA-toolkit
git clone https://github.com/alvin-hsu/optiprime-src.git
```

You can put `optiprime-src` anywhere — the toolkit asks for its path in the
interface. Keeping it inside the toolkit folder matches the default.

**2. Create the environment**

```bash
conda create -y -n optiprime python=3.10
conda activate optiprime
```

**3. The hard packages, via conda**

These fail to build or clash on version if installed with pip. Conda gives you
pre-compiled builds at the versions `rs3` needs.

```bash
conda install -y -c conda-forge lightgbm=3.3.5 pyarrow scikit-learn=1.0.2 numpy=1.26.4
```

**4. JAX, pinned so it does not drag in a newer NumPy**

```bash
pip install jax==0.4.30 jaxlib==0.4.30
pip install numpy==1.26.4
```

**5. The remaining packages**

```bash
pip install scipy matplotlib "pandas<2.2" optax flax networkx h5py biopython tqdm ViennaRNA requests
```

**6. rs3 and sglearn, without their dependencies**

`--no-deps` stops them pulling in a newer NumPy or an older LightGBM, which
would undo the versions set above.

```bash
pip install --no-deps rs3 sglearn seqfold
```

**7. One more JAX helper**

```bash
pip install chex==0.1.86 numpy==1.26.4
```

**8. Verify the three critical versions**

```bash
python -c "import jax; print('jax', jax.__version__)"
python -c "import numpy; print('numpy', numpy.__version__)"
python -c "import sklearn; print('sklearn', sklearn.__version__)"
```

You must see `jax 0.4.30`, `numpy 1.26.4`, `sklearn 1.0.2`. If NumPy shows
2.x, run `pip install numpy==1.26.4` again and re-check.

**9. Smoke test**

```bash
cd ~/Documents/pegRNA/pegRNA-toolkit/optiprime-src
python DESIGN_PE.py smoke
```

No error and a clean exit means the imports work.

**10. Run the official example**

```bash
python DESIGN_PE.py run \
  --run_data COL7A1_R185X.json \
  --graph_rx graphs/pe_model.rx \
  --weight_dirs weights/* \
  --out_path optiprime_test_out/
```

The first run is slow because JAX compiles. Success is a folder
`optiprime_test_out/<job>/` containing `full_results.txt.gz`. Check it:

```bash
cd optiprime_test_out/*/ && gunzip -k full_results.txt.gz && head -5 full_results.txt
```

A table with `pegRNA_name`, ..., `OptiPrime_score` means OptiPrime is ready.

---

## Part 4 — NGS pipeline and toolkit (environment: `snake_libanalysis`)

The pipeline ships an `env.yaml` that already includes snakemake, so build the
environment from that file and then add the extras the toolkit and macOS need.

<!-- TODO: say how a new user obtains ngs_analysis/ — a repo URL, a shared
     drive, or "ask Yuxin". Right now this is the one step someone cannot
     complete on their own. -->

```bash
cd ~/Documents/pegRNA/ngs_analysis

conda env create -f env.yaml
conda activate snake_libanalysis

conda install -y -c conda-forge flask
pip install requests biopython openpyxl

# macOS quirk: BSD gzip cannot read the pipeline's .gz files - install conda's
conda install -y -c conda-forge gzip
```

Test:

```bash
python -c "import flask, pandas, requests, Bio, openpyxl; print('toolkit deps ok')"
snakemake --version
cutadapt --version
```

Three clean lines with no errors means you are ready.

**Why `env.yaml` rather than a fresh environment?** It pins the exact
snakemake and cutadapt versions the pipeline expects. Building from it avoids
version surprises later.

---

---

Next: [Starting the toolkit](starting.md).
