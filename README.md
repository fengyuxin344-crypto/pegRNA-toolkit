# Self-Targeting pegRNA Toolkit

A local web interface for designing self-targeting pegRNA libraries and
analysing the sequencing data that comes back from them.

**Documentation: <https://fengyuxin344-crypto.github.io/pegRNA-toolkit/>**

<!-- TODO: screenshot here -->

## What it does

Takes annotated SnapGene files or a mutation list and returns order-ready
pegRNA constructs — protospacer search, efficiency prediction, ranking,
construct assembly and QC handled end to end. Also drives a Snakemake
pipeline for analysing the sequencing readout.

Two prediction engines behind one interface:

- **PRIDICT2.0** — deep-learning prediction (required)
- **OptiPrime** — alternative model with native bystander handling (optional)

Everything runs locally. No data leaves your computer.

## Install

```bash
mkdir -p ~/Documents/pegRNA && cd ~/Documents/pegRNA
git clone https://github.com/fengyuxin344-crypto/pegRNA-toolkit.git
```

The prediction engines and the NGS pipeline install separately —
see the [installation guide](https://fengyuxin344-crypto.github.io/pegRNA-toolkit/install/).
Budget about an hour the first time.

## Run

```bash
conda activate snake_libanalysis
cd ~/Documents/pegRNA/pegRNA-toolkit
caffeinate -dimsu python app.py
# opens http://localhost:5050
```

## Update

```bash
git pull
```

## Documentation

| | |
|---|---|
| [Installation](https://fengyuxin344-crypto.github.io/pegRNA-toolkit/install/) | One-time setup |
| [Starting the toolkit](https://fengyuxin344-crypto.github.io/pegRNA-toolkit/starting/) | Daily launch, updating |
| [Library design](https://fengyuxin344-crypto.github.io/pegRNA-toolkit/library-design/) | Tab 1 |
| [NGS analysis](https://fengyuxin344-crypto.github.io/pegRNA-toolkit/ngs-analysis/) | Tab 2 |
| [Troubleshooting](https://fengyuxin344-crypto.github.io/pegRNA-toolkit/troubleshooting/) | Known issues and fixes |

## Credits

NGS pipeline: Nicolas Mathis (2024). Interface and automation: Yuxin Feng.

Please also cite the prediction engines if you publish designs made with this
toolkit — see the documentation.

<!-- TODO: add a LICENSE file. MIT matches the surrounding ecosystem. -->
