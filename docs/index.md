# Self-Targeting pegRNA Toolkit

A local web interface for designing self-targeting pegRNA libraries and
analysing the sequencing data that comes back from them.

Everything runs on your own machine. No data leaves your computer.

<!-- TODO: screenshot of Tab 1 here. It does more work than any paragraph. -->

<div class="grid cards" markdown>

-   **[Install](install.md)**

    One-time setup: conda environments, prediction engines, the pipeline.
    Budget about an hour.

-   **[Start using it](starting.md)**

    Three commands to launch the interface each day.

</div>

## What it does

**Library design.** Give it annotated SnapGene files or a mutation list, and
it runs the whole design workflow — protospacer search, efficiency
prediction, ranking, construct assembly, QC — and returns order-ready
constructs.

Two prediction engines sit behind the same interface:

| Engine | Notes |
|--------|-------|
| **PRIDICT2.0** | Deep-learning prediction. Required. |
| **OptiPrime** | Alternative model with native bystander handling. Optional. |

**NGS analysis.** Drives the Snakemake pipeline end to end: merge FASTQs,
generate sample sheets, run the pipeline, aggregate into a control-corrected
results table.

## What you are installing

Four pieces, in three folders, using three conda environments. Understanding
this map first makes the rest much easier.

| Piece | Where it comes from | Environment |
|-------|--------------------|-------------|
| **The toolkit app** | this repository | `snake_libanalysis` |
| **PRIDICT2.0** | cloned from GitHub | `pridict2` |
| **OptiPrime** (optional) | cloned from GitHub | `optiprime` |
| **NGS pipeline** | Nicolas Mathis's Snakemake pipeline | `snake_libanalysis` |

!!! warning "The environments must never be mixed"
    Their dependencies genuinely conflict. The toolkit calls each engine in
    its own environment for you, so you never have to switch by hand during a
    run — but you do have to activate the right one before starting the app.

Target folder layout:

```
~/Documents/pegRNA/
├── pegRNA-toolkit/         <- this repository
│   ├── app.py
│   ├── engine.py
│   ├── analysis.py
│   ├── optiprime_integration.py
│   ├── index.html
│   └── optiprime-src/      <- cloned separately
├── PRIDICT2/               <- cloned separately
└── ngs_analysis/           <- the pipeline working directory
```

Only the toolkit lives in the repository. The engines and the pipeline are
other people's software, installed separately — see
[Installation](install.md).

## Credits

NGS pipeline: Nicolas Mathis (2024). Interface and automation: Yuxin Feng.

Prediction engines are third-party software. Please cite them if you publish
designs made with this toolkit:

- **PRIDICT2.0** — <!-- TODO: citation -->
- **OptiPrime** — <!-- TODO: citation -->
