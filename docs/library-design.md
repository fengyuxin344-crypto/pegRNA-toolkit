# Library design (Tab 1)

[Start the toolkit](starting.md), open <http://localhost:5050>, and use Tab 1.

## Providing mutations

Two ways:

=== "SnapGene files"

    Drag in annotated `.dna` files. Most reliable, and multiple files are
    supported.

    The CDS annotation is what lets the toolkit place silent bystander edits
    in the correct reading frame. Without it, a substitution that is
    synonymous in the protein can still land outside the coding sequence,
    where "silent" is not meaningful.

=== "Mutation list"

    A CSV with a `gDNA_mutation` column. Add a `gene` column if you have more
    than one gene.

## Bystander mode

<!-- TODO: screenshot of the bystander options, and a sentence per mode
     saying when to pick it. This is the least self-explanatory control in
     the interface. -->

## Running

Use the one-click run. When it finishes, click **Run quality check**, then
**Download QC report** at the bottom of the tab — scroll down, it sits below
the one-click run button.

<!-- TODO: screenshot of a completed run, and of a QC report -->

## Reading the QC report

The report lists every construct that was removed and why.

!!! note "Failed constructs are dropped from the CSV, not hidden"
    Constructs that fail a hard filter do not appear in the downloadable CSV,
    but they are still shown in the on-screen report. Check the counts against
    what you expected — a run that returns fewer rows than you thought looks
    like success.

<!-- TODO: a table of the QC checks: what each one means, and whether failing
     it removes the construct. -->

## Runtime

Prediction dominates. Both engines enumerate many PBS × RTT combinations per
target and score each with a neural network, so cost per input mutation is
higher than it looks.

<!-- TODO: measured timings on the lab's own hardware. Anyone planning a large
     library needs this before they start, not after. -->

| Targets | Engine | Approx. runtime |
|---------|--------|-----------------|
| 10 | | |
| 100 | | |
| 1,000 | | |
