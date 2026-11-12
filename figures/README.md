# Figures and reproducibility

These diagrams explain the methods in this repository. They are schematics, not measurements. The final panel is generated from synthetic values and does not report a study result. Blue and orange identify contrasting roles, with text and shape used as additional cues.

| Figure | What it shows | Source for the content | Files |
|---|---|---|---|
| 01 | Counts to tokenisation, embeddings, perturbations, readout and review | [Pipeline](../docs/PIPELINE.md): six steps that summarise its nine stages | [SVG](01-pipeline.svg) · [PNG](01-pipeline.png) |
| 02 | Start-to-goal axis, edit displacement and matched null | [Pipeline, stage 7](../docs/PIPELINE.md#7-embed-and-perturb); [Statistics, section 2](../docs/STATISTICS.md#2-compare-with-a-matched-null-never-with-zero) | [SVG](02-axis-null.svg) · [PNG](02-axis-null.png) |
| 03 | Released and fine-tuned comparison, donor split and gates | [Approaches](../docs/APPROACHES.md), sections 1–2; [Criteria](../docs/CRITERIA.md), sections 1 and 3 | [SVG](03-comparison-design.svg) · [PNG](03-comparison-design.png) |
| 04 | Collapse symptoms, diagnostics and limits on interpretation | [Issues, I-9](../docs/ISSUES.md#i-9-a-fine-tune-collapses-to-a-constant-predictor) | [SVG](04-collapse.svg) · [PNG](04-collapse.png) |
| 05 | Registration, frozen inputs, run and independent reviews | [Criteria](../docs/CRITERIA.md), sections 2–4; [Pipeline, stage 9](../docs/PIPELINE.md#9-record-and-gates) | [SVG](05-review-flow.svg) · [PNG](05-review-flow.png) |
| 06 | Route from two common symptoms to checks and a bounded rerun | [Issues](../docs/ISSUES.md), especially I-9 and I-11 | [SVG](06-troubleshooting.svg) · [PNG](06-troubleshooting.png) |
| 07 | Invented matched-null distribution and an example median | [Statistics, section 2](../docs/STATISTICS.md#2-compare-with-a-matched-null-never-with-zero); synthetic values generated in [draw.py](src/draw.py) | [SVG](07-synthetic-readout.svg) · [PNG](07-synthetic-readout.png) |

Every file is drawn by [figures/src/draw.py](src/draw.py). From the repository root, with Python, NumPy and Matplotlib installed, regenerate all seven offline:

```sh
python figures/src/draw.py
```

The synthetic panel uses a fixed random seed. No dataset download, Geneformer installation or GPU is required. The diagrams summarise the linked text; consult the documents for the full decision rules and limitations.
