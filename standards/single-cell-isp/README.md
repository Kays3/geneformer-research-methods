# Single-cell in silico perturbation standard

Use this workflow to assess a new single-cell or single-nucleus dataset before fine-tuning, run a matched transcriptome comparison, and decide how far the evidence supports a candidate. Use one registration per cell population, target, edit, count layer, split and principal readout. The same schema works for any number of donors and named state axes. The heart benchmark supplies an implementation check for the original model. Its genes, states, donor cutoffs and exploratory endpoints are specific to that experiment.

## Start a private study

Run these commands from the root of this methods repository. Use a new directory in your private study checkout. Python 3.9+ and its standard library suffice. Prefix the commands with `uv run --no-project` for a uv-managed interpreter.

```sh
python3 standards/single-cell-isp/isp.py init ../my-private-study/isp-run-01
# Fill in study.json, donors.tsv and cells.tsv; record scientific choices in the registration.
python3 standards/single-cell-isp/isp.py check ../my-private-study/isp-run-01
python3 standards/single-cell-isp/isp.py lock ../my-private-study/isp-run-01
# Commit the registration, lock and environment before examining test target shifts.
# Run each model adapter on the registered cells and export the table described below.
python3 standards/single-cell-isp/isp.py score ../my-private-study/isp-run-01 ../my-private-study/shifts.tsv
```

The CLI initializes and checks inputs, locks them by SHA-256 and scores an exported shift table. Model inference, reference construction, calibration and experimental validation are upstream tasks. Keep their inputs, adapter code, environment lock, checkpoint hashes, API model identifiers and receipts in the private run package. A lock detects edits relative to its stored hashes. A versioned pre-result commit establishes which design existed at scoring time. Amend via a new version and record results already seen. [Registration](../../templates/preregistration.md) · [Run record](../../templates/methods-record.md)

## Dataset suitability and scope

| Available evidence | Supported use |
|---|---|
| Biological replicate IDs, suitable comparison states and gene-resolved counts | Register a donor-held-out transcriptome experiment. Use organism- and assay-compatible checkpoints. |
| One donor, unknown donor identities, or pooled libraries without recoverable replicates | Run technical and within-sample checks. Report a descriptive result; donor transfer is unassessable. |
| Disease group is inseparable from site, protocol or batch | Seek a balanced external cohort or matched design before attributing a shift to disease. |
| No appropriate state labels or reference population | Specify a separately supported program readout or use a reference with documented transfer validation. A disease-centroid hypothesis requires justified states. |
| No relevant measured positive control | Report an uncalibrated hypothesis-generation screen. Seek context-matched perturbation data before advancing a biological claim. |
| No supported species/build, compatible gene vocabulary or usable target signal | Mark the affected arm ineligible and give its reason. Removing a model changes the scope of any cross-model claim. |
| Defined genomic allele, target gene, supported reference genome and relevant molecular track | Add the separate AlphaGenome 2 cis-expression arm below. A count matrix alone supplies no alternate allele. |

## Before fitting

1. **Fix the estimand.** Register population, intervention, start state, axes and their expected directions, background controls, gene family, exclusion rules and unit of replication. Name the measured positive control, its experimental context and the source of its expected response. Specify what constitutes an informative failure.
2. **Separate development and evaluation.** Assign every biological donor to one partition before fitting. Keep visits, aliquots, paired tissues and technical replicates together. Check class balance by site and assay. Fit references, feature selection, tuning and decision margins on development donors only. Select the checkpoint on validation donors. For grouped cross-validation, register and lock each fold separately, then combine out-of-fold donor predictions once per donor. Disclose known pretraining overlap and uncertainty.
3. **Audit representation.** Store assay and species, raw and corrected layer provenance, stable IDs, gene mapping and duplicate handling, and the denominator at each step by donor: input cells, target detected, mapped, token-retained, controls retained and actually scored. Use composite IDs if barcodes repeat across libraries. Freeze paired cell identities and exact shared control genes before target scoring. Keep the full native eligible cohort as a separate coverage sensitivity. Record both ineligible donors and reasons in the full cohort ledger. The scorer deliberately reports incomplete coverage if a registered test donor has zero paired cells.
4. **Calibrate before screening.** Assess independent no-op forwards, reproducibility of a saved edit, noncollapsed embeddings and distinguishable reference states. For fitted classifiers retain held-out per-class recall, donor balanced accuracy and a simple baseline. Use measured positive perturbations and null genes matched on detection and rank to choose model-specific margins on development data. Run matched null targets through the entire donor and cross-model rule to estimate false positives. Compare predictions of experimentally measured responses with a simple baseline fitted on the same split. Freeze the seed, comparator draw, margins, family and multiplicity procedure before opening test results. [Numerical gate patterns](../../patterns/gates.py) · [Statistical guidance](../../docs/STATISTICS.md)

## Adapter contract

| Arm | Register and verify |
|---|---|
| Geneformer V1 | V1 tokenizer and weights, gene medians, Ensembl mapping, rank order and sequence limit. Verify the native deletion on a small set; record pooling layer and special-token treatment. |
| Geneformer V2 | Pin the specific V2 checkpoint and matching vocabulary, context limit, special tokens and pooling. Embed the same reference nuclei in V2's coordinates. Score matched cell/control IDs and report its additional native coverage separately. |
| scGPT | Pin vocabulary, count layer, normalization, expression bins, input feature order, sampling seed and context length. Keep all other sampled tokens fixed for target and control edits. Record natural versus forced retention. Register omission, expression-zero or another defined edit; compare masking with physical removal if claiming equivalence. |
| AlphaGenome 2 | Register one allele per primary locus, build, REF/ALT orientation, target gene, window, cell/tissue condition, assay track, scorer and model version/fold. Verify reference sequence, replay requests and run applicable no-op/reciprocal controls. Compare ALT direction and uncertainty against an allele-harmonized molecular measurement, ideally a second study. Record missing loci and selection overlap. |

For the DNA arm, keep an expected-locus ledger and include every registered locus in the input, using NA for unavailable comparisons. Run one study, track and condition per table. The cis workflow checks the allele and gene identities against the locked ledger before using the repository locus scorer. Report locus agreement, fold sensitivity, confidence intervals and whether measurements support the predicted direction. A regulatory allele score and a cell embedding shift address different interventions and outcomes. [Cis input schema and command](../../benchmarks/cardiac-isp/README.md#run-the-scoring-fixture)

From the repository root:

```sh
python3 standards/single-cell-isp/cis.py init ../my-private-study/cis-run-01
# Complete cis.json and loci.tsv and verify the reference alleles against the declared genome.
python3 standards/single-cell-isp/cis.py lock ../my-private-study/cis-run-01
python3 standards/single-cell-isp/cis.py score ../my-private-study/cis-run-01 ../my-private-study/alleles.tsv
```

The score table uses the ledger columns plus prediction, beta and se. The molecular comparator is a separate submission for each study. Confidence intervals are normal approximations. The synthetic example also supports this workflow: after copying the example folder, run cis.py lock and cis.py score with its alleles.tsv. One of two registered loci has a measurement and its sign agrees.

## Score and interpret

The private `shifts.tsv` must contain one row per registered cell, model and axis, with columns `cell, model, axis, target_shift, control_mean_shift, controls` (tab-separated). Shifts are edited minus unedited readout on the declared axis. Supply the mean of the corresponding per-gene control shifts in that nucleus and the exact comma-separated control IDs. Retain individual control scores for audit. The scorer verifies manifest hashes, exact input coverage and control identities. It subtracts the control shift per cell, averages within donor and weights donors equally. A model passes a donor if every registered axis exceeds its signed margin. The shared result counts the same donors passing in every included transcriptome model.

Report the continuous donor estimates, target coverage and null-calibrated effect sizes for every model. Add donor-level uncertainty, leave-one-donor-out or grouped resampling where justified, and seed stability. Numerical reproducibility, successful control calibration, agreement across models, and validation against measured biology are distinct findings. The CLI's directional pass addresses only agreement. If the sample is too small to resolve the registered effect, report it as inconclusive. Use independent donors or experimental perturbations for confirmation. Keep exploratory one-axis and model-specific findings explicitly labeled.

## Smoke test

```sh
work=$(mktemp -d)
cp -R standards/single-cell-isp/example "$work/study"
python3 standards/single-cell-isp/isp.py lock "$work/study"
python3 standards/single-cell-isp/isp.py score "$work/study" "$work/study/shifts.tsv"
python3 -B -m unittest discover -s standards/single-cell-isp
```

The fixture uses synthetic values, three test cells, two test donors and two example axes. V1 and V2 pass both donors. scGPT passes X only. X is the one shared passing donor. The example cutoffs illustrate file semantics. Each real experiment establishes its own criteria.

[Geneformer perturbation protocol](https://www.nature.com/articles/s41596-026-01364-8) · [Perturbation-model baseline benchmark](https://www.nature.com/articles/s41592-025-02772-6) · [Cardiac worked benchmark](../../benchmarks/cardiac-isp/README.md)
