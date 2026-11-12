# Pipeline

Nine stages from raw counts to a gated result. Each stage lists its inputs, outputs and checks. Stages 1 to 5 and 8 to 9 need no GPU. The checks are the part most often skipped and the part that saved the most GPU time.

Tags `[D1]` to `[D4]` point to `docs/DATASETS.md`.

```
1 data -> 2 tokenise -> 3 cohort and split -> 4 register -> 5 start checks
      -> 6 fine-tune (optional) -> 7 embed and perturb -> 8 statistics -> 9 record and gates
```

## 1. Prepare data

| | |
|---|---|
| Inputs | Raw integer counts with donor, tissue, disease, assay and cell-type labels |
| Outputs | A provenance table: file name, byte size, sha256, source accession, download time |
| Checks | Hash verified on download. Raw counts present (not normalised values). Donor identifier never missing. If the public file holds more donors than registered, only registered donors are read |

Raw files stay on the compute host and are not committed. Data moves host to host or through shared storage, not through a laptop.

## 2. Tokenise

| | |
|---|---|
| Inputs | Counts with `var["ensembl_id"]`, `obs["n_counts"]`, optional `obs["filter_pass"]` |
| Outputs | A tokenised dataset and its hash; token count per cell |
| Checks | Gene identifiers unique. `collapse_gene_ids` set on purpose (see I-6 in `docs/ISSUES.md`). Any count column the analysis relies on is in the matrix the tokeniser reads, not only in `obs` |

Tokenise once. Every later run checks the stored hash and reuses the same tokenised cells.

## 3. Cohort and split

| | |
|---|---|
| Inputs | Tokenised cells, donor table |
| Outputs | A split table (cell, donor, part or fold) with its sha256; the analysis-cell draw; the training pool |
| Checks | Donor-disjoint parts: all cells of a donor, and both tissues of a paired donor, sit in one part. Classes are not carried by a study or site that supplies no other class. No donor appears under two identifiers across disease groups |

Our balanced design drew a fixed number of analysis cells per donor (100) with a fixed seed, and up to three times that for training. With k-fold donor cross-fitting (k = 5), every donor is scored by a model that never saw it, and evaluation donors are never used for selection.

Cells used to build readout axes, or used as ISP start cells, are held out of training, so a fine-tune cannot learn the readout.

## 4. Register

| | |
|---|---|
| Inputs | Cohort tables, gene statistics of the held-out cells (detection, token rank), the panel |
| Outputs | The registration file; the frozen control table and null-gene draw with seeds; `rules.json`; `required_gates.json` |
| Checks | Eligibility and minimum attainable p computed on the frozen statistics. Controls drawn once. Independent sign-off on named hashes |

See `templates/preregistration.md`. Controls are drawn before any perturbation output exists: genes are grouped into strata by detection and token rank (we used 0.5 log2 detection and 5 rank-percentile points), 20 controls per stratum, panel genes and anchor genes excluded. A stratum with fewer than 20 candidates is `not_estimable`; tolerances are not widened.

## 5. Start checks and costing

| | |
|---|---|
| Inputs | Run package (a code-only checkout of the merged commit), pinned hashes, measured rates |
| Outputs | A start-check record; a projected GPU-hour figure set against the cap |
| Checks | Section 3 of `docs/CRITERIA.md` |

Measure rates on the target cells before costing a full run. A calibration of about 0.05 GPU-h replaced our inferred rates and halved the estimate `[D2]`. Fit cost per call as a fixed part plus a per-cell part, per operation; in one study deletion cost was nearly all per-cell and overexpression mostly fixed (about 82% fixed at 100 cells) `[D2]`. All-gene deletion time rose faster than sequence length: the fitted exponent was 1.44 in a pilot and 1.68 in the screen `[D4]`.

## 6. Fine-tune (optional)

| | |
|---|---|
| Inputs | Training cells, split table, released weights |
| Outputs | Best checkpoint and its sha256; per-step loss and gradient norm; per-epoch evaluation metrics; test predictions |
| Checks | Training-mode memory test. After every epoch: evaluation AUC above 0.5 and the SD of evaluation predictions above a floor. After the base embedding: embedding-collapse check. Classifier gate on held-out donors |

A recipe we ran for a two-class task on the 316M model: the library's sequence-classification head on the pretrained encoder, AdamW, learning rate 5e-5, linear warm-up then linear decay, 3 epochs, class- or stratum-weighted loss with a stated divisor, gradient clipping at max norm 1.0, a checkpoint with optimizer, scheduler and RNG state after each epoch, and the best epoch chosen on the evaluation set (internal record). We do not recommend it as it stands. Without the clipping line it collapsed in the first epoch. With clipping the first epoch trained and the two later epochs collapsed, so the kept model was a one-epoch fine-tune picked by the evaluation rule (I-9). The cause of the instability is not known.

For the T-cell classifiers on public atlases we used one fixed recipe across cohorts (1 epoch, learning rate 5e-5, batch 8, no search) `[D2, D3]`.

Write the best checkpoint to a temporary name and move it into place. Log what the logged loss is: a weighted loss can sit well above ln 2 without divergence.

## 7. Embed and perturb

| | |
|---|---|
| Inputs | Model (released or fine-tuned), held-out cells, state definitions, gene list, controls, null genes |
| Outputs | State embeddings (centroids) per model; raw per-cell shifts per gene, operation and donor; a completion marker per call with its own timer; an output check (count, zeros, NaN, median absolute shift) |
| Checks | No-op gate and resolvability gate before the screen; spot checks during it; the same resolvability check after it. Delete and overexpress scored on the same cells |

Operations:

- **Delete** removes the gene's token from the rank-value encoding. It can act only in cells that hold the token.
- **Overexpress** moves the gene's token to the front of the encoding.

Readouts:

- **Shift toward a goal:** cosine(perturbed cell, goal centroid) minus cosine(original cell, goal centroid). Goal centroids come from training donors only.
- **Magnitude with no goal state:** 1 minus the cosine between the embedding before and after the edit.
- **Direction on an axis:** the projection on the unit vector between two state centroids. Build the axis on embeddings with the perturbed gene's token removed, and cross-fit it so a cell is projected only on an axis built without it. An axis must pass a resolution check before it is read.

State embeddings are specific to the model. The 104M and 316M models differ in embedding width, so centroids are recomputed for each.

Run the screen and its reading in fp32 unless a precision check on the same model and operation passed (I-11). Run Geneformer's dataset mapping with one process once CUDA is initialised (I-7).

The raw output holds per-cell scalar shifts by goal state, not embeddings. An analysis that needs embeddings must plan its own pass.

## 8. Statistics

| | |
|---|---|
| Inputs | Raw per-cell shifts, control table, null genes, `rules.json`, gate results |
| Outputs | One row per registered gene with its status; panel- or set-level results; sensitivity lines |
| Checks | The launcher verifies `required_gates.json`. Analysis code was run on synthetic data before unblinding. Rules are compared byte for byte with the registered file at launch |

See `docs/STATISTICS.md`.

## 9. Record and gates

| | |
|---|---|
| Inputs | Result tables, run logs, ledger |
| Outputs | A methods record (`templates/methods-record.md`), a run record, committed result tables and their hashes, a GPU-hour ledger entry |
| Checks | Three independent reviews before merge: science (design and reading), numbers (every figure recomputed from the tables), privacy (no names, paths, hosts or identifiers). See `templates/gate-checklist.md` |

Large per-cell outputs stay on the compute host. The repository holds tables, configs, ledgers and hashes, and says what cannot be recomputed from the repository alone.

## Driver conventions

- One resumable shell script. Each step leaves a done file; a restart skips finished steps, chunks and units.
- Start detached with output to a log, a pid file and an exit-code file written by a trap, so a dropped connection cannot kill or duplicate the run.
- Host folders and interpreters come from required environment variables. No absolute path is committed.
- One environment per host from a committed lock file. A virtual environment is not copied between hosts. Scripts stop when the library commit or weights hash differ from the pins.
- One heavy job per host. A GPU-sharing check runs between chunks.
- A recurring job on the same host is written into the driver as a blackout window.
- For a chain of multi-hour stages, one driver with a completion marker per stage and one poller for the whole chain.

`patterns/` holds generic examples of the guard, the gates and the control draw.
