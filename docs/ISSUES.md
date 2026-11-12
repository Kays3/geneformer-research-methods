# Failure modes

Each entry: what we saw, the root cause, how to detect it, and the fix. Entries marked "internal record" come from unpublished data and carry no numbers. Tags `[D1]` to `[D4]` point to `docs/DATASETS.md`.

Several of these runs finished with exit code 0. A clean exit says the code ran, not that the output means anything.

## Data and tokenising

### I-1. The weights file is a 135-byte pointer
- **Symptom.** A model loads with random weights, or fails to load, on one host only.
- **Root cause.** The model folder came from a clone without the large-file download. It holds a pointer file in place of the weights.
- **Detect.** Check the weights file's sha256 against the pinned value at every start.
- **Fix.** Pin the folder and the hash in the driver; refuse to run on a mismatch. Also write the hash into the run config: one early run did not, and its model identity could be shown later only by bit-equal reproduction (internal record).

### I-2. A class is carried by one sample site or one study
- **Symptom.** A classifier separates disease classes with high accuracy.
- **Root cause.** One class came mostly from a different sample site (in the small-cell atlas, 74.9% of the held-out cells of one disease class came from a single sample, a pleural effusion `[D1]`) or from studies that supply no other class. The classifier can learn site, chemistry or study.
- **Detect.** Before training, tabulate class by study, site, assay and donor.
- **Fix.** Select datasets with the structural rule in `docs/CRITERIA.md`, section 2. Balancing donors does not repair this.

### I-3. One donor under several identifiers
- **Symptom.** Donor counts do not add up across disease groups.
- **Root cause.** A patient who contributed to two disease groups carried different identifiers.
- **Detect.** Audit the donor table before the split.
- **Fix.** A cross-disease donor guard: such a patient goes to one part of the split, whole.

### I-4. A marker-based cell call does not validate
- **Symptom.** Independent markers are only weakly enriched in the called cells. In heart nuclei the odds ratio was 2.10 (95% CI 1.86 to 2.36) against a registered bar of 3 `[D4]`.
- **Root cause.** Sparse capture in single-nucleus data leaves many true cells of the type uncalled, which dilutes the contrast. The call may also be wrong; the data could not tell.
- **Fix.** Register the validation and the fallback in advance. Ours fell back to the authors' broader class, and every later result used that name.

### I-5. The published example holds different cells from the atlas
- **Symptom.** A reproduction of a published classifier misses the reported score: macro F1 0.644 (95% donor-bootstrap interval 0.332 to 0.741, 5 test donors) against a reported 0.85 `[D4]`.
- **Root cause.** Not established. The example dataset held 8.5% fewer cells of the class than the atlas, by a share that differed by disease, and carried no cell identifier, so its cells could not be selected. The test set also had 5 donors with one donor carrying 46% of the nuclei. A run on the unfiltered class gave 0.764 (0.456 to 0.812).
- **Fix.** Log it as a deviation, label every later value "post-deviation", and do not read later differences as a property of the method under test.

### I-6. The tokeniser sums features that share an identifier
- **Symptom.** Two features that should stay apart behave as one token.
- **Root cause.** `collapse_gene_ids` defaults to True and sums such features without a warning.
- **Detect.** Set it to False once; the tokeniser then raises on non-unique identifiers.
- **Fix.** Decide the mapping on purpose and record it. Related tokeniser facts: genes outside the token dictionary are dropped; each cell is truncated to the model input length, so low-ranked genes of a deep cell are never seen; keys of the custom attribute dictionary must exist in `obs`, or the labels do not reach the tokenised dataset; outputs tokenised with one dictionary family cannot be read with another.

## Library behaviour

### I-7. A crash on the first perturbation call
- **Root cause.** The library's dataset mapping uses a fork-based multiprocessing package. Forking after CUDA is initialised crashes. A spawn setting in your own runner does not reach it.
- **Fix.** Run with one process.

### I-8. Overexpression and deletion cover different cells
- **Symptom.** Every overexpression output holds one value per start cell, while deletion outputs hold one value per token-positive cell. In one study 14,738 of 15,179 overexpression calls held 100 values where deletion held k `[D2]`.
- **Root cause.** With a gene list, the library skips the token filter for overexpression and inserts the gene into every start-state cell. It then sorts cells by length, so output order is length-descending, not dataset order.
- **Detect.** A constant overexpression n across genes while deletion n varies.
- **Fix.** Score both operations on the same token-positive cells. Either pre-filter the input per gene, or reconstruct the subset by replaying the library's filter and sort. Check the reconstruction by count and by an identity run on pairs chosen by a rule frozen in advance (count checks test cardinality only). 80.5% of our gene-donor pairs had a length tie between a token-positive and a token-negative cell, so sort stability matters; test it directly instead of refusing ties `[D2]`.
- **Lesson.** Our synthetic test fixture encoded what we believed the library did. Check fixtures against the real library.

### I-9. A fine-tune collapses to a constant predictor
- **Symptom.** Every test cell gets the same score, evaluation AUC falls to 0.5 or below, the training loss sits near ln 2, and all cell embeddings become the same vector. The run exits 0 (internal record).
- **Root cause.** Not known. The first run had no gradient clipping and collapsed in its first epoch. A re-run that changed one thing, clipping at max norm 1.0, trained through the first epoch and then collapsed in the second and third. Logged gradient norms before clipping were several times the clip value in every epoch, so the clip was active throughout. Clipping delayed the collapse; it did not prevent it, and the re-run does not support missing clipping as the cause. Candidates we have not tested: the learning rate after warm-up, the weighted loss, and the newly initialised classification head.
- **Detect.** After every epoch, check whether evaluation AUC is at or below 0.5 or the SD of evaluation predictions is below 1e-3. Our stop rule looked at the first epoch only, so the later collapse did not stop the run; only the keep-the-best-epoch rule kept the collapsed epochs out of the saved model. After the base embedding, compute the mean pairwise cosine over a seeded, stratified draw of 200 cells; at or above 0.99999 the arm is void. Set that limit far from a healthy model's value: a first limit of 0.999 was too close.
- **Fix.** None established. What protected the result: the best checkpoint chosen on the evaluation set and the embedding-collapse check on the kept model. The per-epoch log of evaluation AUC and prediction SD showed the collapse; a stop on it after every epoch is what we recommend. Log the per-step loss and the pre-clip gradient norm. Report the kept model as what it is (here, a one-epoch fine-tune) and report the collapsed epochs in full.
- **Related.** In a hyperparameter search, one trial at a learning rate of 2.6e-4 predicted one class for every evaluation nucleus `[D4]`.
- **Warning.** A collapsed model still returns ISP numbers: shifts near the float limit and projections on unresolved axes. They are void, not null results.

### I-10. All shifts are exactly zero
- **Symptom.** A full screen returns zeros. Ours cost tens of GPU-hours before anyone looked `[D4]`.
- **Root cause.** The model's embeddings were nearly constant (part of the multi-task model started from random initialisation), and bf16 quantised the tiny differences to zero. In fp32 the same shifts were at rounding-noise level, so fp32 would not have rescued the reading.
- **Detect.** The resolvability gate before the screen, and an output check after the first calls. The no-op gate passes in this state.
- **Fix.** The gate. Count exact zeros by event: one nucleus with three goal states is one event, not three. A true zero is expected when the token is already at the top rank.

### I-11. Reduced precision changes deletion results more than overexpression
- **Symptom.** bf16 or int8 against fp32, 20 probe genes: overexpression medians agreed (Pearson r 0.997 to 0.999) and deletion medians did not reach the 0.99 bar (0.885 to 0.977) `[D4]`.
- **Root cause.** Deletion shifts are small, and bf16 cosine values near 1 move in steps of roughly 2^-8.
- **Fix.** Run the screen and its reading in fp32, or report agreement per operation on the same model. Report the absolute difference of medians beside their spread, not Pearson r alone. A "top-20 overlap" clause with 20 probe genes is always met and decides nothing.

### I-12. Training in bf16 produces exact logit ties
- **Symptom.** The library's evaluation crashes with a mixed-type error.
- **Root cause.** The vote function returns the string "tie" when two class logits are equal, which turns the prediction array into strings. Ties are rare in fp32 and plausible in bf16. We saw at least one in an evaluation split and 0 in 9,376 test predictions `[D1]`.
- **Fix.** Patch the vote function from outside the library file: break ties by the first index and count them. The checkpoint was saved before the crash, so no retraining was needed.

### I-13. `bf16=True` is not a bf16 model
- The trainer flag keeps fp32 weights and computes under bf16 autocast; the saved checkpoint is fp32. A model cast to bf16 is a different regime. Say which one ran.

### I-14. "Inference only" does not quantise a plain classifier
- In the public library source we read, the perturber quantises only for the model types named as quantised. Check the loader before you report a "quantised against full precision" comparison.

## Design and statistics

### I-15. A goal axis that encodes the perturbed gene
- **Symptom.** Deleting the gene moves cells along the axis, strongly.
- **Root cause.** The axis was built on embeddings that still held the gene's token, so the move is by construction (internal record).
- **Fix.** Remove the token from every cell before building axes.

### I-16. Opposite signs in the two operations look like specificity
- **Symptom.** A gene moves cells toward the goal when deleted and away when overexpressed.
- **Root cause.** Random genes do this too. Across random eligible genes the rank correlation of deletion and overexpression shifts was -0.593 in lung T cells and -0.245 in colon T cells `[D2, D3]`.
- **Fix.** Read concordance against a random-gene null, with matched controls. See `docs/STATISTICS.md`, section 5.

### I-17. A null that is too narrow
- **Symptom.** Any named gene clears the 95th percentile of the null by a wide margin.
- **Root cause.** Each null replicate deleted a different random token per cell, so replicate medians averaged out gene-to-gene spread (internal record).
- **Fix.** A per-gene null: the same random gene deleted across cells.

### I-18. A circular ambient-RNA diagnostic
- **Symptom.** A cross-validated AUC of 1.0.
- **Root cause.** The diagnostic was trained and validated on its own anchor gene lists.
- **Fix.** Leave-one-anchor-out scoring with an AUC bar. If it fails, every non-anchor gene is undetermined and cannot be a plain positive.
- **Related.** In a small panel, detection and ambient risk were almost perfectly rank-correlated, so an ambient test could not be told from a detection effect `[D1]`.

### I-19. Cell-weighted and donor-weighted scores disagree
- **Symptom.** A pooled score orders two disease groups one way and the donor-level score the other way `[D1]`.
- **Root cause.** One donor with many cells carries the pooled value.
- **Fix.** The donor is the unit. Also check denominators: two pipelines that normalised by different totals (a stored count column against the current matrix row sum) differed by a near-constant ratio.

### I-20. A critical value from a sorted null
- **Root cause.** A critical value taken as a quantile index into a sorted null differed from the value in the set of attainable p values (internal record). A library "exact" rank test that ignores ties also gave a different p from full enumeration.
- **Fix.** Enumerate. See `patterns/attainable_p.py`.

### I-21. Calibration checks that cannot fail
- **Root cause.** A rule that fails only when the lower end of the interval passes alpha tolerates about twice the nominal hit count with 100 null genes.
- **Fix.** Report the share of null genes below alpha with its exact interval, and size the null set so the interval can separate the nominal rate from a clearly inflated one.
- **Related.** A "simulated false-positive rate" that is an algebraic identity of the test (for a rank p among n, only the top k can ever pass) confirms nothing when reproduced. A calibration figure must be able to come out otherwise (internal record). A gate built on a measured noise floor degenerates when fp32 inference is bit-identical on repeat, because the floor is zero; give such a gate an explicit epsilon.

## Compute and operations

### I-22. Cost estimated from wall clock
- **Symptom.** A planning rate about three times too high: 6.59 and 6.95 s per call from wall clock against 1.846 and 4.859 s from per-call timers `[D2]`.
- **Root cause.** Wall-clock totals include start and stop overhead, log flushes and checksum time.
- **Fix.** Write a timer around the perturbation call into each completion marker and cost from that. State what the timer includes (model load) and excludes (dataset build).

### I-23. Predictions carried across cohorts
- A second cohort's ISP cost about 2.5 times its prediction, probably because its null genes were more highly detected, so each call had more token-positive cells `[D3]`. The cause was not tested. Recost on the target cells.

### I-24. A hyperparameter search that does not fit
- The search framework ran at about 0.4 to 0.5 times the rate of the same fine-tune outside it; host and precision were ruled out by measurement and the cause was not found. 4 of 30 registered trials fit under the ceiling `[D4]`. Best of four on five donors is optimistic, and the top two trials differed by 0.00006 in the selection metric.

### I-25. The fp32 to bf16 cost ratio is not one number
- Across our runs fp32 cost roughly 1.5 to 5 times as much as bf16, depending on the arm and the batch shape. This is a compute lesson from internal run notes across several runs, not a result tied to one dataset. Reprice each arm from its own gate run.

### I-26. Small units waste the GPU
- Chunks of same-length sequences gave full batches and several times the rate of many small units (internal record).

### I-27. Operational traps
- gzip files of identical tables hash differently across builds, because the header stores the write time. Pin the hash of the uncompressed content.
- A reproduction check at 1e-6 fails on batch-shape noise when the re-embedding uses a different batch size. Use a tolerance relative to the effect.
- A driver stopped after the GPU work because a numpy bool is not JSON-serialisable. Cast before writing; keep the units already written.
- A budget monitor that found its process with a pattern match on the command line matched the wrong process. Target the run's process group.
- A launch over ssh hung because the whole command list was backgrounded. Launch the detached job as its own statement with stdin closed.
- A crash in evaluation after training is not a training failure. Check the saved checkpoint before spending GPU time again.
- A partial clone lacks old blobs, so "the bytes that ran" cannot be read from it.
- In zsh, `"$C:path"` is parsed as a modifier; write `"${C}:path"`.
- A gate that finished was not polled, and the GPU sat idle for 91 minutes before the next stage started `[D4]`. Chain stages in the driver.
