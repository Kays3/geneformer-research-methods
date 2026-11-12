# Approaches: what worked, what did not, when to use which

Tags `[D1]` to `[D4]` point to `docs/DATASETS.md`. "Internal record" marks a lesson from unpublished data.

## 1. Zero-shot or fine-tuned

| | Zero-shot (released weights) | Fine-tuned classifier |
|---|---|---|
| Use when | One sample or few donors; no label worth learning; a first look at whether an edit moves the embedding at all | Many donors per class; the question is about a state the classifier can separate on held-out donors |
| Goal state | A direction between centroids of the released model's embeddings, or none (magnitude only) | Class centroids of the fine-tuned model, from training donors |
| Main risk | Axes that do not resolve; effects that track token rank and sequence length | A classifier that learned site or study; collapse; a readout the fine-tune was trained on |
| Required checks | Axis resolution; size check against rank and length; rank-matched null | Donor-disjoint split; classifier gate on held-out donors; collapse checks; held-out axis and start cells |

A fine-tuned model's classification score and its usefulness for ISP are separate questions. Report both arms, released and fine-tuned, on the same readouts.

With one animal or one library, any result holds for that sample only. Say so in the first paragraph (internal record).

## 2. One split or donor cross-fitting

| | Single donor-disjoint split | k-fold donor cross-fitting |
|---|---|---|
| Donors in the statistics | The test donors only (15 of 43 in our costing `[D2]`) | Every donor once (43 of 43) |
| What it buys | Low cost | Precision: the interval on a mean within-donor shift narrows by the square root of the donor ratio (1.69 times here). Coverage of all studies in their real proportions |
| What it costs | A result that depends on which donors landed in the test set; stratify the split by study | k fine-tunes, about 2.9 times the ISP cost here, and between-model variance in the donor values |
| p-value floor | Not binding in either design at these donor counts | |

We chose k = 5. Use one recipe and seed for every fold, record the fold as a column, and downgrade a positive whose deletion median flips sign in more than one fold.

## 3. Paired within-donor or between-donor contrasts

Paired tumour against adjacent normal from the same donors removes between-donor, study, site and chemistry variation `[D2, D3]`. It does not remove between-tissue ambient RNA: tumour tissue sheds more epithelial RNA than adjacent normal tissue in the same donor. Register that limit with the claims and apply an ambient rule to every gene.

Disease as a between-donor label (the heart atlas `[D4]`) confounds any donor effect with disease. The sign test over held-out donors is the guard, and the classifier gate carries more weight.

## 4. Which null

| Null | What it controls | Use when | Caveat |
|---|---|---|---|
| Rank-matched token in the same cell, many replicates | Token rank, sequence length | Single-gene deletion, zero-shot | Compare gene with gene: a null that changes the token per cell is too narrow (I-17) |
| Matched control genes (detection and token rank), 20 per stratum | Detection, rank, the cells the gene is in | A curated panel on a fine-tuned classifier | Genes in one stratum share controls, so panel-level tests are dependent |
| Random eligible genes (100 to 200), raw shift | The generic behaviour of the model, including the delete and overexpress anti-correlation | Any claim about concordance or direction | Costs GPU time for genes of no interest; worth it |
| Matched gene sets (same size, matched rank and presence), 150 sets | Set size and composition | Set-level and program-level edits | The smallest p is 1/(N+1): 20 sets give 1/21, a screening resolution |
| Label permutation | Whether the classifier pipeline can pass on noise | Once per design | Ours gave pooled macro F1 0.419 against 0.675 on real labels, chance 0.333, and failed the gate as it should `[D4]` |

An edit always moves the embedding. The governing comparison is the observed statistic against its matched null, never against zero.

A model-based prediction should also beat a baseline that uses no model (expression level, a differential-expression rank, plain co-expression). In a separate bulk RNA-seq project that used a linear network in place of Geneformer, random genes predicted the targets about as well as the named regulators, and a co-expression baseline did as well as the perturbation model on public knockout sets (internal record; that project did not use Geneformer, and its numbers are not quoted here). We have not run the corresponding baseline comparison for Geneformer ISP.

## 5. Single genes or gene sets

With 11 held-out donors, one gene reaches a small p only if every donor agrees. A set pools genes and can reach a label where no single gene does. In one study a 7-gene set was labelled against random 7-sets (primary p 0.0082, calibrated p 0.024, donor companion p 0.00098, the floor for 11 donors) while 0 of 7 genes reached the supported level `[D4]`. The simulated false-positive rate of the label rule was 5.9%, so a single awarded label is weaker evidence than its p value suggests.

Use sets when the hypothesis is about a program. Fix an eligibility rule (a minimum number of set members present in the cell). Report the set and the single-gene readings together.

A design we wrote but did not run: perturb a gene set derived from a genome-wide association study and compare it with sets matched on length, expression and locus structure, in a cell type where the trait should act and one where it should not, with an unrelated trait and a second gene-mapping method.

## 6. All genes or a panel

An all-gene deletion screen in the protocol's style is feasible on a small number of cells: 1,100 nuclei in bf16 took 35.5 GPU-h (that run returned all zeros, see I-10) and 220 nuclei in fp32 took 27.5 GPU-h on our hardware `[D4]`. It ranks genes on pooled cells. Among 14,375 screened genes, 28 reached FDR below 0.05 in cardiomyocytes, 14 toward the goal state and 14 away, the split chance gives; in macrophages the split was 7 toward and 27 away `[D4]`. We read a screen as a source of candidates, at screen level.

A panel with matched controls costs GPU time per gene and donor, and controls dominate that cost: in our costing, dropping the larger of two panels (36 of 51 genes) saved about 21%, because the control strata remained `[D2]`. Use a panel for the donor-level test of named candidates.

## 7. Model size and precision

- **Size.** 104M and 316M classifiers reached similar held-out accuracy on one task (macro F1 0.903 and 0.898 `[D1]`), and their ISP rankings still disagreed: rank correlation 0.489 and 9 of 20 top genes shared on a targeted panel `[D1]`; no measurable agreement on 20 probe genes in heart cells, with hyperparameters tuned at 316M only `[D4]`. Similar classifier quality does not imply similar ISP readings. Read a model-size difference beside the fit-to-fit variation of one size.
- **Precision.** fp32 for the screen and its reading. bf16 is several times faster and reproduced overexpression but not deletion (I-11).
- **Run-to-run.** Two fine-tunes with the same seed, code and data agreed on 95.9% of held-out labels, with per-donor balanced accuracy moving by up to 0.045 `[D2]`; in the heart atlas a same-seed repeat agreed on 87.9% of nuclei `[D4]`. Measure this band and report any gate result that falls inside it as within run-to-run noise.
- **Cross-host.** In fp32, one probe gene gave bit-identical shifts on two hosts with different drivers `[D4]`. Compare against the same-host repeat before blaming a host.

## 8. Replication

Run the design unchanged on an independent cohort, with the differences listed in advance. The lung classifier result replicated in colon (pooled held-out balanced accuracy 0.825 and 0.904, every donor above chance), the random-gene anti-correlation replicated in sign, and the curated lung hits did not: 3 of 10 kept their sign `[D2, D3]`.

## 9. What did not work

- Reading opposite signs under deletion and overexpression as evidence of specific regulation (I-16).
- Treating a cross-validated AUC of 1.0 from an anchor-trained diagnostic as validation (I-18).
- Balancing donors to repair a class carried by one sample site (I-2).
- A no-op gate as the only numerical check (I-10).
- bf16 for deletion screens (I-11).
- A hyperparameter search inside a GPU-hour ceiling sized for something else (I-24).
- Fine-tuning without gradient clipping and without a check on the predictions (I-9).
- Reproducing the protocol's classifier score on its own atlas with our pipeline (I-5). We do not know whether the gap is ours.
- Adopting a merged criteria set: the comparison that should have justified it was not completed on its own terms.

## 10. What we have not tested

- Recovery of experimentally validated targets at a stated power. The one validated positive available to us was not scored in the criteria comparison.
- Any wet-lab validation of a candidate from these screens.
- Tissues beyond tumour-infiltrating T cells, heart cell types and one single-nucleus study.
- Library versions after our pinned commit.
