# Public datasets behind the numbers

The [illustrated guide](../figures/README.md) uses schematics and a clearly labelled synthetic example; it does not plot these datasets.

Every number in this repository comes from one of the datasets below. Each number in the other files carries a tag such as `[D2]` that points here. The analysis code for these studies is not public. Public status was taken from our download records; check the hosting site before you rely on an accession.

| Tag | Dataset | Accession | Cells used | Models |
|---|---|---|---|---|
| D1 | Small-cell lung cancer atlas, T-cell object (lung adenocarcinoma and normal lung included) | CELLxGENE collection `62e8f058-9c37-48bc-9200-e767f318a8ec` | T cells, three disease classes | Geneformer V2 104M (fp32 fine-tune); V2 316M (bf16 fine-tune) as a benchmark |
| D2 | Extended single-cell lung cancer atlas | CELLxGENE dataset version `33165751-ae33-4a65-94ee-52d9bc38f97e` | T cells of 43 donors with both tumour and adjacent-normal samples; 100 analysis cells per donor and tissue | Geneformer V2 316M, bf16, 5-fold donor cross-fitting |
| D3 | Colorectal cancer atlas | GEO `GSE178341` (processed counts; the controlled-access raw reads were not used) | T cells of 19 donors, same design as D2 | Geneformer V2 316M, bf16, 5-fold donor cross-fitting |
| D4 | Heart atlas of non-failing, hypertrophic and dilated cardiomyopathy donors, the example dataset of the Geneformer protocol | Human Cell Atlas project `85c0d6fa-f117-4d76-b01a-5d5e8f5f9188`; also the `human_dcm_hcm_nf.dataset` in the Genecorpus-30M collection on Hugging Face | Cardiomyocytes, macrophages, lymphocytes | Geneformer V2 316M (fine-tuned, single-task fold models and one multi-task model); V2 104M for one comparison |

## Model and library

- Checkpoints: the public `Geneformer-V2-104M` and `Geneformer-V2-316M` weights from the Geneformer model repository on Hugging Face. Record the sha256 of the weights file you load. A clone made without the large-file download holds a pointer file of about 135 bytes in place of the weights (see `docs/ISSUES.md`, I-1).
- Library: one pinned commit of the public Geneformer repository (December 2025). Statements about library behaviour in this repository were read at that commit.
- Hardware: single-GPU workstations. GPU-hour figures quoted here are for that hardware and are useful as ratios, not as absolute costs.

## What is not here

No dataset from an unpublished study is named or described. Lessons from such studies appear as "internal record" with no number.

Named candidate genes from our screens are not listed, including for the public datasets. They are unpublished candidate lists and the methods do not need them.
