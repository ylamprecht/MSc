# Ancestry Inference

Global ancestry inference for the study cohorts, using a reference panel of
1000 Genomes, Nama and Malay populations.

Ancestry was inferred separately for each cohort (UV-HE, UV-MD, UV-SD, V-MD)
to maintain a balanced ratio of target to reference individuals.

Steps follow the TB Host Genetics
[PopGen workshop](https://github.com/TBHostGenetics/PopGen-workshop).

VCFs were merged using `merge_batch.pbs` in `00_variant_calling_qc/`.

## Scripts

| File | Purpose |
|---|---|
| `admixture.pbs` | ADMIXTURE runs across K values |
| `ancestry_pca.R` | PCA of target and reference populations |
| `ancestry_comparison.R` | Compare ancestry proportions between cohorts |
| `admixture_sex_compare.R` | Compare ancestry proportions between sexes |
| `ancestry_inference_steps.md` | Full command log for the interactive steps |

## Input files

| File | Format |
|---|---|
| `cohorts.txt` | Two columns, no header: IID, population or cohort |
| `ind2keep.txt` | Individual IDs to retain after pruning |
| `ind2pop.txt` | Population label per individual, same order as the `.fam` file |
| `pop_order.txt` | Population display order for pong |
