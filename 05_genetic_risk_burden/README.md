# Genetic risk burden

Comparison of genetic risk burden across the study cohorts, using published
COVID-19 risk variants and the top-ranked variants from this study's own
association analysis.

## Scripts

| File | Description |
| --- | --- |
| `extract_burden_calls.sh` | Genotype calls at the published risk variant positions |
| `wgrs.R` | Weighted genetic risk score and allele counts per individual |
| `wgrs_regression.R` | Firth regression of case status on the severity wGRS |
| `wgrs_sum_comparison.R` | Descriptive statistics, Wilcoxon tests and box plots |
| `extract_top100_calls.sh` | Genotype calls at the top 100 association positions |
| `top100_clump_count.R` | LD clumping and risk/protective allele counts |
| `top100_comparison.R` | Descriptive statistics, Wilcoxon tests and box plots |
| `top100_vmd_adjusted.R` | Binomial GLM adjusted for age, sex and ancestry |

## Input files

Not tracked in this repository. Copy them into this folder before running.

| File | Source |
| --- | --- |
| `published_risk_variants.csv` | exported in `04_COVID_19_risk_variant_db` |
| `top100_annotated.tsv` | produced in `02_GRAFX` |
| `meta_data.txt` | prepared manually |
| `meta_data_with_wgrs.txt` | `meta_data.txt` with the wGRS columns appended |
| `meta_data_UVMD_VMD_with_wgrs.txt` | the UV-MD and V-MD subset |
| `reinfected_samples.txt` | one sample ID per line, prepared manually |

Per-sample VCFs are read from `../CVD_vcfs/<cohort>/`. All scripts are run from
inside this folder.

## Software

| Tool | Version |
| --- | --- |
| bcftools | 1.19 |
| R | 4.6.0 |

## 1. Published risk variants

```bash
bash extract_burden_calls.sh
Rscript wgrs.R
Rscript wgrs_regression.R
Rscript wgrs_sum_comparison.R
```

| Output | Contents |
| --- | --- |
| `burden_calls_long.tsv` | Genotype calls per variant and sample |
| `wgrs_scores.tsv` | wGRS and allele counts per individual |
| `wgrs_variant_audit.tsv` | Variants retained, with orientation and weight |
| `firth_wgrs_severity_all_models_UVMD_VMD.tsv` | All regression coefficients |
| `firth_wgrs_severity_model_comparison_UVMD_VMD.tsv` | The wGRS term across models |
| `wgrs_descriptive_summary.tsv` | Means, medians and IQRs by cohort |
| `wgrs_wilcoxon.tsv` | UV-SD versus V-MD rank-sum tests |
| `wgrs_boxplots.svg` | Allele counts and scores by cohort |

## 2. Study-identified top 100

```bash
bash extract_top100_calls.sh
Rscript top100_clump_count.R
Rscript top100_comparison.R
Rscript top100_vmd_adjusted.R
```

| Output | Contents |
| --- | --- |
| `top100_calls_long.tsv` | Genotype calls per variant and sample |
| `top100_allele_counts.tsv` | Risk and protective allele counts per individual |
| `top100_clumped_variants.tsv` | Loci retained after clumping, with direction |
| `top100_count_summary.tsv` | Means, medians and IQRs by cohort |
| `top100_count_wilcoxon.tsv` | V-MD against each unvaccinated cohort, BH corrected |
| `top100_count_boxplots.svg` | Allele counts by cohort |
| `vmd_adjusted_models.tsv` | Binomial GLM, V-MD versus the combined controls |

Together these compare the burden of published and study-identified risk alleles
between the vaccinated and unvaccinated cohorts.
