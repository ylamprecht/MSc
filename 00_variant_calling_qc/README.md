# Variant Calling and Quality Control

Variant calling, batch merging, and assessment of sequencing batch effects.

Samples were sequenced in two batches:

| Batch | Platform | Coverage | n |
|---|---|---|---|
| 1 | DNBSeq-G400 | 30X | 85 |
| 2 | DNBSeq-T7 | 24X | 42 |

## Scripts

| File | Purpose |
|---|---|
| `prep_genome.pbs` | Reference genome preparation and indexing |
| `variant_calling.pbs` | Per-sample variant calling |
| `merge_batch.pbs` | Merge, index and normalise VCFs (bcftools v1.19) |
| `batch_effect.R` | PCA visualisation by sequencing batch |
| `parent_multiqc_report.html` | Sequencing QC summary, batch 1 |
| `beta_multiqc_report.html` | Sequencing QC summary, batch 2 |

## Batch effect assessment

Steps below were run interactively and are recorded here rather than as scripts.

Multi-allelic sites were split into biallelic records during merging
(`merge_batch.pbs`).

Convert to PLINK format, autosomes only:

    plink2 --vcf merged_norm.vcf.gz \
           --autosome \
           --set-all-var-ids '@:#$r$a' \
           --snps-only just-acgt \
           --allow-extra-chr \
           --make-bed \
           --out merged_batches_autosomes

Filter on missingness and allele frequency:

    plink2 --bfile merged_batches_autosomes \
           --maf 0.003937 \
           --geno 0.05 \
           --allow-extra-chr \
           --make-bed \
           --out merged_batches_autosomes_filtered

    30 498 256 variants removed (missingness)
    42 486 variants removed (allele frequency)

LD pruning:

    plink2 --bfile merged_batches_autosomes_filtered \
           --indep-pairwise 50 10 0.1 \
           --out LD-pruned

    plink2 --bfile merged_batches_autosomes_filtered \
           --extract LD-pruned.prune.in \
           --allow-extra-chr \
           --make-bed \
           --out merged_batches_autosomes_filtered_pruned

    72 811 variants remaining

Principal component analysis:

    plink2 --bfile merged_batches_autosomes_filtered_pruned \
           --pca 10 \
           --out merged_batches_autosomes_filtered_pruned_PCA

Principal components were plotted by sequencing batch (`batch_effect.R`).
No clustering by batch was observed.
