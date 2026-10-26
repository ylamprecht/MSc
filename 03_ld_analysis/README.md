# LD Analysis

Linkage disequilibrium across three loci, to establish how many independent
signals each contains and how that differs between populations.

## Loci

| Locus | Region | Question |
|---|---|---|
| SCN1A | chr2:166.07 Mb | How many independent signals among the five top-ranked variants in this window |
| CDKAL1 | chr6:20.63-21.06 Mb | Whether rs4710962 tags the known type 2 diabetes block |
| OAS1 | chr12:112.91-112.93 Mb | Whether rs10774671 is separable from the surrounding haplotype |

## Datasets

| Label | Description | Build |
|---|---|---|
| `cases` | Study cohort, severe disease | GRCh38 |
| `controls` | Study cohort, mild disease and highly exposed | GRCh38 |
| `GBR` | 1000 Genomes British | GRCh38 |
| `MEGA` | South African Coloured MEGA controls | GRCh37 |

Lead variant positions are given in GRCh38. MEGA is on GRCh37, where positions
are shifted by a fixed amount at each locus; this is applied through the
`offset` column of the `datasets` table in each script.

SCN1A was run on cases only. The five variants have no control carriers, so
the locus is monomorphic in controls and produces an empty heatmap.

## Scripts

| File | Purpose |
|---|---|
| `ld_analysis_SCN1A.R` | LD across the SCN1A block, cases only |
| `ld_analysis_CDKAL1.R` | LD across the CDKAL1 locus, four datasets |
| `ld_analysis_OAS1.R` | LD across the OAS1 locus, four datasets |

Each script writes one SVG heatmap and one r² matrix per dataset.

## Extracting the windows

Each script expects a VCF covering its locus, extracted a little wider than the
plotted window so the block edges remain visible:

    bcftools view -r <region> -Oz -o <locus>_<dataset>.vcf.gz <input>.vcf.gz
    bcftools index <locus>_<dataset>.vcf.gz

| Locus | GRCh38 region | GRCh37 region (MEGA) |
|---|---|---|
| SCN1A | chr2:166020000-166130000 | n/a |
| CDKAL1 | chr6:20630000-21060000 | 6:20630231-21060231 |
| OAS1 | chr12:112850000-113000000 | 12:113287805-113437805 |
