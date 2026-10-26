# Ancestry Inference: Command Log

Steps performed interactively on the CHPC. Run separately for each cohort;
paths below are from the UV-HE run.

Software: bcftools 1.21, PLINK 1.90 and 2.00a, KING, ADMIXTURE, pong.

Interactive session:

    qsub -I -l select=1:ncpus=24:mpiprocs=24 -q smp -P HEAL1360 -l walltime=10:00:00

## 1. Merge VCFs and convert to PLINK

Performed using `merge_batch.pbs` in `00_variant_calling_qc/`. Per-sample VCFs
were compressed, renamed, indexed, merged and normalised to split multi-allelic
sites into biallelic records.

## 2. Genotype QC

    plink2 --bfile 1000G_Nama_Malay_reference \
           --geno 0.05 --mind 0.1 --maf 0.0010373 --hwe 0.00001 \
           --make-bed --out Reference_QC

    plink2 --bfile Target \
           --geno 0.05 --mind 0.1 --maf 0.0147 --hwe 0.00001 \
           --allow-extra-chr \
           --make-bed --out Target_QC

    Note: maf is set to 1/2n for each dataset, where n is the number of
    individuals (wc -l on the .fam file).

## 3. Merge reference and target

Only intersecting SNPs can be merged, so these are identified and extracted
from both datasets first.

    plink2 --bfile Target_QC --set-all-var-ids '@:#' --allow-extra-chr \
           --make-bed --out Target_QC_ChrBp
    plink2 --bfile Reference_QC --set-all-var-ids '@:#' \
           --make-bed --out Reference_QC_ChrBp

    awk '{print $2}' Target_QC_ChrBp.bim    | sort > Target_SNP_list.txt
    awk '{print $2}' Reference_QC_ChrBp.bim | sort > Reference_SNP_list.txt
    comm -12 Target_SNP_list.txt Reference_SNP_list.txt > Intersecting_SNPs.txt

    plink2 --bfile Target_QC_ChrBp --extract Intersecting_SNPs.txt \
           --allow-extra-chr --make-bed --out Target_IntersectingSNPs
    plink2 --bfile Reference_QC_ChrBp --extract Intersecting_SNPs.txt \
           --make-bed --out Reference_IntersectingSNPs

    plink --bfile Target_IntersectingSNPs \
          --bmerge Reference_IntersectingSNPs.bed \
                   Reference_IntersectingSNPs.bim \
                   Reference_IntersectingSNPs.fam \
          --make-bed --out Merged_files

    Note: if the merge fails with "variants with 3+ alleles present",
    exclude the offending variants from both datasets and merge again:

        plink2 --bfile Target_IntersectingSNPs \
               --exclude Merged_files-merge.missnp \
               --make-bed --out Target_IntersectingSNPs_excl_multi
        plink2 --bfile Reference_IntersectingSNPs \
               --exclude Merged_files-merge.missnp \
               --make-bed --out Reference_IntersectingSNPs_excl_multi

QC the merged file to remove residual genotyping errors:

    plink2 --bfile Merged_files \
           --geno 0.05 --mind 0.1 --hwe 0.00001 \
           --make-bed --out Merged_files_QC

Imputation was not required, as the target data are whole-genome sequences.

## 4. Relatedness

    king -b Merged_files_QC.bed --unrelated --degree 2

    plink2 --bfile Merged_files_QC \
           --remove kingunrelated_toberemoved.txt \
           --make-bed --out Merged_files_QC_RemoveRelateds

## 5. LD pruning

r² > 0.1 within a 50 SNP window, advanced 10 SNPs at a time.

    plink2 --bfile Merged_files_QC_RemoveRelateds \
           --indep-pairwise 50 10 0.1 --out LD-pruned

    plink2 --bfile Merged_files_QC_RemoveRelateds \
           --extract LD-pruned.prune.in \
           --make-bed --out Merged_files_QC_RemoveRelateds_pruned

    plink2 --bfile Merged_files_QC_RemoveRelateds_pruned \
           --keep ind2keep.txt \
           --make-bed --out UVHE_Ancestry_final

## 6. PCA

    plink2 --bfile Merged_files_QC_RemoveRelateds_pruned \
           --pca 10 --out Merged_files_QC_RemoveRelateds_pruned_PCA

Plotted with `ancestry_pca.R`.

## 7. ADMIXTURE

Performed in `admixture.pbs`, for K = 2 to 10. The optimal K was selected as
the value with the lowest cross-validation error.

## 8. Visualisation

Run locally:

    pong -m pong_filemap.txt -i ind2pop.txt -n pop_order.txt

    Note: ind2pop.txt must be in the same order as the .fam file.
