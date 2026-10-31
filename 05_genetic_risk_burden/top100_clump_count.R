# ============================================================
# Clump Top 100 Loci and Count Alleles
#
# Input:  top100_annotated.tsv
#         top100_calls_long.tsv
# Output: top100_allele_counts.tsv
#         top100_clumped_variants.tsv
# ============================================================

library(data.table)

# Settings
LD_WINDOW_KB <- 250
R2_MAX       <- 0.50
LD_COHORTS   <- c("UV-HE", "UV-MD", "UV-SD")

# File paths
top_file   <- "top100_annotated.tsv"
calls_file <- "top100_calls_long.tsv"
out_counts <- "top100_allele_counts.tsv"
out_loci   <- "top100_clumped_variants.tsv"

# Input
top   <- fread(top_file)
calls <- fread(calls_file)

top[, c("chrom", "pos") := .(tstrsplit(SNP, ":", fixed = TRUE)[[1]],
                             as.integer(tstrsplit(SNP, ":", fixed = TRUE)[[2]]))]
cat(sprintf("variants supplied: %d\n", nrow(top)))

# Genotype dosages
dose_alt <- c("0/0" = 0L, "0|0" = 0L,
              "0/1" = 1L, "1/0" = 1L, "0|1" = 1L, "1|0" = 1L,
              "1/1" = 2L, "1|1" = 2L,
              "./." = NA_integer_, ".|." = NA_integer_, "." = NA_integer_)

samples <- unique(calls[, .(sample, cohort)])
missing_v <- setdiff(top$SNP, unique(calls$variant))
if (length(missing_v) == nrow(top))
  stop("none of the supplied variants are in the genotype file - ",
       "run extract_top100_calls.sh first")
if (length(missing_v))
  cat(sprintf("note: %d variants have no genotype record in any sample ",
              length(missing_v)),
      "(monomorphic in this cohort) and are dropped\n", sep = "")
top <- top[!SNP %in% missing_v]

grid <- CJ(variant = top$SNP, sample = samples$sample, unique = TRUE)
grid <- merge(grid, calls[, .(variant, sample, GT)], by = c("variant", "sample"), all.x = TRUE)
grid[, alt_dose := fifelse(is.na(GT), 0L, dose_alt[GT])]

# Clumping on genotype correlation
ld_s <- samples[cohort %in% LD_COHORTS, sample]
gm   <- dcast(grid[sample %in% ld_s], sample ~ variant, value.var = "alt_dose")
gmat <- as.matrix(gm[, -1]); rownames(gmat) <- gm$sample
cat(sprintf("LD estimated in %d unvaccinated individuals\n", nrow(gmat)))

setorder(top, p, chrom, pos)
keep <- character(); dropped <- 0L
for (ch in unique(top$chrom)) {
  w <- top[chrom == ch]
  taken <- data.table(SNP = character(), pos = integer())
  for (i in seq_len(nrow(w))) {
    near <- taken[abs(pos - w$pos[i]) <= LD_WINDOW_KB * 1000]
    collapse <- FALSE
    for (j in near$SNP) {
      r2 <- suppressWarnings(cor(gmat[, w$SNP[i]], gmat[, j],
                                 use = "pairwise.complete.obs"))^2
      if (!is.na(r2) && r2 >= R2_MAX) { collapse <- TRUE; break }
    }
    if (collapse) dropped <- dropped + 1L else {
      taken <- rbind(taken, data.table(SNP = w$SNP[i], pos = w$pos[i]))
      keep  <- c(keep, w$SNP[i])
    }
  }
}
top <- top[SNP %in% keep]
cat(sprintf("clumping at %d kb, r2 >= %.2f: %d collapsed, %d loci retained\n",
            LD_WINDOW_KB, R2_MAX, dropped, nrow(top)))

# Allele counts by direction of effect
top[, direction := fifelse(OR > 1, "risk", "protective")]
cat(sprintf("  risk loci (OR > 1): %d | protective loci (OR < 1): %d\n",
            sum(top$direction == "risk"), sum(top$direction == "protective")))

grid <- grid[variant %in% top$SNP]
grid <- merge(grid, top[, .(variant = SNP, direction)], by = "variant")

cnt <- dcast(grid[, .(alleles  = sum(alt_dose, na.rm = TRUE),
                      loci     = .N,
                      missing  = sum(is.na(alt_dose))), by = .(sample, direction)],
             sample ~ direction, value.var = c("alleles", "loci", "missing"))
setnames(cnt,
         c("alleles_risk", "alleles_protective", "loci_risk", "loci_protective",
           "missing_risk", "missing_protective"),
         c("n_risk_alleles", "n_protective_alleles", "n_risk_loci",
           "n_protective_loci", "n_missing_risk", "n_missing_protective"),
         skip_absent = TRUE)

out <- merge(samples, cnt, by = "sample")
setnames(out, "sample", "SampleID")
setorder(out, cohort, SampleID)
fwrite(out, out_counts, sep = "\t")
fwrite(top[, .(SNP, chrom, pos, p, OR, car_case, car_ctrl, SYMBOL, direction)],
       out_loci, sep = "\t")

cat(sprintf("\n%d individuals scored over %d loci\n", nrow(out), nrow(top)))
print(out[, .(n = .N,
              mean_risk       = round(mean(n_risk_alleles), 1),
              mean_protective = round(mean(n_protective_alleles), 1)), by = cohort])
cat(sprintf("\nwritten: %s, %s\n", out_counts, out_loci))
