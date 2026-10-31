# ============================================================
# Pooled Benjamini-Hochberg Correction
#
# Input:  chr{1-22}_firth_results.tsv
#         chrX_combined_gwama_results.tsv
# Output: top100_hits_pooled_BH.tsv
#         p_BH overwritten in place in every input file
# ============================================================

library(data.table)

# Settings
Q_THRESHOLD <- 0.05
TOP_N       <- 100

# File paths
autosome_files <- sprintf("chr%d_firth_results.tsv", 1:22)
chrX_file      <- "chrX_combined_gwama_results.tsv"
all_files      <- c(autosome_files, chrX_file)

missing <- all_files[!file.exists(all_files)]
if (length(missing))
  stop("missing input files: ", paste(missing, collapse = ", "))

# Pooled p-values
cat("Pass 1: reading p-values\n")
pooled <- rbindlist(lapply(all_files, function(f) {
  d <- fread(f, select = c("SNP", "p"))
  cat(sprintf("  %-35s %9d variants (%d with NA p)\n", f, nrow(d), sum(is.na(d$p))))
  data.table(SNP = d$SNP, p = d$p, file = f)
}))

if (anyDuplicated(pooled$SNP))
  warning(sum(duplicated(pooled$SNP)), " duplicated SNP IDs across files - ",
          "p_BH will be written to every matching row")

# BH correction
tested <- !is.na(pooled$p)
pooled[, p_BH := NA_real_]
pooled[tested, p_BH := p.adjust(p, method = "BH")]

m <- sum(tested)
cat(sprintf("\nBH correction over m = %d tests (%d excluded as NA)\n",
            m, sum(!tested)))
cat(sprintf("smallest p = %.3g | smallest p_BH = %.3g\n",
            min(pooled$p, na.rm = TRUE), min(pooled$p_BH, na.rm = TRUE)))
cat(sprintf("BH threshold for the top-ranked variant (q/m) = %.3g\n",
            Q_THRESHOLD / m))
cat(sprintf("variants with p_BH < %.2f: %d\n",
            Q_THRESHOLD, sum(pooled$p_BH < Q_THRESHOLD, na.rm = TRUE)))

# Write p_BH back
setorder(pooled, p, na.last = TRUE)
top_snps <- pooled$SNP[seq_len(min(TOP_N, nrow(pooled)))]

bh <- pooled[, .(SNP, file, p_BH)]
setkey(bh, SNP)

cat("\nPass 2: writing p_BH back\n")
top_rows <- vector("list", length(all_files))

for (i in seq_along(all_files)) {
  f <- all_files[i]
  d <- fread(f)

  d[bh[file == f], on = "SNP", p_BH := i.p_BH]

  tmp <- paste0(f, ".tmp")
  fwrite(d, tmp, sep = "\t")
  if (!file.rename(tmp, f)) stop("could not replace ", f)
  cat(sprintf("  %-35s p_BH written\n", f))

  hit <- d[SNP %in% top_snps]
  if (nrow(hit)) {
    hit[, source := f]
    top_rows[[i]] <- hit
  }
}

# Top hits
top <- rbindlist(top_rows, fill = TRUE)
setorder(top, p, na.last = TRUE)
top <- top[seq_len(min(TOP_N, nrow(top)))]

front <- intersect(c("SNP", "source", "p", "p_BH"), names(top))
setcolorder(top, front)
fwrite(top, "top100_hits_pooled_BH.tsv", sep = "\t")

cat(sprintf("\n%d rows written to top100_hits_pooled_BH.tsv\n", nrow(top)))
print(top[, .N, by = source])
