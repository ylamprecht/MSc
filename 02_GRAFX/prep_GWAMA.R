# ============================================================
# GWAMA Input: Chromosome X
#
# Input:  chrX_female_firth_results.tsv
#         chrX_male_firth_results.tsv
# Output: chrX_female.gwama.txt
#         chrX_male.gwama.txt
#         gwama.in
# ============================================================

library(data.table)

# File paths
female_file <- "chrX_female_firth_results.tsv"
male_file   <- "chrX_male_firth_results.tsv"
female_out  <- "chrX_female.gwama.txt"
male_out    <- "chrX_male.gwama.txt"
list_out    <- "gwama.in"

# Firth results, dropping unbounded confidence intervals
f <- fread(female_file)
m <- fread(male_file)

usable <- function(d) d[is.finite(OR_95L) & is.finite(OR_95U) & OR_95L > 0 & OR > 0]
f <- usable(f)
m <- usable(m)

# Intersecting variants
shared <- intersect(f$SNP, m$SNP)
cat(sprintf("%d shared | %d female-only | %d male-only\n",
            length(shared), sum(!f$SNP %in% m$SNP), sum(!m$SNP %in% f$SNP)))

# GWAMA input
gwama_fmt <- function(d) d[SNP %in% shared,
  .(MARKERNAME = SNP, EA, NEA, OR, OR_95L, OR_95U, N, EAF, STRAND = "+")]

fwrite(gwama_fmt(f), female_out, sep = "\t")
fwrite(gwama_fmt(m), male_out,   sep = "\t")
writeLines(c(paste0(female_out, "\tF"), paste0(male_out, "\tM")), list_out)
