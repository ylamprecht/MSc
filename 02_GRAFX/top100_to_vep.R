# ============================================================
# VEP Input: Pooled Top 100
#
# Input:  top100_hits_pooled_BH.tsv
# Output: top100_for_vep.txt
# ============================================================

library(data.table)

# File paths
in_file  <- "top100_hits_pooled_BH.tsv"
out_file <- "top100_for_vep.txt"

# VCF-format records
d <- fread(in_file)

parts <- tstrsplit(d$SNP, ":", fixed = TRUE)

vep <- data.table(CHROM  = sub("^chr", "", parts[[1]]),
                  POS    = as.integer(parts[[2]]),
                  ID     = d$SNP,
                  REF    = parts[[3]],
                  ALT    = parts[[4]],
                  QUAL   = ".",
                  FILTER = ".",
                  INFO   = ".")

# Chromosome and position order
vep[, ord := match(CHROM, c(as.character(1:22), "X", "Y", "MT"))]
setorder(vep, ord, POS)
vep[, ord := NULL]

fwrite(vep, out_file, sep = "\t", col.names = FALSE, quote = FALSE)

cat(sprintf("%d variants -> %s\n", nrow(vep), out_file))
cat(sprintf("SNVs: %d | indels: %d\n",
            sum(nchar(vep$REF) == 1 & nchar(vep$ALT) == 1),
            sum(nchar(vep$REF) > 1 | nchar(vep$ALT) > 1)))
cat(sprintf("chromosomes: %s\n", paste(unique(vep$CHROM), collapse = ", ")))
