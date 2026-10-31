# ============================================================
# Extract GWAMA Meta-Analysis Results
#
# Input:  chrX_meta.out
# Output: chrX_combined_gwama_results.tsv
# ============================================================

library(data.table)

# File paths
in_file  <- "chrX_meta.out"
out_file <- "chrX_combined_gwama_results.tsv"

# Columns of interest
meta <- fread(in_file)

res <- meta[, .(SNP       = rs_number,
                p         = `gender_differentiated_p-value`,
                p_het     = `gender_heterogeneity_p-value`,
                female_OR = female_OR,
                male_OR   = male_OR)]

setorder(res, p, na.last = TRUE)
fwrite(res, out_file, sep = "\t")

cat(sprintf("%d variants written to %s\n", nrow(res), out_file))
cat(sprintf("missing p: %d | smallest p: %.3g\n",
            sum(is.na(res$p)), min(res$p, na.rm = TRUE)))
