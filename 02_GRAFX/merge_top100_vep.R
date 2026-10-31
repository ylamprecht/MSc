# ============================================================
# Merge Top 100 with VEP Annotation
#
# Input:  top100_hits_pooled_BH.tsv
#         vep_annotated_top100.txt
# Output: top100_annotated.tsv
# ============================================================

library(data.table)

# File paths
res_file <- "top100_hits_pooled_BH.tsv"
vep_file <- "vep_annotated_top100.txt"
out_file <- "top100_annotated.tsv"

# VEP columns retained
keep <- c("SNP", "Consequence", "SYMBOL", "Gene", "BIOTYPE", "Existing_variation",
          "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF",
          "gnomADg_ASJ_AF", "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF",
          "gnomADg_NFE_AF", "gnomADg_REMAINING_AF", "gnomADg_SAS_AF")

res <- fread(res_file)
vep <- fread(vep_file)

setnames(vep, "#Uploaded_variation", "SNP")

missing_cols <- setdiff(keep, names(vep))
if (length(missing_cols))
  stop("VEP file is missing: ", paste(missing_cols, collapse = ", "))

vep <- vep[, ..keep]
stopifnot(!anyDuplicated(vep$SNP))

# Merge
if ("source" %in% names(res)) res[, source := NULL]

out <- merge(res, vep, by = "SNP", all.x = TRUE, sort = FALSE)

freq_cols <- grep("^gnomADg_", names(out), value = TRUE)
out[, (freq_cols) := lapply(.SD, function(x)
      suppressWarnings(as.numeric(ifelse(x == "-", NA, x)))), .SDcols = freq_cols]

setorder(out, p, na.last = TRUE)
fwrite(out, out_file, sep = "\t")

# Diagnostics
n_noannot <- sum(is.na(out$Consequence))
cat(sprintf("%d variants written to %s\n", nrow(out), out_file))
cat(sprintf("columns: %s\n", paste(names(out), collapse = ", ")))
if (n_noannot > 0)
  cat(sprintf("WARNING: %d variants had no VEP annotation\n", n_noannot))
cat(sprintf("no gnomAD AFR frequency: %d\n", sum(is.na(out$gnomADg_AFR_AF))))
