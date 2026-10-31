# ============================================================
# Firth Regression: Severity wGRS
#
# Cases:    V-MD  (phenotype = 1)
# Controls: UV-MD (phenotype = 0)
#
# Model 1: phenotype ~ wGRS_severity_z
# Model 2: phenotype ~ wGRS_severity_z + age + sex
# Model 3: phenotype ~ wGRS_severity_z + age + sex + GBR + YRI + NAMA
#
# Input:  meta_data_UVMD_VMD_with_wgrs.txt
# Output: firth_wgrs_severity_all_models_UVMD_VMD.tsv
#         firth_wgrs_severity_model_comparison_UVMD_VMD.tsv
# ============================================================

library(data.table)
library(logistf)

# File paths
meta_file <- "meta_data_UVMD_VMD_with_wgrs.txt"
out_all   <- "firth_wgrs_severity_all_models_UVMD_VMD.tsv"
out_wgrs  <- "firth_wgrs_severity_model_comparison_UVMD_VMD.tsv"

# Metadata
meta <- fread(meta_file)
meta$sex <- factor(meta$sex)

cat(sprintf("Samples: %d cases, %d controls, %d total\n",
            sum(meta$phenotype == 1),
            sum(meta$phenotype == 0),
            nrow(meta)))

required <- c("phenotype", "wGRS_severity_z", "age", "sex", "GBR", "YRI", "NAMA")

missing_cols <- setdiff(required, names(meta))
if (length(missing_cols))
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

if (any(!complete.cases(meta[, ..required])))
  stop("Missing values present in variables used for regression.")

# Models
model1 <- logistf(phenotype ~ wGRS_severity_z, data = meta)
model2 <- logistf(phenotype ~ wGRS_severity_z + age + sex, data = meta)
model3 <- logistf(phenotype ~ wGRS_severity_z + age + sex + GBR + YRI + NAMA, data = meta)

# Results
extract_results <- function(fit, model_name) {
  data.table(model     = model_name,
             predictor = names(coef(fit)),
             beta      = unname(coef(fit)),
             OR        = unname(exp(coef(fit))),
             CI_lower  = unname(exp(fit$ci.lower)),
             CI_upper  = unname(exp(fit$ci.upper)),
             p         = unname(fit$prob))
}

result <- rbindlist(list(
  extract_results(model1, "Model 1: wGRS"),
  extract_results(model2, "Model 2: + age + sex"),
  extract_results(model3, "Model 3: + age + sex + ancestry")
))

result <- result[predictor != "(Intercept)"]
fwrite(result, out_all, sep = "\t")

cat("\nFull Firth logistic regression results:\n")
print(result)

# Summary
wgrs_result <- result[predictor == "wGRS_severity_z",
                      .(model, beta, OR, CI_lower, CI_upper, p)]

cat("\nSeverity wGRS across models:\n")
print(wgrs_result)

cat("\nFormatted wGRS results:\n")
for (i in seq_len(nrow(wgrs_result))) {
  cat(sprintf("%s: OR = %.3f (95%% CI %.3f-%.3f), p = %.4g\n",
              wgrs_result$model[i],
              wgrs_result$OR[i],
              wgrs_result$CI_lower[i],
              wgrs_result$CI_upper[i],
              wgrs_result$p[i]))
}

fwrite(wgrs_result, out_wgrs, sep = "\t")

cat("\nWritten to:\n")
cat(sprintf("  %s\n", out_all))
cat(sprintf("  %s\n", out_wgrs))
