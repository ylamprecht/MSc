# ============================================================
# Binomial GLM: Allele Count and Vaccinated-Cohort Membership
#
# Comparison: V-MD vs the combined UV-HE and UV-MD control group
#
# logit(E[Y_i / n]) = b0 + b1*V-MD + b2*age + b3*sex
#                        + b4*GBR + b5*YRI + b6*NAMA
#
# Input:  top100_allele_counts.tsv
#         top100_clumped_variants.tsv
#         meta_data_with_wgrs.txt
#         meta_data_UVMD_VMD_with_wgrs.txt
# Output: vmd_adjusted_models.tsv
# ============================================================

library(data.table)

# File paths
counts_file <- "top100_allele_counts.tsv"
loci_file   <- "top100_clumped_variants.tsv"
meta_files  <- c("meta_data_with_wgrs.txt", "meta_data_UVMD_VMD_with_wgrs.txt")
out_file    <- "vmd_adjusted_models.tsv"

# Allele counts and trials
cnt <- fread(counts_file)
cnt[, SampleID := as.character(SampleID)]

loci <- fread(loci_file)
n_trials <- c(n_risk_alleles       = 2L * sum(loci$direction == "risk"),
              n_protective_alleles = 2L * sum(loci$direction == "protective"))
cat("loci retained:\n"); print(loci[, .N, by = direction])
cat(sprintf("trials per individual: %d risk, %d protective\n",
            n_trials[["n_risk_alleles"]], n_trials[["n_protective_alleles"]]))

# Metadata
anc <- rbindlist(lapply(meta_files, function(f) {
  a <- fread(f); a[, SampleID := as.character(SampleID)]
  a[, .(SampleID, age, sex, GBR, YRI, NAMA)]
}), fill = TRUE)
anc <- unique(anc, by = "SampleID")

d <- merge(cnt, anc, by = "SampleID")
d[, is_vmd := as.integer(cohort == "V-MD")]
cat("\nindividuals per cohort:\n"); print(d[, .N, by = cohort])

# Models
forms <- c("unadjusted"    = "~ is_vmd",
           "+ age"         = "~ is_vmd + age",
           "+ age and sex" = "~ is_vmd + age + sex",
           "+ ancestry"    = "~ is_vmd + age + sex + GBR + YRI + NAMA")

groups <- list("V-MD vs UV-HE and UV-MD" = d[cohort != "UV-SD"])

fit_one <- function(dt, y, f, label, grp) {
  n <- n_trials[[y]]
  dt <- copy(dt)
  dt[, `:=`(success = get(y), failure = n - get(y))]
  m <- glm(as.formula(paste("cbind(success, failure)", f)),
           family = binomial, data = dt)

  disp <- sum(residuals(m, type = "pearson")^2) / df.residual(m)
  fam  <- "binomial"
  if (disp < 0.5 || disp > 1.5) {
    m   <- update(m, family = quasibinomial)
    fam <- "quasibinomial"
  }

  co <- summary(m)$coefficients["is_vmd", ]
  ci <- suppressMessages(confint(m)["is_vmd", ])
  data.table(comparison = grp,
             outcome    = if (y == "n_risk_alleles") "Risk alleles" else "Protective alleles",
             model      = label,
             family     = fam,
             dispersion = round(disp, 2),
             beta       = round(co[1], 3),
             OR         = round(exp(co[1]), 3),
             CI_lower   = round(exp(ci[1]), 3),
             CI_upper   = round(exp(ci[2]), 3),
             p          = signif(co[4], 3))
}

# Results
res <- rbindlist(lapply(names(groups), function(g)
  rbindlist(lapply(names(n_trials), function(y)
    rbindlist(lapply(names(forms), function(f)
      fit_one(groups[[g]], y, forms[[f]], f, g)))))))

cat("\n=== binomial GLM: effect of vaccinated-cohort membership on allele count ===\n")
print(res)
fwrite(res, out_file, sep = "\t")

cat("\nOdds ratios are per-allele: the odds that any given allele is a risk\n",
    "(or protective) allele, in the vaccinated cohort relative to the\n",
    "comparison cohort.\n", sep = "")
cat(sprintf("\nwritten: %s\n", out_file))
