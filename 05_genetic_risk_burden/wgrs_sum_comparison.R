# ============================================================
# Summary and Comparison of the Weighted Genetic Risk Scores
#
# Input:  wgrs_scores.tsv
#         reinfected_samples.txt
# Output: wgrs_descriptive_summary.tsv
#         wgrs_wilcoxon.tsv
#         wgrs_boxplots.svg
# ============================================================

library(data.table)

# File paths
in_file    <- "wgrs_scores.tsv"
reinf_file <- "reinfected_samples.txt"
out_summ   <- "wgrs_descriptive_summary.tsv"
out_wilcox <- "wgrs_wilcoxon.tsv"
out_fig    <- "wgrs_boxplots.svg"

# Scores
dat <- fread(in_file)
dat[, SampleID := as.character(SampleID)]
dat[, cohort := factor(cohort, levels = c("UV-HE", "UV-MD", "UV-SD", "V-MD"))]

# Vaccinated participants with documented multiple infections
reinfected <- trimws(readLines(reinf_file))
reinfected <- reinfected[nzchar(reinfected)]
missing_reinfected <- setdiff(reinfected, dat$SampleID)
if (length(missing_reinfected))
  cat("WARNING: reinfected samples not found: ",
      paste(missing_reinfected, collapse = ", "), "\n", sep = "")

cat("samples per cohort:\n"); print(table(dat$cohort))

# Measures
measures <- c(
  n_risk_alleles_severity             = "Severity risk alleles",
  n_protective_alleles_severity       = "Severity protective alleles",
  n_risk_alleles_susceptibility       = "Susceptibility risk alleles",
  n_protective_alleles_susceptibility = "Susceptibility protective alleles",
  wGRS_severity                       = "Severity wGRS",
  wGRS_susceptibility                 = "Susceptibility wGRS")
measures <- measures[names(measures) %in% names(dat)]

# Descriptive statistics
describe <- function(dt, grp_col, grp_label) {
  rbindlist(lapply(names(measures), function(k)
    dt[, .(measure = measures[[k]],
           group   = grp_label(get(grp_col)),
           n       = sum(!is.na(get(k))),
           mean    = round(mean(get(k), na.rm = TRUE), 2),
           median  = round(as.double(median(get(k), na.rm = TRUE)), 2),
           IQR     = sprintf("%.2f-%.2f",
                             quantile(get(k), 0.25, na.rm = TRUE),
                             quantile(get(k), 0.75, na.rm = TRUE))),
       by = c(grp_col)][, !..grp_col]))
}

summ <- describe(dat, "cohort", as.character)

dat[, reinf_grp := fifelse(cohort == "V-MD" & SampleID %in% reinfected,
                           "V-MD (multiple infections)", NA_character_)]
summ_reinf <- describe(dat[!is.na(reinf_grp)], "reinf_grp", as.character)

summ <- rbind(summ, summ_reinf)
summ[, measure := factor(measure, levels = unname(measures))]
summ[, group := factor(group, levels = c("UV-HE", "UV-MD", "UV-SD", "V-MD",
                                         "V-MD (multiple infections)"))]
setorder(summ, measure, group)
setcolorder(summ, c("measure", "group"))

cat("\n=== descriptive statistics ===\n"); print(summ)
fwrite(summ, out_summ, sep = "\t")

# Wilcoxon rank-sum tests
score_cols <- intersect(c("wGRS_severity", "wGRS_susceptibility"), names(dat))

comparisons <- list(
  "UV-SD vs V-MD" = list(a = dat[cohort == "UV-SD"], b = dat[cohort == "V-MD"]),
  "UV-SD vs V-MD (multiple infections)" =
    list(a = dat[cohort == "UV-SD"],
         b = dat[cohort == "V-MD" & SampleID %in% reinfected]))

wilcox_results <- rbindlist(lapply(names(comparisons), function(cmp) {
  g <- comparisons[[cmp]]
  rbindlist(lapply(score_cols, function(k) {
    x <- g$a[[k]][!is.na(g$a[[k]])]
    y <- g$b[[k]][!is.na(g$b[[k]])]
    wt <- wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
    data.table(comparison = cmp, measure = measures[[k]],
               n_UV_SD = length(x), n_comparison = length(y),
               median_UV_SD = round(as.double(median(x)), 2),
               median_comparison = round(as.double(median(y)), 2),
               W = unname(wt$statistic), p = round(wt$p.value, 4))
  }))
}))

cat("\n=== Wilcoxon rank-sum tests ===\n"); print(wilcox_results)
fwrite(wilcox_results, out_wilcox, sep = "\t")

# Box plots
panel_cols <- c(n_risk_alleles_severity       = "Severity risk alleles",
                wGRS_severity                 = "Severity wGRS",
                n_risk_alleles_susceptibility = "Susceptibility risk alleles",
                wGRS_susceptibility           = "Susceptibility wGRS")
panel_cols <- panel_cols[names(panel_cols) %in% names(dat)]

panel <- function(y, grp, main, ylab) {
  yl <- range(pretty(range(y, na.rm = TRUE)))
  boxplot(y ~ grp, main = main, xlab = "Cohort", ylab = ylab,
          ylim = yl, outline = FALSE, border = NA, col = NA)
  abline(h = pretty(yl), col = "grey90", lwd = 1)
  boxplot(y ~ grp, add = TRUE, col = "white", border = "black",
          outline = FALSE, xaxt = "n", yaxt = "n")
  stripchart(y ~ grp, vertical = TRUE, method = "jitter",
             pch = 16, cex = 0.5, add = TRUE, col = "black")
  box()
}

svglite::svglite(out_fig, width = 10, height = 9)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
set.seed(1)
for (k in names(panel_cols))
  panel(dat[[k]], dat$cohort, panel_cols[[k]],
        if (grepl("^wGRS", k)) "wGRS" else "Alleles (count)")
invisible(dev.off())

cat("\nWritten to:\n")
cat(sprintf("  %s\n  %s\n  %s\n", out_summ, out_wilcox, out_fig))
