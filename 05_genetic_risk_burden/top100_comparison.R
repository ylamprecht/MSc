# ============================================================
# Compare Top 100 Allele Counts Across Cohorts
#
# Input:  top100_allele_counts.tsv
# Output: top100_count_summary.tsv
#         top100_count_wilcoxon.tsv
#         top100_count_boxplots.svg
# ============================================================

library(data.table)

# File paths
in_file    <- "top100_allele_counts.tsv"
out_summ   <- "top100_count_summary.tsv"
out_wilcox <- "top100_count_wilcoxon.tsv"
out_fig    <- "top100_count_boxplots.svg"

# Allele counts
d <- fread(in_file)
d[, cohort := factor(cohort, levels = c("UV-HE", "UV-MD", "UV-SD", "V-MD"))]

measures <- c(n_risk_alleles       = "Risk alleles",
              n_protective_alleles = "Protective alleles")

cat("samples per cohort:\n"); print(table(d$cohort))

# Descriptive statistics
summ <- rbindlist(lapply(names(measures), function(k)
  d[, .(measure = measures[[k]], n = .N,
        mean    = round(mean(get(k)), 2),
        median  = round(as.double(median(get(k))), 2),
        IQR     = sprintf("%.2f-%.2f", quantile(get(k), .25), quantile(get(k), .75))),
    by = cohort]))
summ[, measure := factor(measure, levels = unname(measures))]
setorder(summ, measure, cohort)
setcolorder(summ, c("measure", "cohort"))
cat("\n=== allele counts per cohort ===\n"); print(summ)
fwrite(summ, out_summ, sep = "\t")

# Wilcoxon rank-sum tests
comparisons <- c("UV-SD", "UV-HE", "UV-MD")

res <- rbindlist(lapply(comparisons, function(g)
  rbindlist(lapply(names(measures), function(k) {
    x <- d[cohort == "V-MD"][[k]]
    y <- d[cohort == g][[k]]
    wt <- wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
    data.table(comparison = paste("V-MD vs", g), measure = measures[[k]],
               n_VMD = length(x), n_other = length(y),
               median_VMD   = round(as.double(median(x)), 2),
               median_other = round(as.double(median(y)), 2),
               W = unname(wt$statistic), p = round(wt$p.value, 4))
  }))))
res[, p_BH := round(p.adjust(p, "BH"), 4), by = measure]

cat("\n=== Wilcoxon rank-sum tests ===\n"); print(res)
fwrite(res, out_wilcox, sep = "\t")

# Box plots
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

svglite::svglite(out_fig, width = 10, height = 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
set.seed(1)
for (k in names(measures))
  panel(d[[k]], d$cohort, measures[[k]], "Alleles (count)")
invisible(dev.off())

cat(sprintf("\nwritten: %s, %s,\n  %s\n", out_summ, out_wilcox, out_fig))
