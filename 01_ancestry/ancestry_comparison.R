# ============================================================
# Global Ancestry Comparison Between Cohorts
#
# Input:  ancestry_by_cohort.tsv
# Output: ancestry_summary.tsv
#         ancestry_kruskal.tsv
#         ancestry_boxplots.svg
# ============================================================

library(data.table)

# Ancestry components
comps <- c(NAMA  = "KhoeSan",
           YRI   = "West African",
           GBR   = "European",
           ASIAN = "Asian")

# Load ancestry proportions
anc <- fread("ancestry_by_cohort.tsv")
anc[, cohort := factor(cohort, levels = c("UV-HE", "UV-MD", "UV-SD", "V-MD"))]

cat("individuals per cohort:\n"); print(table(anc$cohort))

# Median (IQR) per cohort, as percentages
fmt <- function(x) sprintf("%.1f (%.1f-%.1f)",
                           100 * median(x),
                           100 * quantile(x, 0.25),
                           100 * quantile(x, 0.75))

summ <- anc[, lapply(.SD, fmt), by = cohort, .SDcols = names(comps)][order(cohort)]
setnames(summ, names(comps), paste0(comps, " %"))
cat("\nmedian (IQR) ancestry proportion, %:\n"); print(summ)

# Kruskal-Wallis per component
kw <- rbindlist(lapply(names(comps), function(k) {
  kt <- kruskal.test(anc[[k]] ~ anc$cohort)
  data.table(component = comps[[k]],
             H         = round(kt$statistic, 2),
             df        = kt$parameter,
             p         = kt$p.value)
}))
kw[, p_BH := p.adjust(p, method = "BH")]
cat("\nKruskal-Wallis:\n")
print(kw[, .(component, H, df, p = signif(p, 3), p_BH = signif(p_BH, 3))])

# Pairwise Wilcoxon where Kruskal-Wallis is significant
for (k in names(comps)[kw$p_BH < 0.05]) {
  pw <- pairwise.wilcox.test(anc[[k]], anc$cohort,
                             p.adjust.method = "BH",
                             exact           = FALSE)
  cat(sprintf("\npairwise Wilcoxon, %s (BH-adjusted p):\n", comps[[k]]))
  print(signif(pw$p.value, 3))
}

# Save tables
fwrite(summ, "ancestry_summary.tsv", sep = "\t")
fwrite(kw,   "ancestry_kruskal.tsv", sep = "\t")

# Boxplots
svglite::svglite("ancestry_boxplots.svg", width = 10, height = 9)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
set.seed(1)
for (k in names(comps)) {
  y <- 100 * anc[[k]]
  boxplot(y ~ anc$cohort,
          main    = comps[[k]],
          xlab    = "Cohort",
          ylab    = "Ancestry (%)",
          ylim    = c(0, 100),
          outline = FALSE,
          border  = NA,
          col     = NA)
  abline(h = seq(0, 100, 20), col = "grey90", lwd = 1)
  boxplot(y ~ anc$cohort, add = TRUE, col = "white", border = "black",
          outline = FALSE, xaxt = "n", yaxt = "n")
  stripchart(y ~ anc$cohort, vertical = TRUE, method = "jitter",
             pch = 16, cex = 0.5, add = TRUE, col = "black")
  box()
}
invisible(dev.off())