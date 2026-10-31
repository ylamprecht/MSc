# ============================================================
# QQ and Manhattan Plots: Autosomes
#
# Input:  chr{1-22}_firth_results.tsv
# Output: plots/qq_plot.png
#         plots/manhattan_plot.png
# ============================================================

library(data.table)
library(ggplot2)

setDTthreads(0)

# File paths
in_dir  <- "."
out_dir <- "plots"
dir.create(out_dir, showWarnings = FALSE)

# Plot settings
BASE_FAMILY <- "sans"
GRID_X      <- 3000
GRID_Y      <- 1200

theme_grafx <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      axis.text  = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.line  = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4)
    )
}

# Autosomal results
read_chr <- function(i) {
  f <- file.path(in_dir, sprintf("chr%d_firth_results.tsv", i))
  d <- fread(f, select = c("SNP", "p", "p_BH"), showProgress = FALSE)
  d[, BP  := as.numeric(tstrsplit(SNP, ":", fixed = TRUE, keep = 2)[[1]])]
  d[, CHR := i]
  d[, SNP := NULL]
  d[!is.na(p) & p > 0 & p <= 1 & !is.na(BP)]
}

res <- rbindlist(lapply(1:22, read_chr))
cat("Variants plotted:", format(nrow(res), big.mark = " "), "\n")
cat("Object size     :", format(object.size(res), units = "MB"), "\n")

# Genomic inflation
lambda <- median(qchisq(res$p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
cat("Genomic inflation factor (lambda):", round(lambda, 3), "\n")

n <- nrow(res)

# QQ plot
ord    <- order(res$p)
p_sort <- res$p[ord]
n_tail <- sum(p_sort < 0.01)
idx    <- unique(c(seq_len(max(n_tail, 1000)),
                   round(exp(seq(log(max(n_tail, 1000) + 1),
                                 log(n), length.out = 5000)))))
idx    <- idx[idx <= n]

qq <- data.table(
  obs = -log10(p_sort[idx]),
  exp = -log10(ppoints(n)[idx])
)

p_qq <- ggplot(qq, aes(exp, obs)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "#B22222", linewidth = 0.5) +
  geom_point(size = 0.7, colour = "#2E5496", alpha = 0.7) +
  annotate("text", x = min(qq$exp), y = max(qq$obs),
           hjust = 0, vjust = 1, parse = TRUE, family = BASE_FAMILY,
           label = paste0("lambda == ", round(lambda, 3))) +
  labs(x = expression(Expected ~ -log[10](italic(p))),
       y = expression(Observed ~ -log[10](italic(p)))) +
  theme_grafx()

ggsave(file.path(out_dir, "qq_plot.png"), p_qq,
       width = 5, height = 5, dpi = 300)
cat("QQ plot written\n")

# Manhattan plot
chr_len <- res[, .(len = max(BP)), by = CHR][order(CHR)]
chr_len[, offset := cumsum(as.numeric(len)) - len]
res[chr_len, BPcum := BP + i.offset, on = "CHR"]

res[, `:=`(xr = round(BPcum / max(BPcum) * GRID_X),
           yr = round(-log10(p) / max(-log10(p)) * GRID_Y))]

plot_dt <- unique(res, by = c("xr", "yr"))[, .(CHR, BPcum, p)]
setorder(plot_dt, CHR, BPcum)

axis_df <- res[, .(centre = (max(BPcum) + min(BPcum)) / 2), by = CHR][order(CHR)]

bh_line <- 0.05 / n

p_man <- ggplot(plot_dt, aes(BPcum, -log10(p), colour = factor(CHR %% 2))) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = -log10(bh_line),
             linetype = "dashed", colour = "#B22222", linewidth = 0.4) +
  scale_colour_manual(values = c("0" = "#2E5496", "1" = "#9DB4D4"),
                      guide = "none") +
  scale_x_continuous(breaks = axis_df$centre, labels = axis_df$CHR,
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +
  theme_grafx() +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(file.path(out_dir, "manhattan_plot.png"), p_man,
       width = 10, height = 4.5, dpi = 300)
cat("Manhattan plot written\n")

# Summary
cat("\nSmallest raw p       :", format(min(res$p), scientific = TRUE), "\n")
cat("Smallest BH-adj p    :", format(min(res$p_BH), scientific = TRUE), "\n")
cat("Variants at BH < 0.05:", sum(res$p_BH < 0.05), "\n")
cat("Points on Manhattan  :", format(nrow(plot_dt), big.mark = " "),
    sprintf("(%.1f%% of total)\n", 100 * nrow(plot_dt) / n))
