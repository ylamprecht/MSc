# ============================================================
# LD Analysis: SCN1A Block (chr2:166.07 Mb)
#
# Input:  SCN1A_cases.vcf.gz
# Output: LD_SCN1A_block.svg
#         LD_SCN1A_cases_r2.csv
# ============================================================

library(gaston)

# File paths
vcf_file <- "SCN1A_cases.vcf.gz"
out_svg  <- "LD_SCN1A_block.svg"
out_csv  <- "LD_SCN1A_cases_r2.csv"

# Window and filtering
WIN_START <- 166055000
WIN_END   <- 166095000
MAF_MIN   <- 0.05
MAX_SNPS  <- 30

# Associated variants
lead_pos <- c(166067779, 166072841, 166075785, 166076007, 166079462)

# Load
x <- read.vcf(vcf_file, convert.chr = FALSE)
x <- set.stats(x, verbose = FALSE)

# Subset to window, retaining lead variants regardless of MAF
keep <- x@snps$pos >= WIN_START & x@snps$pos <= WIN_END &
  (x@snps$maf >= MAF_MIN | x@snps$pos %in% lead_pos)
x <- x[, which(keep)]

# Thin for legibility, retaining lead variants
if (ncol(x) > MAX_SNPS) {
  idx  <- unique(round(seq(1, ncol(x), length.out = MAX_SNPS)))
  must <- which(x@snps$pos %in% lead_pos)
  idx  <- sort(unique(c(idx, must)))
  x <- x[, idx]
}

# Axis labels
x@snps$id <- ifelse(x@snps$pos %in% lead_pos,
                    paste0("* ", x@snps$pos),
                    as.character(x@snps$pos))

# LD matrix
ld  <- LD(x, lim = c(1, ncol(x)), measure = "r2")
idx <- match(lead_pos, x@snps$pos)
idx <- idx[!is.na(idx)]

m <- round(ld[idx, idx], 3)
dimnames(m) <- list(x@snps$pos[idx], x@snps$pos[idx])

cat("\nr2 among the associated variants:\n")
print(m)
cat(sprintf("\nminimum pairwise r2: %.3f\n", min(m[upper.tri(m)])))
write.csv(m, out_csv)

# Plot
svglite::svglite(out_svg, width = 11, height = 9)
LD.plot(ld,
        snp.positions = x@snps$pos,
        write.ld      = function(ld) if (ld > 0.02) round(ld, 2),
        write.snp.id  = TRUE,
        draw.chr      = TRUE,
        cex.ld        = 0.5,
        cex.snp       = 0.55,
        above.space   = 14,
        polygon.par   = list(border = "white"))
dev.off()