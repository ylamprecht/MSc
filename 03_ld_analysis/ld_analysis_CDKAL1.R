# ============================================================
# LD Analysis: CDKAL1 Locus (chr6:20.63-21.06 Mb)
#
# Input:  cases_CDKAL1.vcf.gz
#         controls_CDKAL1.vcf.gz
#         GBR_CDKAL1.vcf.gz
#         MEGA_CDKAL1.vcf.gz
# Output: LD_CDKAL1_<dataset>.svg
#         LD_CDKAL1_<dataset>_r2.csv
# ============================================================

library(gaston)

# Datasets
datasets <- data.frame(
  label            = c("cases", "controls", "GBR", "MEGA"),
  vcf              = c("cases_CDKAL1.vcf.gz",
                       "controls_CDKAL1.vcf.gz",
                       "GBR_CDKAL1.vcf.gz",
                       "MEGA_CDKAL1.vcf.gz"),
  offset           = c(0, 0, 0, 231),
  stringsAsFactors = FALSE
)

# Window and filtering
WIN_START <- 20630000
WIN_END   <- 21060000
MAF_MIN   <- 0.05
MAX_SNPS  <- 40

# Lead variants: rs4710962, then the CDKAL1 T2D leads
lead_pos <- c(21031204,
              20686342, 20687890, 20686765, 20680447,
              20682933, 20679478, 20661019, 20660803)

lead_name <- c("rs4710962",
               "rs7766070", "rs10440833", "rs9368222", "rs9348441",
               "rs9350271", "rs7756992",  "rs7754840", "rs10946398")

for (i in seq_len(nrow(datasets))) {
  
  label  <- datasets$label[i]
  offset <- datasets$offset[i]
  pos    <- lead_pos + offset
  
  cat(sprintf("\n=== %s ===\n", label))
  
  # Load
  x <- read.vcf(datasets$vcf[i], convert.chr = FALSE)
  x <- set.stats(x, verbose = FALSE)
  
  # Subset to window, retaining lead variants regardless of MAF
  keep <- x@snps$pos >= (WIN_START + offset) &
    x@snps$pos <= (WIN_END + offset) &
    (x@snps$maf >= MAF_MIN | x@snps$pos %in% pos)
  x <- x[, which(keep)]
  
  # Thin for legibility, retaining lead variants
  if (ncol(x) > MAX_SNPS) {
    idx  <- unique(round(seq(1, ncol(x), length.out = MAX_SNPS)))
    must <- which(x@snps$pos %in% pos)
    idx  <- sort(unique(c(idx, must)))
    x <- x[, idx]
  }
  
  # Axis labels
  lbl <- as.character(x@snps$pos)
  hit <- match(x@snps$pos, pos)
  lbl[!is.na(hit)] <- paste0("* ", lead_name[hit[!is.na(hit)]])
  x@snps$id <- lbl
  
  # LD matrix
  ld  <- LD(x, lim = c(1, ncol(x)), measure = "r2")
  idx <- match(pos, x@snps$pos)
  nm  <- lead_name[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  
  m <- round(ld[idx, idx], 3)
  dimnames(m) <- list(nm, nm)
  cat("\nr2 among the lead variants:\n")
  print(m)
  write.csv(m, sprintf("LD_CDKAL1_%s_r2.csv", label))
  
  # rs4710962 against each T2D lead
  v <- m["rs4710962", setdiff(nm, "rs4710962")]
  cat("\nr2 between rs4710962 and each T2D lead:\n")
  print(v)
  cat(sprintf("\nmax r2 = %.3f\n", max(v)))
  
  # Plot
  svglite::svglite(sprintf("LD_CDKAL1_%s.svg", label), width = 11, height = 9)
  LD.plot(ld,
          snp.positions = x@snps$pos,
          write.ld      = function(ld) if (ld > 0.02) round(ld, 2),
          write.snp.id  = TRUE,
          draw.chr      = TRUE,
          cex.ld        = 0.45,
          cex.snp       = 0.55,
          above.space   = 14,
          polygon.par   = list(border = "white"))
  dev.off()
}
