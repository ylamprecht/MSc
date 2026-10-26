# ============================================================
# LD Analysis: OAS1 Locus (chr12:112.91-112.93 Mb)
#
# Input:  OAS1_cases.vcf.gz
#         OAS1_controls.vcf.gz
#         OAS1_GBR.vcf.gz
#         OAS1_MEGA.vcf.gz
# Output: LD_OAS1_<dataset>.svg
#         LD_OAS1_<dataset>_r2.csv
# ============================================================

library(gaston)

# Datasets
datasets <- data.frame(
  label            = c("cases", "controls", "GBR", "MEGA"),
  vcf              = c("OAS1_cases.vcf.gz",
                       "OAS1_controls.vcf.gz",
                       "OAS1_GBR.vcf.gz",
                       "OAS1_MEGA.vcf.gz"),
  offset           = c(0, 0, 0, 437805),
  stringsAsFactors = FALSE
)

# Window and filtering
FLANK_LEFT  <- 10000
FLANK_RIGHT <- 5500
MAF_MIN     <- 0.05
MAX_SNPS    <- 50

# The nine OAS1 variants
oas1 <- data.frame(
  RSID  = c("rs10850094", "rs4767027", "rs1051042", "rs2660", "rs1131476",
            "rs4766664", "rs4767028", "rs10774671", "rs1131454"),
  pos38 = c(112922758, 112921352, 112919432, 112919637, 112919404,
            112925192, 112921383, 112919388, 112911065)
)

for (i in seq_len(nrow(datasets))) {
  
  label  <- datasets$label[i]
  offset <- datasets$offset[i]
  
  cat(sprintf("\n=== %s ===\n", label))
  
  # Load
  x <- read.vcf(datasets$vcf[i], convert.chr = FALSE, verbose = FALSE)
  x <- set.stats(x, verbose = FALSE)
  
  # Locate the OAS1 variants by rsID, falling back to position
  by_id    <- match(oas1$RSID, x@snps$id)
  by_pos   <- match(oas1$pos38 + offset, x@snps$pos)
  lead_idx <- ifelse(!is.na(by_id), by_id, by_pos)
  found    <- !is.na(lead_idx)
  
  cat(sprintf("OAS1 variants found: %d of %d\n", sum(found), nrow(oas1)))
  if (any(!found))
    cat("  not found: ", paste(oas1$RSID[!found], collapse = ", "), "\n", sep = "")
  
  # Subset to the locus, retaining OAS1 variants regardless of MAF
  lab <- rep("", ncol(x))
  lab[lead_idx[found]] <- paste0("* ", oas1$RSID[found])
  
  lead_pos <- x@snps$pos[lead_idx[found]]
  lo <- min(lead_pos) - FLANK_LEFT
  hi <- max(lead_pos) + FLANK_RIGHT
  
  keep <- which((x@snps$pos >= lo & x@snps$pos <= hi & x@snps$maf >= MAF_MIN) |
                  lab != "")
  x <- x[, keep]; lab <- lab[keep]
  
  # Thin for legibility, retaining OAS1 variants
  if (ncol(x) > MAX_SNPS) {
    idx <- unique(round(seq(1, ncol(x), length.out = MAX_SNPS)))
    idx <- sort(unique(c(idx, which(lab != ""))))
    x <- x[, idx]; lab <- lab[idx]
  }
  
  # Axis labels
  x@snps$id <- ifelse(lab != "", lab, as.character(x@snps$pos))
  
  # LD matrix
  ld <- LD(x, lim = c(1, ncol(x)), measure = "r2")
  li <- which(lab != "")
  m  <- round(ld[li, li], 3)
  nm <- sub("^\\* ", "", lab[li])
  dimnames(m) <- list(nm, nm)
  cat("\nr2 among the OAS1 variants:\n")
  print(m)
  write.csv(m, sprintf("LD_OAS1_%s_r2.csv", label))
  
  # rs10774671 against the other OAS1 variants
  key    <- "rs10774671"
  v      <- m[key, nm != key]
  others <- m[nm != key, nm != key]
  cat(sprintf("\nrs10774671 vs the other OAS1 variants: median r2 %.2f (range %.2f-%.2f)\n",
              median(v), min(v), max(v)))
  cat(sprintf("among the other OAS1 variants:          median r2 %.2f\n",
              median(others[upper.tri(others)])))
  
  # Plot
  svglite::svglite(sprintf("LD_OAS1_%s.svg", label),
                   width = 11, height = 11 * 0.5 + 2.4)
  LD.plot(ld,
          snp.positions = x@snps$pos,
          write.ld      = function(ld) if (ld > 0.02) round(ld, 2),
          write.snp.id  = TRUE,
          draw.chr      = TRUE,
          cex.ld        = 0.4,
          cex.snp       = 0.5,
          above.space   = 10,
          polygon.par   = list(border = "white"))
  dev.off()
}
