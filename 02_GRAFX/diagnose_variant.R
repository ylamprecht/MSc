# ============================================================
# Variant Diagnostics
#
# Input:  chr1_gt.tsv
#         samples.txt
#         PRF_meta_data.txt
# Output: console only
# ============================================================

library(data.table)
library(logistf)

# Settings
SNP          <- "chr1:48753712:G:A"
gt_file      <- "chr1_gt.tsv"
samples_file <- "samples.txt"
meta_file    <- "BETA_meta_data.txt"
covariates   <- c("age", "sex", "GBR", "YRI", "NAMA")

# Samples and metadata
samples <- trimws(readLines(samples_file))
samples <- samples[nzchar(samples)]

meta <- fread(meta_file, colClasses = list(character = "SampleID"))
meta[, SampleID := trimws(SampleID)]
meta <- meta[match(samples, meta$SampleID), ]
stopifnot(identical(samples, meta$SampleID))

# Genotypes
gt  <- fread(gt_file, header = FALSE)
row <- which(gt$V1 == SNP)
if (length(row) != 1) stop("expected exactly 1 row for ", SNP, ", found ", length(row))

g_chr <- as.character(unlist(gt[row, -1]))
geno  <- rep(0L, length(g_chr))
geno[g_chr == "0/1"] <- 1L
geno[g_chr == "1/1"] <- 2L
geno[g_chr == "./."] <- 0L
rm(gt)
invisible(gc())

dat <- data.frame(phenotype = meta$phenotype,
                  genotype  = geno,
                  meta[, ..covariates])

# Genotype by phenotype
cat("\n================ ", SNP, " ================\n", sep = "")
cat("\ngenotype x phenotype:\n")
print(table(genotype = dat$genotype, phenotype = dat$phenotype))

alt_case <- sum(dat$genotype[dat$phenotype == 1])
alt_ctrl <- sum(dat$genotype[dat$phenotype == 0])
n_case   <- sum(dat$phenotype == 1)
n_ctrl   <- sum(dat$phenotype == 0)
f1 <- alt_case / (2 * n_case)
f0 <- alt_ctrl / (2 * n_ctrl)
cat(sprintf("\nALT allele freq: cases %.3f | controls %.3f\n", f1, f0))
cat(sprintf("crude per-allele OR: %.2f\n", (f1 / (1 - f1)) / (f0 / (1 - f0))))

# Full model
full_form <- as.formula(paste("phenotype ~ genotype +",
                              paste(covariates, collapse = " + ")))
m <- logistf(full_form, data = dat, pl = TRUE,
             control = logistf.control(maxit = 1000))
i <- match("genotype", names(coef(m)))

cat(sprintf("\nFULL MODEL  OR = %.4g   95%% CI %.3g - %.3g   p = %.3g\n",
            exp(coef(m)[i]), exp(m$ci.lower[i]), exp(m$ci.upper[i]), m$prob[i]))

pr <- m$predict
cat(sprintf("fitted probabilities: min %.5f  max %.5f\n", min(pr), max(pr)))
cat(sprintf("pinned at 0 or 1 (<0.001 or >0.999): %d of %d  <-- quasi-separation if > 0\n",
            sum(pr < 0.001 | pr > 0.999), length(pr)))

# Dropping each term
cat("\n---- refit dropping each term ----\n")
cat(sprintf("%-14s %14s %14s %14s\n", "model", "OR", "CI low", "CI high"))

report <- function(label, f) {
  mm <- tryCatch(logistf(f, data = dat, pl = TRUE,
                         control = logistf.control(maxit = 1000)),
                 error = function(e) NULL)
  if (is.null(mm)) {
    cat(sprintf("%-14s %14s\n", label, "fit failed"))
    return(invisible(NULL))
  }
  k <- match("genotype", names(coef(mm)))
  cat(sprintf("%-14s %14.4g %14.4g %14.4g\n", label,
              exp(coef(mm)[k]), exp(mm$ci.lower[k]), exp(mm$ci.upper[k])))
}

report("full", full_form)
report("genotype only", as.formula("phenotype ~ genotype"))
for (term in covariates) {
  keep <- setdiff(covariates, term)
  f <- as.formula(paste("phenotype ~ genotype +", paste(keep, collapse = " + ")))
  report(paste0("drop ", term), f)
}

# Genotype and covariate correlation
cat("\n---- correlation of genotype with each covariate ----\n")
for (term in covariates) {
  cat(sprintf("  %-5s r = %+.3f\n", term, cor(dat$genotype, dat[[term]])))
}

cat("\nInterpretation:\n")
cat(" - OR collapses when one term is dropped  -> that term drives it\n")
cat(" - OR stays huge even for 'genotype only' -> the genotype vector itself\n")
cat(" - fitted probabilities pinned at 0/1     -> quasi-separation\n")
