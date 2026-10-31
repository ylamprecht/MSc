# ============================================================
# Firth Logistic Regression: Autosomes
#
# Model:  phenotype ~ genotype + age + sex + GBR + YRI + NAMA
#
# Input:  chr{1-22}_gt.tsv
#         samples.txt
#         meta_data.txt
# Output: chr{1-22}_firth_results.tsv
# Usage:  Rscript firthLR_perchr.R chr22
# ============================================================

library(data.table)
library(logistf)

# Chromosome
args <- commandArgs(trailingOnly = TRUE)
CHR  <- args[1]

# File paths
gt_file  <- paste0(CHR, "_gt.tsv")
out_file <- paste0(CHR, "_firth_results.tsv")

# Genotype dosage matrix
samples <- readLines("samples.txt")
gt <- fread(gt_file, header = FALSE)
G  <- as.matrix(gt[, -1])

D <- matrix(0L, nrow(G), ncol(G), dimnames = list(gt$V1, samples))
D[G == "0/1"] <- 1L
D[G == "1/1"] <- 2L
X <- t(D)

genotype_var <- colMeans(X^2) - colMeans(X)^2
X <- X[, genotype_var > 0, drop = FALSE]
cat(sprintf("%s: %d variants, %d testable after removing constant sites\n",
            CHR, length(genotype_var), ncol(X)))

# Metadata
meta <- fread("meta_data.txt")
meta$sex <- factor(meta$sex)
meta <- meta[match(rownames(X), meta$SampleID), ]

# Firth regression
is_case <- meta$phenotype == 1

fit_variant <- function(genotype) {
  carrier  <- genotype > 0
  car_case <- sum(carrier &  is_case)
  car_ctrl <- sum(carrier & !is_case)

  dat <- data.frame(phenotype = meta$phenotype,
                    genotype  = genotype,
                    age       = meta$age,
                    sex       = meta$sex,
                    GBR       = meta$GBR,
                    YRI       = meta$YRI,
                    NAMA      = meta$NAMA)

  m <- tryCatch(
    logistf(phenotype ~ genotype + age + sex + GBR + YRI + NAMA, data = dat),
    error = function(e) NULL)

  if (is.null(m))
    return(c(OR = NA_real_, p = NA_real_, car_case = car_case, car_ctrl = car_ctrl))

  i <- match("genotype", names(coef(m)))

  c(OR       = unname(exp(coef(m)[i])),
    p        = unname(m$prob[i]),
    car_case = car_case,
    car_ctrl = car_ctrl)
}

# Results
variant_fits <- t(apply(X, 2, fit_variant))

result <- data.table(SNP      = colnames(X),
                     OR       = variant_fits[, "OR"],
                     p        = variant_fits[, "p"],
                     car_case = as.integer(variant_fits[, "car_case"]),
                     car_ctrl = as.integer(variant_fits[, "car_ctrl"]))
result[, p_BH := p.adjust(p, "BH")]
setorder(result, p)
fwrite(result, out_file, sep = "\t")

cat(sprintf("%s: %d variants tested, %d could not be fitted; smallest p = %.3g\n",
            CHR, nrow(result), sum(is.na(result$p)), min(result$p, na.rm = TRUE)))
