# ============================================================
# Firth Logistic Regression: Chromosome X, Females
#
# Model:  phenotype ~ genotype + age + GBR + YRI + NAMA
#
# Input:  chrX_female_gt.tsv
#         samples_female.txt
#         meta_females.txt
# Output: chrX_female_firth_results.tsv
# ============================================================

library(data.table)
library(logistf)

# File paths
gt_file  <- "chrX_female_gt.tsv"
out_file <- "chrX_female_firth_results.tsv"

# Genotype dosage matrix
samples <- readLines("samples_female.txt")
gt <- fread(gt_file, header = FALSE)
G  <- as.matrix(gt[, -1])

D <- matrix(0L, nrow(G), ncol(G), dimnames = list(gt$V1, samples))
D[G == "0/0"] <- 0L
D[G == "0/1"] <- 1L
D[G == "1/1"] <- 2L
D[G == "./."] <- 0L
X <- t(D)

genotype_var <- colMeans(X^2) - colMeans(X)^2
X <- X[, genotype_var > 0, drop = FALSE]
cat(sprintf("chrX (females): %d variants, %d testable after removing constant sites\n",
            length(genotype_var), ncol(X)))

# Metadata
meta <- fread("meta_females.txt")
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
                    GBR       = meta$GBR,
                    YRI       = meta$YRI,
                    NAMA      = meta$NAMA)

  m <- tryCatch(
    logistf(phenotype ~ genotype + age + GBR + YRI + NAMA,
            data = dat, pl = TRUE, plconf = 2),
    error = function(e) NULL)

  if (is.null(m))
    return(c(OR = NA_real_, OR_95L = NA_real_, OR_95U = NA_real_, p = NA_real_,
             car_case = car_case, car_ctrl = car_ctrl))

  i <- match("genotype", names(coef(m)))

  c(OR       = unname(exp(coef(m)[i])),
    OR_95L   = unname(exp(m$ci.lower[i])),
    OR_95U   = unname(exp(m$ci.upper[i])),
    p        = unname(m$prob[i]),
    car_case = car_case,
    car_ctrl = car_ctrl)
}

# Results
variant_fits <- t(apply(X, 2, fit_variant))
snp_parts    <- tstrsplit(colnames(X), ":", fixed = TRUE)

result <- data.table(SNP      = colnames(X),
                     EA       = snp_parts[[4]],
                     NEA      = snp_parts[[3]],
                     OR       = variant_fits[, "OR"],
                     OR_95L   = variant_fits[, "OR_95L"],
                     OR_95U   = variant_fits[, "OR_95U"],
                     p        = variant_fits[, "p"],
                     N        = nrow(meta),
                     EAF      = unname(colMeans(X)[colnames(X)] / 2),
                     car_case = as.integer(variant_fits[, "car_case"]),
                     car_ctrl = as.integer(variant_fits[, "car_ctrl"]))
setorder(result, p)
fwrite(result, out_file, sep = "\t")

# Diagnostics
n_badci <- sum(!is.na(result$p) &
               (!is.finite(result$OR_95L) | !is.finite(result$OR_95U) |
                result$OR_95L <= 0))

cat(sprintf("chrX (females): %d variants tested, %d could not be fitted; smallest p = %.3g\n",
            nrow(result), sum(is.na(result$p)), min(result$p, na.rm = TRUE)))
if (sum(duplicated(result$SNP)) > 0)
  cat(sprintf("WARNING: %d duplicated marker names - GWAMA requires unique IDs\n",
              sum(duplicated(result$SNP))))
if (n_badci > 0)
  cat(sprintf("WARNING: %d variants have an unbounded CI - GWAMA will reject these\n",
              n_badci))
