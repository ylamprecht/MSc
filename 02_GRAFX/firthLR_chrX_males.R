# ============================================================
# Firth Logistic Regression: Chromosome X, Males
#
# Model:  phenotype ~ genotype + age + GBR + YRI + NAMA
#
# Input:  chrX_male_gt.tsv
#         samples_male.txt
#         meta_males.txt
# Output: chrX_male_firth_results.tsv
# ============================================================

library(data.table)
library(logistf)

# File paths
gt_file  <- "chrX_male_gt.tsv"
out_file <- "chrX_male_firth_results.tsv"

# Filtering
HEMI_MIN <- 0.9

# Genotype dosage matrix
samples <- readLines("samples_male.txt")
gt <- fread(gt_file, header = FALSE)
G  <- as.matrix(gt[, -1])

D <- matrix(0L, nrow(G), ncol(G), dimnames = list(gt$V1, samples))
D[G == "0/0"] <- 0L
D[G == "0/1"] <- NA_integer_
D[G == "1/1"] <- 2L
D[G == "./."] <- 0L
X <- t(D)

cat(sprintf("chrX (males): %d heterozygous calls set to NA\n", sum(is.na(X))))

hemizygous_rate <- colMeans(!is.na(X))
genotype_var    <- colMeans(X^2, na.rm = TRUE) - colMeans(X, na.rm = TRUE)^2
keep <- !is.na(genotype_var) & genotype_var > 0 & hemizygous_rate >= HEMI_MIN

cat(sprintf("chrX (males): %d variants below the %.1f hemizygous-call threshold\n",
            sum(hemizygous_rate < HEMI_MIN), HEMI_MIN))
print(table(cut(hemizygous_rate, c(0, 0.5, 0.8, 0.9, 0.95, 0.99, 1))))

X <- X[, keep, drop = FALSE]
cat(sprintf("chrX (males): %d variants, %d testable after filtering\n",
            length(genotype_var), ncol(X)))

# Metadata
meta <- fread("meta_males.txt")
meta <- meta[match(rownames(X), meta$SampleID), ]

# Firth regression
is_case <- meta$phenotype == 1

fit_variant <- function(genotype) {
  carrier  <- genotype > 0
  car_case <- sum(carrier &  is_case, na.rm = TRUE)
  car_ctrl <- sum(carrier & !is_case, na.rm = TRUE)

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
             N = NA_real_, car_case = car_case, car_ctrl = car_ctrl))

  i <- match("genotype", names(coef(m)))

  c(OR       = unname(exp(coef(m)[i])),
    OR_95L   = unname(exp(m$ci.lower[i])),
    OR_95U   = unname(exp(m$ci.upper[i])),
    p        = unname(m$prob[i]),
    N        = m$n,
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
                     N        = as.integer(variant_fits[, "N"]),
                     EAF      = unname(colMeans(X > 0, na.rm = TRUE)[colnames(X)]),
                     car_case = as.integer(variant_fits[, "car_case"]),
                     car_ctrl = as.integer(variant_fits[, "car_ctrl"]))

n_nofit <- sum(is.na(result$p))
result  <- result[!is.na(p)]
setorder(result, p, na.last = TRUE)
fwrite(result, out_file, sep = "\t")

# Diagnostics
n_badci <- sum(!is.finite(result$OR_95L) | !is.finite(result$OR_95U) |
               result$OR_95L <= 0)

cat(sprintf("chrX (males): %d variants written, %d could not be fitted and were dropped; smallest p = %.3g\n",
            nrow(result), n_nofit, min(result$p)))
cat(sprintf("chrX (males): per-variant N ranges %d - %d\n",
            min(result$N), max(result$N)))
if (sum(duplicated(result$SNP)) > 0)
  cat(sprintf("WARNING: %d duplicated marker names - GWAMA requires unique IDs\n",
              sum(duplicated(result$SNP))))
if (n_badci > 0)
  cat(sprintf("WARNING: %d variants have an unbounded CI - GWAMA will reject these\n",
              n_badci))
