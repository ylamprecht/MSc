# ============================================================
# Ancestry Proportions by Sex
#
# Input:  meta_data.txt
# Output: printed summary and Wilcoxon test results
# ============================================================

library(data.table)

# Ancestry components
comps <- c("GBR", "YRI", "NAMA", "ASIAN")

# Load metadata
meta <- fread("meta_data.txt")
meta <- meta[SampleID != "CVD719"]
meta$sex <- factor(meta$sex, levels = c(0, 1), labels = c("female", "male"))

# Counts
cat("males:", sum(meta$sex == "male"),
    " females:", sum(meta$sex == "female"), "\n")

# Mean ancestry by sex
print(meta[, lapply(.SD, function(x) round(mean(x), 3)), by = sex, .SDcols = comps])

# Wilcoxon test per component
results <- rbindlist(lapply(comps, function(comp) {
  w <- wilcox.test(meta[[comp]] ~ meta$sex)
  data.table(component     = comp,
             median_female = round(median(meta[sex == "female"][[comp]]), 6),
             median_male   = round(median(meta[sex == "male"][[comp]]), 6),
             p_value       = round(w$p.value, 4))
}))
print(results)
