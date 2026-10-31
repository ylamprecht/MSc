# ============================================================
# Literature-Derived Weighted Genetic Risk Score
#
# wGRS = sum over retained loci of ln(OR) x dosage of the
#        published associated allele
#
# Input:  published_risk_variants.csv
#         burden_calls_long.tsv
#         ld_overrides.tsv (optional)
# Output: wgrs_scores.tsv
#         wgrs_variant_audit.tsv
# ============================================================

library(data.table)

# Settings
LD_WINDOW_KB <- 250
R2_MAX       <- 0.50
CAP_BETA     <- NA
LD_COHORTS   <- c("UV-HE", "UV-MD", "UV-SD")
MIN_VARIANT_CALL <- 0.90
MIN_SAMPLE_CALL  <- 0.90
INDEL_INFERENCE  <- TRUE
INDEL_WINDOW     <- 25

SEVERITY_STATUS       <- c("severe COVID-19", "COVID-19 mortality", "ICU admission")
SUSCEPTIBILITY_STATUS <- c("infection susceptibility", "symptomatic COVID-19", "unknown")

NON_ADDITIVE <- c("rs17047200", "rs738409", "rs10741657", "rs12785878",
                  "rs34418657", "rs117011822", "rs7208524")

SCORE_SETS <- c("severity", "susceptibility")

# Input
v     <- fread("published_risk_variants.csv", sep = "\t", colClasses = "character")
calls <- fread("burden_calls_long.tsv")

audit <- data.table(stage = character(), reason = character(), n = integer())
log_drop <- function(stage, reason, n) if (n > 0)
  audit <<- rbind(audit, data.table(stage = stage, reason = reason, n = n))

# Duplicate records
dup_cohort <- calls[, .(n = uniqueN(cohort)), by = sample][n > 1]
if (nrow(dup_cohort)) {
  cat("WARNING: samples found in more than one cohort folder:\n"); print(dup_cohort)
  calls <- calls[calls[, .I[1], by = .(sample, variant)]$V1]
}
if (sum(duplicated(calls, by = c("variant", "sample")))) {
  cat("WARNING: duplicate variant/sample records collapsed\n")
  calls <- unique(calls, by = c("variant", "sample"))
}

# Score sets
unlisted <- setdiff(unique(v$Disease_status), c(SEVERITY_STATUS, SUSCEPTIBILITY_STATUS))
if (length(unlisted)) {
  cat("\nSTOP: disease-status values not in either predefined list:\n")
  print(v[Disease_status %in% unlisted, .N, by = Disease_status])
  cat("\nClassify each of these into SEVERITY_STATUS or SUSCEPTIBILITY_STATUS,\n",
      "or remove those variants, then rerun.\n", sep = "")
  quit(status = 1)
}
v[, score_set := fifelse(Disease_status %in% SEVERITY_STATUS, "severity", "susceptibility")]
missing_na <- setdiff(NON_ADDITIVE, v$RSID)
if (length(missing_na))
  cat("NOTE: listed as non-additive but not in the variant file: ",
      paste(missing_na, collapse = ", "), "\n", sep = "")

cat("phenotypes mapped to score sets:\n")
print(v[, .N, by = .(score_set, Disease_status)][order(score_set, -N)])

cat(sprintf("\nweighted scores will be produced for: %s\n",
            paste(SCORE_SETS, collapse = ", ")))
cat("allele counts are reported for every score set\n")

# Usability filters
n0 <- nrow(v)
f <- function(cond, reason) {
  k <- sum(cond); log_drop("filter", reason, k); v <<- v[!cond]
}
f(v$RSID %in% NON_ADDITIVE,              "genotype-model OR (not per-allele)")
f(v$Association == "not classified",     "no published direction")
f(v$Chromosome == "" | v$Position == "", "no position")
f(v$Chromosome == "X",                   "chromosome X")
f(v$Associated_allele == "",             "no associated allele")

v[, OR := suppressWarnings(as.numeric(Odds_ratio))]
f(is.na(v$OR),                           "no odds ratio")
f(v$OR <= 0 | v$OR == 1,                 "odds ratio <= 0 or exactly 1")

v[, beta  := log(OR)]
v[, p_pub := suppressWarnings(as.numeric(p_value))]
v[, pos   := as.integer(Position)]
v[, chrom := paste0("chr", Chromosome)]

# Direction check
bad <- v[(Association == "risk" & OR < 1) | (Association == "protective" & OR > 1)]
if (nrow(bad)) {
  cat("\nSTOP: published OR disagrees with the stored direction for these variants.\n")
  cat("The score is only valid if the OR refers to the associated allele.\n\n")
  print(bad[, .(RSID, Gene, Association, OR, p_value)])
  quit(status = 1)
}
cat("\ndirection check: all", nrow(v), "variants consistent with their OR\n")

# Allele harmonisation
calls[, c("chrom", "pos", "ref", "alt") := tstrsplit(variant, ":", fixed = TRUE)]
calls[, pos := as.integer(pos)]
recs <- unique(calls[, .(chrom, pos, ref, alt, variant)])

comp <- c(A = "T", T = "A", C = "G", G = "C")
flip <- function(x) if (all(strsplit(x, "")[[1]] %in% names(comp)))
  paste(comp[strsplit(x, "")[[1]]], collapse = "") else NA_character_

harmonise <- function(i) {
  vi <- v[i]
  a  <- vi$Associated_allele; r <- vi$Reference_allele
  none <- function(o) list(variant = NA_character_, orientation = o, method = NA_character_)

  is_snv <- vi$Mutation_type == "SNV" ||
            (nchar(a) == 1 && nchar(r) == 1 && a %in% names(comp) && r %in% names(comp))

  if (is_snv) {
    if (nchar(a) != 1 || nchar(r) != 1 || !a %in% names(comp) || !r %in% names(comp))
      return(none("allele mismatch"))
    if (comp[[a]] == r) return(none("palindromic"))
    cand <- recs[chrom == vi$chrom & pos == vi$pos & nchar(ref) == 1 & nchar(alt) == 1]
    if (!nrow(cand)) return(none("not observed"))
    fa <- flip(a)
    if (nrow(h <- cand[alt == a]))
      return(list(variant = h$variant[1], orientation = "ALT", method = "allele match"))
    if (nrow(h <- cand[ref == a]))
      return(list(variant = h$variant[1], orientation = "REF", method = "allele match"))
    if (!is.na(fa) && nrow(h <- cand[alt == fa]))
      return(list(variant = h$variant[1], orientation = "ALT", method = "allele match (complement)"))
    if (!is.na(fa) && nrow(h <- cand[ref == fa]))
      return(list(variant = h$variant[1], orientation = "REF", method = "allele match (complement)"))
    return(none("allele mismatch"))
  }

  exact <- recs[chrom == vi$chrom & pos == vi$pos]
  for (k in seq_len(nrow(exact))) {
    if (identical(exact$ref[k], r) && identical(exact$alt[k], a))
      return(list(variant = exact$variant[k], orientation = "ALT", method = "exact"))
    if (identical(exact$ref[k], a) && identical(exact$alt[k], r))
      return(list(variant = exact$variant[k], orientation = "REF", method = "exact"))
  }
  if (!INDEL_INFERENCE) return(none("indel unresolved"))

  alen <- function(x) if (x %in% c("", "-")) 0L else nchar(x)
  d <- alen(a) - alen(r)
  if (d == 0) return(none("indel unresolved"))

  wide <- recs[chrom == vi$chrom & abs(pos - vi$pos) <= INDEL_WINDOW &
               nchar(ref) != nchar(alt)]
  if (!nrow(wide)) return(none("not observed"))
  wide[, rdiff := nchar(alt) - nchar(ref)]
  wide <- wide[sign(rdiff) == sign(d)]
  if (!nrow(wide)) return(none("indel unresolved"))
  same_size <- wide[abs(rdiff) == abs(d)]
  how  <- if (nrow(same_size)) "inferred (length and size)" else "inferred (length direction)"
  if (nrow(same_size)) wide <- same_size
  hit <- wide[order(abs(pos - vi$pos), ref, alt)][1]
  list(variant = hit$variant, orientation = "ALT", method = how)
}

v <- cbind(v, rbindlist(lapply(seq_len(nrow(v)), harmonise), fill = TRUE))
for (o in c("not observed", "palindromic", "allele mismatch", "indel unresolved"))
  log_drop("harmonisation", o, sum(v$orientation == o))
v <- v[orientation %in% c("ALT", "REF") & !is.na(variant)]

cat("\nallele harmonisation:\n"); print(v[, .N, by = method][order(-N)])
inf <- v[grepl("^inferred", method)]
if (nrow(inf)) {
  cat(sprintf("\n%d indels oriented by allele length; check these individually:\n", nrow(inf)))
  print(inf[, .(RSID, Gene, Reference_allele, Associated_allele, variant, method,
                OR = round(OR, 2))])
}

if (any(duplicated(v$variant))) {
  d <- v$variant[duplicated(v$variant)]
  cat(sprintf("\nNOTE: %d RSIDs resolve to %d shared record(s); keeping the smallest p\n",
              sum(v$variant %in% d), uniqueN(d)))
  v[, abs_beta := abs(beta)]
  setorder(v, p_pub, -abs_beta, na.last = TRUE)
  log_drop("harmonisation", "duplicate record (same chrom:pos:ref:alt)", sum(duplicated(v$variant)))
  v <- unique(v, by = "variant")
}

# Genotype dosages
samples <- unique(calls[, .(sample, cohort)])

dose_alt <- c("0/0" = 0L, "0|0" = 0L,
              "0/1" = 1L, "1/0" = 1L, "0|1" = 1L, "1|0" = 1L,
              "1/1" = 2L, "1|1" = 2L,
              "./." = NA_integer_, ".|." = NA_integer_, "." = NA_integer_)

grid <- CJ(variant = v$variant, sample = samples$sample, unique = TRUE)
grid <- merge(grid, calls[, .(variant, sample, GT)], by = c("variant", "sample"), all.x = TRUE)
grid[, alt_dose := fifelse(is.na(GT), 0L, dose_alt[GT])]
if (any(!is.na(grid$GT) & is.na(grid$alt_dose) & !grid$GT %in% names(dose_alt)))
  cat("WARNING: unrecognised genotype strings were left missing\n")

cat(sprintf("\nexplicit no-call genotypes: %d of %d (%.2f%%)\n",
            sum(is.na(grid$alt_dose)), nrow(grid),
            100 * mean(is.na(grid$alt_dose))))

vcall <- grid[, .(call_rate = mean(!is.na(alt_dose))), by = variant]
lowv  <- vcall[call_rate < MIN_VARIANT_CALL, variant]
log_drop("missingness", sprintf("variant call rate < %.2f", MIN_VARIANT_CALL), length(lowv))
v    <- v[!variant %in% lowv]
grid <- grid[!variant %in% lowv]

# LD matrix
ld_samples <- samples[cohort %in% LD_COHORTS, sample]
cat(sprintf("\nLD estimated in %d unvaccinated individuals (%s)\n",
            length(ld_samples), paste(LD_COHORTS, collapse = ", ")))

gm   <- dcast(grid[sample %in% ld_samples], sample ~ variant, value.var = "alt_dose")
gmat <- as.matrix(gm[, -1]); rownames(gmat) <- gm$sample

# Published conditional overrides
ov <- if (file.exists("ld_overrides.tsv"))
  fread("ld_overrides.tsv", colClasses = "character") else
  data.table(RSID_A = character(), RSID_B = character(), relationship = character())
cat(sprintf("published conditional overrides loaded: %d\n", nrow(ov)))

rel <- function(a, b) {
  hit <- ov[(RSID_A == a & RSID_B == b) | (RSID_A == b & RSID_B == a)]
  if (nrow(hit)) hit$relationship[1] else NA_character_
}

# LD-aware clumping
clump <- function(vs) {
  keep <- character(); dropped <- 0L
  vs[, abs_beta := abs(beta)]
  setorder(vs, p_pub, -abs_beta, chrom, pos, na.last = TRUE)
  for (ch in unique(vs$chrom)) {
    w <- vs[chrom == ch]
    taken <- data.table(variant = character(), RSID = character(), pos = integer())
    for (i in seq_len(nrow(w))) {
      near <- taken[abs(pos - w$pos[i]) <= LD_WINDOW_KB * 1000]
      collapse <- FALSE
      for (j in seq_len(nrow(near))) {
        r <- rel(w$RSID[i], near$RSID[j])
        if (!is.na(r) && r == "independent")  next
        if (!is.na(r) && r == "same_signal") { collapse <- TRUE; break }
        r2 <- suppressWarnings(cor(gmat[, w$variant[i]], gmat[, near$variant[j]],
                                   use = "pairwise.complete.obs"))^2
        if (is.na(r2) || r2 < R2_MAX) next
        collapse <- TRUE; break
      }
      if (collapse) { dropped <- dropped + 1L } else {
        taken <- rbind(taken, data.table(variant = w$variant[i], RSID = w$RSID[i],
                                         pos = w$pos[i]))
        keep  <- c(keep, w$variant[i])
      }
    }
  }
  attr(keep, "dropped") <- dropped
  keep
}

kept <- lapply(split(v, v$score_set), clump)
for (s in names(kept)) log_drop(paste0("clumping (", s, ")"), "collapsed with a correlated variant",
                                attr(kept[[s]], "dropped"))
v    <- v[variant %in% unlist(kept)]
grid <- grid[variant %in% v$variant]

if (!is.na(CAP_BETA)) {
  cat(sprintf("\nSENSITIVITY ANALYSIS: |ln(OR)| capped at %.2f for %d loci\n",
              CAP_BETA, sum(abs(v$beta) > CAP_BETA)))
  v[, beta := pmax(pmin(beta, CAP_BETA), -CAP_BETA)]
}

# Scores
grid <- merge(grid, v[, .(variant, beta, orientation, score_set)], by = "variant")
grid[, dose := fifelse(orientation == "ALT", alt_dose, 2L - alt_dose)]
grid[, term := beta * dose]

grid <- merge(grid, v[, .(variant, Association)], by = "variant")

score <- grid[, .(wGRS      = sum(term, na.rm = TRUE),
                  n_risk_alleles       = sum(dose[Association == "risk"], na.rm = TRUE),
                  n_protective_alleles = sum(dose[Association == "protective"], na.rm = TRUE),
                  n_loci    = .N,
                  n_scored  = sum(!is.na(dose)),
                  n_missing = sum(is.na(dose))), by = .(sample, score_set)]
score[, call_rate := round(n_scored / n_loci, 4)]
score[!score_set %in% SCORE_SETS, wGRS := NA_real_]

low <- score[call_rate < MIN_SAMPLE_CALL]
if (nrow(low)) {
  cat(sprintf("\nWARNING: %d individual/score combinations below a call rate of %.2f\n",
              nrow(low), MIN_SAMPLE_CALL))
  print(low[, .(sample, score_set, n_missing, call_rate)])
}

wide <- dcast(score, sample ~ score_set,
              value.var = c("wGRS", "n_risk_alleles", "n_protective_alleles",
                            "n_loci", "n_missing", "call_rate"))
for (k in grep("^wGRS_", names(wide), value = TRUE))
  if (all(is.na(wide[[k]]))) wide[, (k) := NULL]
out  <- merge(samples, wide, by = "sample")
setnames(out, "sample", "SampleID")
for (k in grep("^wGRS_", names(out), value = TRUE))
  out[, (paste0(k, "_z")) := round(as.numeric(scale(get(k))), 4)]
for (k in grep("^wGRS_", names(out), value = TRUE)) out[, (k) := round(get(k), 4)]
setorder(out, cohort, SampleID)
fwrite(out, "wgrs_scores.tsv", sep = "\t")

# Variant audit
cat(sprintf("\n=== variant composition (from %d published variants) ===\n", n0))
print(audit)
cat("\nretained loci per score:\n")
print(v[, .(loci = .N, risk = sum(Association == "risk"),
            protective = sum(Association == "protective")), by = score_set])
if ("abs_beta" %in% names(v)) v[, abs_beta := NULL]
fwrite(v[, .(RSID, Gene, score_set, Disease_status, Association, OR, p_value,
             variant, orientation, method, beta = round(beta, 4))],
       "wgrs_variant_audit.tsv", sep = "\t")

# Sensitivity: five largest weights per score
cat("\n=== sensitivity: five largest |ln(OR)| within each score ===\n")
for (s in unique(v$score_set)) {
  top <- v[score_set == s][order(-abs(beta))][1:min(5, .N)]
  cat(sprintf("\n%s:\n", s)); print(top[, .(RSID, Gene, Association, OR, beta = round(beta, 2))])
  g <- merge(grid[score_set == s], samples, by.x = "sample", by.y = "sample")
  cat("mean score by cohort, all loci vs excluding those five:\n")
  print(g[, .(all_loci  = round(sum(term, na.rm = TRUE) / uniqueN(sample), 2),
              excl_top5 = round(sum(term[!variant %in% top$variant], na.rm = TRUE) /
                                uniqueN(sample), 2)), by = cohort])
}

# Summary
cat("\nallele counts per cohort:\n")
cnt <- grep("^n_(risk|protective)_alleles_", names(out), value = TRUE)
print(out[, lapply(.SD, function(x) round(mean(x), 1)), by = cohort, .SDcols = cnt])

cat(sprintf("\n%d individuals scored\n", nrow(out)))
for (k in grep("^wGRS_[a-z]+$", names(out), value = TRUE)) {
  cat(sprintf("\n%s:\n", k))
  print(out[, .(n = .N, mean = round(mean(get(k)), 2),
                median = round(as.double(median(get(k))), 2),
                IQR = sprintf("%.2f-%.2f", quantile(get(k), .25), quantile(get(k), .75))),
            by = cohort])
}
cat("\nwritten: wgrs_scores.tsv, wgrs_variant_audit.tsv\n")
