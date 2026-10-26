# ============================================================
# Firth Results for Published Risk Variants
#
# Input:  published_risk_variants.csv
#         chr{1-22}_firth_results.tsv
#         chrX_combined_gwama_results.tsv
# Output: published_risk_variants_results.tsv
#         not_concordant.tsv
# ============================================================

library(data.table)

# File paths
in_file   <- "published_risk_variants.csv"
out_file  <- "published_risk_variants_results.tsv"
chrX_file <- "chrX_combined_gwama_results.tsv"

# Variant set
vd <- fread(in_file, sep = "\t", colClasses = "character", na.strings = NULL)
vd[, row_id := .I]
vd[, pos := suppressWarnings(as.integer(Position))]

# Alleles
comp <- c(A = "T", T = "A", C = "G", G = "C")
palindromic <- function(a, b) isTRUE(comp[a] == b)
alen <- function(x) if (x %in% c("", "-")) 0L else nchar(x)

# Association results near the required positions
load_chr <- function(f, cols, positions) {
  r <- fread(f, select = cols)
  r[, c("rchr", "rpos", "ref", "alt") := tstrsplit(SNP, ":", fixed = TRUE)]
  r[, rpos := as.integer(rpos)]
  r[rpos %in% unique(unlist(lapply(positions, function(p) (p - 30):(p + 30))))]
}

res <- list()
for (ch in unique(vd$Chromosome[!is.na(vd$pos)])) {
  p <- vd[Chromosome == ch & !is.na(pos), pos]
  if (ch == "X") {
    if (!file.exists(chrX_file)) next
    res[[ch]] <- load_chr(chrX_file, c("SNP", "p", "female_OR", "male_OR"), p)
  } else {
    f <- sprintf("chr%s_firth_results.tsv", ch)
    if (!file.exists(f)) { cat("missing results file:", f, "\n"); next }
    res[[ch]] <- load_chr(f, c("SNP", "OR", "p", "car_case", "car_ctrl"), p)
  }
}

# Matching and orientation
match_one <- function(v) {
  out <- list(res_SNP     = NA_character_,
              OR          = NA_real_,
              p           = NA_real_,
              car_case    = NA_integer_,
              car_ctrl    = NA_integer_,
              orientation = NA_character_,
              OR_assoc    = NA_real_,
              note        = "")
  
  if (is.na(v$pos)) { out$note <- "no position in source"; return(out) }
  r <- res[[v$Chromosome]]
  if (is.null(r))   { out$note <- "no results file";       return(out) }
  
  ra <- v$Reference_allele
  a  <- v$Associated_allele
  is_snv <- v$Mutation_type == "SNV" ||
    (nchar(ra) == 1 && nchar(a) == 1 && ra %in% names(comp) && a %in% names(comp))
  oriented <- FALSE
  flip     <- FALSE
  
  if (is_snv) {
    cand <- r[rpos == v$pos & nchar(ref) == 1 & nchar(alt) == 1]
    if (!nrow(cand)) {
      out$note <- "not in results (monomorphic, not called or QC-removed)"
      return(out)
    }
    if (a == "") {
      hit <- cand[1]
      out$orientation <- "no associated allele given"
    } else if (nrow(h <- cand[alt == a])) {
      hit <- h[1]
      out$orientation <- "ALT"
      oriented <- TRUE
    } else if (nrow(h <- cand[ref == a])) {
      hit <- h[1]
      out$orientation <- "REF (OR inverted)"
      oriented <- TRUE
      flip     <- TRUE
    } else if (!is.na(comp[a]) && palindromic(ra, a)) {
      hit <- cand[1]
      out$orientation <- "palindromic, strand ambiguous"
    } else if (!is.na(comp[a]) && nrow(h <- cand[alt == comp[a]])) {
      hit <- h[1]
      out$orientation <- "ALT (opposite strand)"
      oriented <- TRUE
    } else if (!is.na(comp[a]) && nrow(h <- cand[ref == comp[a]])) {
      hit <- h[1]
      out$orientation <- "REF (opposite strand, OR inverted)"
      oriented <- TRUE
      flip     <- TRUE
    } else {
      hit <- cand[1]
      out$orientation <- "allele mismatch"
    }
    if (nrow(cand) > 1) out$note <- sprintf("%d records at this position", nrow(cand))
  } else {
    w    <- max(2, nchar(ra)) + 2
    cand <- r[abs(rpos - v$pos) <= w & nchar(ref) != nchar(alt)]
    if (!nrow(cand)) {
      out$note <- "indel not in results within window"
      return(out)
    }
    cand[, rdiff := nchar(alt) - nchar(ref)]
    d <- if (a == "") NA_integer_ else alen(a) - alen(ra)
    if (!is.na(d) && d != 0 && nrow(same <- cand[abs(rdiff) == abs(d)])) cand <- same
    hit <- cand[order(abs(rpos - v$pos))][1]
    if (is.na(d) || d == 0) {
      out$orientation <- "indel: associated allele not resolvable"
    } else if (sign(hit$rdiff) == sign(d)) {
      out$orientation <- "indel ALT (by allele length)"
      oriented <- TRUE
    } else {
      out$orientation <- "indel REF (by allele length, OR inverted)"
      oriented <- TRUE
      flip     <- TRUE
    }
    out$note <- sprintf("window match at %+d bp; %d candidate record(s)",
                        hit$rpos - v$pos, nrow(cand))
  }
  
  out$res_SNP <- hit$SNP
  out$p       <- hit$p
  
  if (v$Chromosome == "X") {
    fo <- if (flip) 1 / hit$female_OR else hit$female_OR
    mo <- if (flip) 1 / hit$male_OR   else hit$male_OR
    out$note <- trimws(paste(out$note, sprintf(
      "chrX sex-differentiated p; female_OR=%.3g, male_OR=%.3g %s", fo, mo,
      if (oriented) "(per associated allele)" else "(per ALT)")))
  } else {
    out$OR       <- hit$OR
    out$car_case <- as.integer(hit$car_case)
    out$car_ctrl <- as.integer(hit$car_ctrl)
    if (oriented) out$OR_assoc <- if (flip) 1 / hit$OR else hit$OR
  }
  out
}

added <- rbindlist(lapply(seq_len(nrow(vd)), function(i) match_one(vd[i])))
out   <- cbind(vd[, !c("row_id", "pos")], added)

# Multiple testing and concordance
out[!is.na(p), q_candidate := p.adjust(p, method = "BH")]

out[, expected := fcase(Association == "risk",        1L,
                        Association == "protective", -1L,
                        default = NA_integer_)]
out[, concordant := fifelse(is.na(OR_assoc) | is.na(expected), NA,
                            sign(log(OR_assoc)) == expected)]
out[, expected := NULL]
fwrite(out, out_file, sep = "\t")

# Report
cat(sprintf("%d of %d variants present in this study\n",
            sum(!is.na(out$p)), nrow(out)))

cc <- out[!is.na(concordant)]
cat(sprintf("\nconcordance: %d / %d (%.1f%%)\n",
            sum(cc$concordant), nrow(cc), 100 * mean(cc$concordant)))

disc <- cc[concordant == FALSE,
           .(RSID, Gene, Association, Associated_allele,
             OR_assoc = round(OR_assoc, 2), p = signif(p, 3))][order(Gene, RSID)]
cat(sprintf("\nnot concordant (%d):\n", nrow(disc)))
print(disc, nrows = Inf)
fwrite(disc, "not_concordant.tsv", sep = "\t")

cat(sprintf("\nresults: %s | not concordant list: not_concordant.tsv\n", out_file))