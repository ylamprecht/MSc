#!/usr/bin/env bash
# ============================================================
# Extract Calls at Published Risk Variant Positions
#
# Input:  published_risk_variants.csv
#         ../CVD_vcfs/<cohort>/*.vcf.gz
# Output: burden_regions.txt
#         burden_calls_long.tsv
# ============================================================
set -euo pipefail

# Settings
COHORTS=(UV-HE UV-MD UV-SD V-MD)
VARIANTS="published_risk_variants.csv"
REGIONS="burden_regions.txt"
OUT="burden_calls_long.tsv"

# Query regions
awk -F'\t' 'NR>1 && $3 != "" && $3 != "X" && $4 != "" && $8 != "not classified" {
              print "chr"$3"\t"$4 }' "$VARIANTS" | sort -k1,1 -k2,2n -u > "$REGIONS"
echo "positions queried: $(wc -l < "$REGIONS")"

# Genotype calls
printf 'variant\tsample\tcohort\tGT\n' > "$OUT"

n=0
for cohort in "${COHORTS[@]}"; do
  dir="../CVD_vcfs/$cohort"
  [ -d "$dir" ] || { echo "MISSING DIR: $dir" >&2; continue; }

  for vcf in "$dir"/*.vcf.gz; do
    [ -e "$vcf" ] || continue
    s=$(bcftools query -l "$vcf" | head -1)
    [ -z "$s" ] && { echo "WARNING: no sample name in $vcf" >&2; continue; }
    [ -e "${vcf}.csi" ] || [ -e "${vcf}.tbi" ] || bcftools index "$vcf"

    bcftools query -R "$REGIONS" -f '%CHROM:%POS:%REF:%ALT\t[%GT]\n' "$vcf" \
      | awk -v s="$s" -v c="$cohort" 'BEGIN{OFS="\t"} {print $1, s, c, $2}' >> "$OUT"
    n=$((n+1))
  done
  echo "  $cohort: done"
done

echo "samples processed: $n"
echo "call records: $(( $(wc -l < "$OUT") - 1 ))"
