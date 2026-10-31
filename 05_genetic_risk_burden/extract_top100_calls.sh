#!/usr/bin/env bash
# ============================================================
# Extract Calls at Top 100 Association Positions
#
# Input:  top100_annotated.tsv
#         ../CVD_vcfs/<cohort>/*.vcf.gz
# Output: top100_regions.txt
#         top100_calls_long.tsv
# ============================================================
set -euo pipefail

# Settings
COHORTS=(UV-HE UV-MD UV-SD V-MD)
VARIANTS="top100_annotated.tsv"
REGIONS="top100_regions.txt"
OUT="top100_calls_long.tsv"

# Query regions
awk -F'\t' 'NR>1 {split($1,a,":"); print a[1]"\t"a[2]}' "$VARIANTS" \
  | sort -k1,1 -k2,2n -u > "$REGIONS"
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
