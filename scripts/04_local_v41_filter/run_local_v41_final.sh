#!/usr/bin/env bash
# Re-finalize from the PREVIOUS COMPLETED run; no reannotation, and no old API calls.
set -euo pipefail
ROOT="${PROJECT_ROOT:?Set PROJECT_ROOT to the private TS_ANNOVAR_analysis directory}"
SUB="$ROOT/filtering/family_batches11_14"
HERE="$(cd "$(dirname "$0")" && pwd)"
SOURCE="${SOURCE_RUN:-}"
if [[ -z "$SOURCE" ]];then
  while IFS= read -r d;do
    if [[ -s "$d/final_analysis/ALL_14124_COMPLETE_AUDIT.tsv" && -s "$d/local_annotation/local_gnomad_flag_evidence.tsv" ]];then
      SOURCE="$d";break
    fi
  done < <(find "$SUB" -maxdepth 1 -type d -name 'final_jointcall_20*' | sort -r)
fi
[[ -n "$SOURCE" ]] || { echo 'ERROR: Could not find a completed source run. Set SOURCE_RUN to its full directory.' >&2;exit 2; }
[[ -s "$SOURCE/final_analysis/ALL_14124_COMPLETE_AUDIT.tsv" ]] || { echo 'Missing prior full audit' >&2;exit 2; }
[[ -s "$SOURCE/local_annotation/local_gnomad_flag_evidence.tsv" ]] || { echo 'Missing prior echtvar/vcfanno evidence' >&2;exit 2; }
OUT="$SUB/local_v41_final_$(date +%Y%m%d_%H%M%S)"
[[ ! -e "$OUT" ]] || { echo "Output path exists: $OUT" >&2;exit 2; }
VCF="$SUB/lof_joint_20261026/results/batch1114_lof_QCsite_candidates.with_genotypes.vcf.gz"
ORIG="$SUB/scripts/02_lof_family_screen.py"
QC="${SAMPLE_QC:?Set SAMPLE_QC to the local sample QC table}"
PHENO="${PHENOTYPE_XLSX:?Set PHENOTYPE_XLSX to the private phenotype workbook}"
for f in "$VCF" "$ORIG" "$QC" "$PHENO";do [[ -s "$f" ]] || { echo "Missing required file: $f" >&2;exit 2; }; done
printf 'SOURCE RUN: %s\nOUTPUT: %s\n' "$SOURCE" "$OUT"
python3 -u "$HERE/refinalize_local_gnomad_v41.py" \
 --audit "$SOURCE/final_analysis/ALL_14124_COMPLETE_AUDIT.tsv" \
 --local-evidence "$SOURCE/local_annotation/local_gnomad_flag_evidence.tsv" \
 --vcf "$VCF" --qc-code "$ORIG" --sample-qc "$QC" --phenotype "$PHENO" \
 --family-summary "$SOURCE/final_analysis/family_summary.tsv" \
 --outdir "$OUT"
echo "DONE: $OUT"
