#!/usr/bin/env bash
set -euo pipefail

# Batch11–14 quality and ANNOVAR preparation (non-destructive).
# Usage: bash 01_site_qc_and_annovar.sh qc
#        ANNOVAR_DIR=/path/to/annovar HUMANDB=/path/to/humandb bash 01_site_qc_and_annovar.sh annotate
MODE="${1:-}"
PROJECT="${PROJECT:-${PROJECT_ROOT:-}}"
DATA="${DATA:-${SOURCE_DATA_DIR:-}}"
[[ -n "$PROJECT" && -n "$DATA" ]] || { echo "Set PROJECT_ROOT and SOURCE_DATA_DIR." >&2; exit 2; }
EXTRACT="$PROJECT/filtering/family_batches11_14/targets_twist20/exome_calls.twist20.all_filters.vcf.gz"
OUT="$PROJECT/filtering/family_batches11_14/lof_joint_20261026"
VQ="$OUT/01_vqsr997_950.all_ycga_filters.vcf.gz"
STRICT="$OUT/02_vqsr997_950.originalPASS.vcf.gz"
NORM="$OUT/03_vqsr997_950.originalPASS.splitALT.vcf.gz"
ANNO_INPUT="$OUT/04_for_annovar.one_sample.vcf"
ANNO_PREFIX="$OUT/04_refgene"
LOGDIR="$OUT/logs"
mkdir -p "$LOGDIR"

need() { test -s "$1" || { echo "Missing or empty input: $1" >&2; exit 1; }; }
refuse() { if test -e "$1"; then echo "REFUSING overwrite: $1" >&2; exit 1; fi; }
sizecheck() {
  # Do not proceed on congested shared storage
  local avail
  avail=$(df -Pk "$PROJECT" | awk 'NR==2 {print $4}')
  test "$avail" -ge 104857600 || { echo 'Less than 100 GiB free, stop for lab review' >&2; exit 1; }
}
count() { bcftools view -H "$1" | wc -l; }

if [[ "$MODE" == qc ]]; then
  need "$EXTRACT"
  need "$DATA/filtering/vqsr_only.sites.vcf.gz"
  need "$DATA/filtering/vqsr_only.sites.vcf.gz.tbi"
  command -v bcftools >/dev/null
  for f in "$VQ" "$STRICT" "$NORM"; do refuse "$f"; refuse "$f.tbi"; done
  sizecheck
  echo "[1/3] Annotate VQSR scoring model, retain SNP TS=99.7% and indel TS=95.0%."
  echo "      Preserve original FILTER field: VQSR-only selection is NOT final QC."
  bcftools annotate \
      -a "$DATA/filtering/vqsr_only.sites.vcf.gz" \
      -c INFO/VQSR_MODEL "$EXTRACT" -Ou \
    | bcftools view \
      -i '(INFO/VQSR_MODEL="SNP" && INFO/VQSLOD>=-3.4427) || (INFO/VQSR_MODEL="INDEL" && INFO/VQSLOD>=2.4465)' \
      -Oz -o "$VQ"
  bcftools index -t "$VQ"
  echo "[2/3] Intersection with Yale's original FILTER=PASS for conservative main analysis."
  bcftools view -f PASS "$VQ" -Oz -o "$STRICT"
  bcftools index -t "$STRICT"
  echo "[3/3] Split multiallelic calls into ALT-specific records. No reference normalization."
  bcftools norm -m -any "$STRICT" -Oz -o "$NORM"
  bcftools index -t "$NORM"
  echo "==== VALIDATION / COUNTS ===="
  for f in "$EXTRACT" "$VQ" "$STRICT" "$NORM"; do
    echo "$(basename "$f")  sites=$(count "$f")  samples=$(bcftools query -l "$f" | wc -l)"
  done
  echo "== Number of originally LowQual or otherwise non-PASS calls among VQSR-retained =="
  bcftools view -H "$VQ" | awk -F '\t' '$7!="PASS" && $7!="." {n++} END {print n+0}'
  echo 'QC_STAGE_SUCCESS (not yet genotype-screened or annotated)'

elif [[ "$MODE" == annotate ]]; then
  need "$NORM"
  [[ -n "${ANNOVAR_DIR:-}" && -n "${HUMANDB:-}" ]] || {
    echo "Set ANNOVAR_DIR=/path/to/annovar and HUMANDB=/path/to/humandb." >&2; exit 1;
  }
  need "$ANNOVAR_DIR/table_annovar.pl"
  need "$HUMANDB/hg38_refGene.txt"
  need "$HUMANDB/hg38_refGeneMrna.fa"
  refuse "$ANNO_INPUT"
  refuse "$ANNO_PREFIX.hg38_multianno.vcf"
  refuse "$ANNO_PREFIX.hg38_multianno.txt"
  sizecheck
  echo '[1/2] Generate a sites-only VCF with one explicitly ARTIFICIAL carrier genotype for ANNOVAR.'
  # Do not choose a real study sample: rare sites may have 0/0 in that sample,
  # which can cause some ANNOVAR input conversion modes to omit their ALT alleles.
  # This synthetic ALT genotype is used only for obtaining annotations and is
  # NEVER included in the full 1,245-sample analysis or family screen.
  bcftools view -G -Ov "$NORM" | awk 'BEGIN{OFS="\t"} /^#CHROM/ {print $0,"FORMAT","ANNOTATION_ONLY_SYNTHETIC_ALT"; next} /^#/ {print; next} {print $0,"GT","0/1"}' > "$ANNO_INPUT"
  echo '[2/2] ANNOVAR refGene annotation of each ALT allele (synthetic annotation-only sample).'
  expected_sites=$(grep -vc '^#' "$ANNO_INPUT")
  echo "Expected annotation input sites: $expected_sites"
  # ANNOVAR's VCF output retains the VCF coordinate/REF/ALT representation.
  (cd "$ANNOVAR_DIR" && perl table_annovar.pl \
     "$ANNO_INPUT" "$HUMANDB" \
     -buildver hg38 -out "$ANNO_PREFIX" \
     -protocol refGene -operation g -nastring . -vcfinput) \
     > "$LOGDIR/annovar.stdout.log" 2> "$LOGDIR/annovar.stderr.log"
  need "$ANNO_PREFIX.hg38_multianno.vcf"
  need "$ANNO_PREFIX.hg38_multianno.txt"
  actual_sites=$(grep -vc '^#' "$ANNO_PREFIX.hg38_multianno.vcf")
  echo "ANNOVAR output variants: $actual_sites"
  if [[ "$actual_sites" -ne "$expected_sites" ]]; then
    echo "ANNOTATION COUNT MISMATCH: input $expected_sites vs annotated $actual_sites; stop and investigate." >&2
    exit 1
  fi
  echo 'ANNOVAR_STAGE_SUCCESS' 
  echo "Annotated (synthetic-only-sample) VCF: $ANNO_PREFIX.hg38_multianno.vcf"
else
  echo 'Usage: bash 01_site_qc_and_annovar.sh qc|annotate' >&2; exit 2
fi
