#!/usr/bin/env bash
# Annotate existing 14,124-QC-passing site VCF using lab-managed reference resources.
# Requires locally maintained gnomAD v4.1 archive and/or vcfanno configuration.
# Never applies filters to remove rows; outputs annotation for inspection/reconciliation.
set -euo pipefail
ECH="${ECHTVAR_EXE:-echtvar}"
VCFANNO="${VCFANNO_EXE:-vcfanno}"
ANNOTATION_ROOT="${ANNOTATION_ROOT:-}"
LUA="${VCFANNO_LUA:-${ANNOTATION_ROOT:+$ANNOTATION_ROOT/scripts/annotation/vcfanno.lua}}"
RESOURCE_BASE="${VCFANNO_RESOURCE_BASE:-${ANNOTATION_ROOT:+$ANNOTATION_ROOT/resources}}"
BCF="${BCFTOOLS_EXE:-bcftools}"
SRC=''
ARCHIVE=''
CONFIG=''
OUT=''
usage(){
cat <<'USAGE'
Usage:
  bash annotate_local_gnomad_flags.sh --vcf FILTERED_LOF_SITE_VCF --outdir NEW_EMPTY_DIR \
      --echtvar-archive ABS_PATH_TO_GNOMAD_V4_ECHTVAR_ZIP \
      --vcfanno-config ABS_PATH_TO_GNOMAD_V4_VCFANNO_TOML

At least one annotation source is required. Source paths and INFO field names must
come from the local reference configuration.
The flag-extraction script then needs --filter-key EXACT_INFO_NAME (repeatable).

Does not contact gnomAD API, overwrite prior files, or use an include/exclude
expression: all records are preserved for QC reconciliation.
USAGE
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vcf) SRC="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --echtvar-archive) ARCHIVE="$2"; shift 2;;
    --vcfanno-config) CONFIG="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown argument: $1" >&2;usage;exit 2;;
  esac
done
[[ -n "$SRC" && -n "$OUT" && ( -n "$ARCHIVE" || -n "$CONFIG" ) ]] || { usage;exit 2; }
[[ -s "$SRC" ]] || { echo "Missing VCF: $SRC" >&2; exit 1; }
command -v "$BCF" >/dev/null 2>&1 || { echo "Missing bcftools: $BCF" >&2;exit 1; }
if [[ -n "$ARCHIVE" ]];then
  [[ -s "$ARCHIVE" ]] && command -v "$ECH" >/dev/null 2>&1 || { echo "Missing echtvar archive or executable" >&2;exit 1; }
fi
if [[ -n "$CONFIG" ]];then
  [[ -s "$CONFIG" && -s "$LUA" && -d "$RESOURCE_BASE" ]] && command -v "$VCFANNO" >/dev/null 2>&1 || { echo "Missing vcfanno config, executable, Lua helper or resource path" >&2;exit 1; }
fi
if [[ -e "$OUT" ]];then
  [[ -d "$OUT" && -z "$(ls -A "$OUT")" ]] || { echo "Refusing to use nonempty output location $OUT" >&2;exit 1; }
else mkdir -p "$OUT";fi
# Use variant-only VCF for faster execution and avoid producing patient-level data.
"$BCF" view -G -Ov -o "$OUT/query_site_only.vcf" "$SRC"
CURRENT="$OUT/query_site_only.vcf"
if [[ -n "$ARCHIVE" ]];then
  "$ECH" anno -e "$ARCHIVE" "$CURRENT" "$OUT/echtvar_annotated.bcf"
  "$BCF" view -Ov -o "$OUT/echtvar_annotated.vcf" "$OUT/echtvar_annotated.bcf"
  CURRENT="$OUT/echtvar_annotated.vcf"
fi
if [[ -n "$CONFIG" ]];then
  "$BCF" view -Oz -o "$OUT/vcfanno_input.vcf.gz" "$CURRENT"
  "$BCF" index -t "$OUT/vcfanno_input.vcf.gz"
  "$VCFANNO" -p 4 -base-path "$RESOURCE_BASE" -lua "$LUA" "$CONFIG" "$OUT/vcfanno_input.vcf.gz" > "$OUT/local_gnomad_annotations.vcf"
else
  cp "$CURRENT" "$OUT/local_gnomad_annotations.vcf"
fi
# Detect common representation problems prior to joining evidence by exact allele.
awk -F '\t' '$1 !~ /^#/ && $5 ~ /,/ {print "MULTI_ALT: " $1 ":" $2 > "/dev/stderr"; exit 1}' "$OUT/local_gnomad_annotations.vcf"
python3 "$(dirname "$0")/gnomad_flag_evidence_from_vcf.py" --vcf "$OUT/local_gnomad_annotations.vcf" --show-info-ids > "$OUT/INFO_FIELDS_TO_MAP.txt"
echo "Local gnomAD annotations saved: $OUT/local_gnomad_annotations.vcf"
echo "Review $OUT/INFO_FIELDS_TO_MAP.txt against the reference configuration."
echo "Next: python3 gnomad_flag_evidence_from_vcf.py --vcf '$OUT/local_gnomad_annotations.vcf' --filter-key gnomad_filter --output '$OUT/local_gnomad_flag_evidence.tsv'"
