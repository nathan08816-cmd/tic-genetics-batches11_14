#!/usr/bin/env bash
# Safety-first entry point: no implicit reanalysis or overwriting.
set -euo pipefail
REPO="$(cd "$(dirname "$0")" && pwd)"
help() { cat <<'USAGE'
Usage:
  bash workflow.sh --show       List all stages and report their scope
  bash workflow.sh --check      Verify required paths from config/paths.env
  bash workflow.sh --clinvar    Run only the current ClinVar P/LP extraction (requires --check)
  bash workflow.sh --help

No automatic full pipeline run. Historic scripts require human review and
confirmed local paths. For ClinVar stage, FINAL_V41_AUDIT, FAMILY_EVENTS,
CLINVAR_OUTDIR must be set in the ignored config/paths.env.
USAGE
}
show() { cat <<'STAGES'
01  Site QC: Yale PASS + VQSR thresholds + Twist targets + split ALT + ANNOVAR
02  Predicted LoF/coding-splicing + genotype QC + affected/family evidence
03  Local gnomAD v4.1 echtvar; vcfanno ClinVar/region annotations
04  Current local v4.1 decision: PASS, REVIEW, EXCLUDE
05  Exact ClinVar Pathogenic/Likely_pathogenic and gene aggregation
06  Targeted trio source VCF QC (preliminary, NOT read-level confirmation)
07  Historical NDD overlap, recurrence and gene-level reviews
See docs/WORKFLOW.md for historical commands and paths requiring verification.
STAGES
}
load_config(){
  if [[ ! -f "$REPO/config/paths.env" ]]; then
    echo 'Missing config/paths.env; copy config/paths.env.example and fill verified paths.' >&2;return 2
  fi
  # shellcheck source=/dev/null
  source "$REPO/config/paths.env"
}
check_clinvar(){
  local x
  for x in FINAL_V41_AUDIT FAMILY_EVENTS CLINVAR_OUTDIR; do
    [[ -n "${!x:-}" ]] || { echo "Missing $x in config/paths.env" >&2; return 2; }
  done
  [[ -s "$FINAL_V41_AUDIT" ]] || { echo "Missing/empty audit: $FINAL_V41_AUDIT" >&2;return 2; }
  [[ -s "$FAMILY_EVENTS" ]] || { echo "Missing/empty family events: $FAMILY_EVENTS" >&2;return 2; }
  case "$CLINVAR_OUTDIR" in
    "$REPO"|"$REPO"/*) echo 'Output must be outside Git repository.' >&2;return 2;;
  esac
  [[ ! -e "$CLINVAR_OUTDIR" ]] || { echo 'Output directory already exists; refusing to overwrite.' >&2;return 2; }
  echo "ClinVar input files exist; output location is new: $CLINVAR_OUTDIR"
}
case "${1:---help}" in
  --help|-h)help;;
  --show)show;;
  --check) load_config;check_clinvar;;
  --clinvar)
    load_config;check_clinvar
    python3 "$REPO/scripts/05_clinvar_priority/filter_clinvar_plp_pass.py" \
      --audit "$FINAL_V41_AUDIT" --events "$FAMILY_EVENTS" --outdir "$CLINVAR_OUTDIR"
    ;;
  *) help;exit 2;;
esac
