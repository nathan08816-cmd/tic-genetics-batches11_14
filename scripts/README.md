# Scripts

Stages follow `README.md`. Original production scripts are listed in `docs/SOURCE_MAP.md`. Portable copies omit local names/paths and private trio IDs; confirm file paths and rerun before asserting exact reproduction. Historical scripts are for reference, not the current gnomAD v4.1 final filter.

`08_diagnostics/audit_annotation_summary.py` reports aggregate gnomAD,
ClinVar, and representation flags from the private final audit without printing
variant or sample IDs. Historical data files must not be committed.
