# Source verification

Source script checksums before portability/privacy adjustments. File paths are relative to the private project root.

| Server source | Repository copy | Source SHA-256 prefix |
|---|---|---|
| `JointCall_Configured_Local_gnomAD_20261110/annotate_local_gnomad_flags.sh` | `scripts/03_local_annotation/annotate_local_gnomad_flags.sh` | `7c2e51a9ea14c023` |
| `JointCall_Configured_Local_gnomAD_20261110/gnomad_flag_evidence_from_vcf.py` | `scripts/03_local_annotation/gnomad_flag_evidence_from_vcf.py` | `1c97c4c6ad9ed333` |
| `JointCall_LocalV41_Final_20261110/refinalize_local_gnomad_v41.py` | `scripts/04_local_v41_filter/refinalize_local_gnomad_v41.py` | `6ab7a4f6f637deb3` |
| `JointCall_LocalV41_Final_20261110/run_local_v41_final.sh` | `scripts/04_local_v41_filter/run_local_v41_final.sh` | `3401e225967fb0f9` |
| `archive/legacy_root_20261110_171554/batch1114_denovo_priority.py` | `scripts/07_gene_analysis_historical/batch1114_denovo_priority.py` | `2616e12485509738` |
| `archive/legacy_root_20261110_171554/batch1114_full_jointcall_gene_review.py` | `scripts/07_gene_analysis_historical/batch1114_full_jointcall_gene_review.py` | `c795ba0424735e2c` |
| `family_structure/scripts/full_phenotype_family_audit_20261009.py` | `scripts/02_lof_and_family/full_phenotype_family_audit_20261009.py` | `ef54f781f1df36bc` |
| `filtering/family_batches11_14/batch11_style_review/batch1114_clinvar_ndd_gene_review_20261027.py` | `scripts/07_gene_analysis_historical/batch1114_clinvar_ndd_gene_review_20261027.py` | `5b67c16d2927cdfa` |
| `filtering/family_batches11_14/scripts/01_site_qc_and_annovar.sh` | `scripts/01_site_qc/01_site_qc_and_annovar.sh` | `cd37774d9824e70c` |
| `filtering/family_batches11_14/scripts/02_lof_family_screen.py` | `scripts/02_lof_and_family/02_lof_family_screen.py` | `83778034671f496e` |
| `filtering/family_batches11_14/trio_validation_tools_20261110/validate_five_trios_from_vcf.py` | `scripts/06_trio_validation/validate_five_trios_from_vcf.py` | `e811699e1f39d26e` |

Portability edits: path defaults, nonessential names, generated email text, and private trio sample identifiers. The gnomAD v4.1 candidate decision function was not intentionally modified. **End-to-end equivalence has not been confirmed on the server.**

The original `family_analysis_td_updated_20261009.py` remains lab-only because it includes individual family identifiers and case notes; its logic is not replaced by this repository. The ClinVar shortlist script came from the separately retained analysis package and needs a confirmed server run/location.
