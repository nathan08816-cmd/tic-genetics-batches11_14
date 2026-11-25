# Server paths to confirm

Relative to `$PROJECT_ROOT` (`TS_ANNOVAR_analysis`):

| Stage | Known server location |
|---|---|
| QC and family screen scripts | `filtering/family_batches11_14/scripts/` |
| Intermediate LoF results | `filtering/family_batches11_14/lof_joint_20261026/results/` |
| Earlier full audit/local flags | `filtering/family_batches11_14/final_jointcall_20261110_171610/` |
| Local gnomAD annotation scripts | `JointCall_Configured_Local_gnomAD_20261110/` |
| Final gnomAD v4.1 code | `JointCall_LocalV41_Final_20261110/` |
| Final gnomAD v4.1 results | `filtering/family_batches11_14/local_v41_final_20261110_174914/` |
| Family classification | `family_structure/scripts/` and `family_structure/results/` |
| Earlier NDD/ClinVar review | `filtering/family_batches11_14/batch11_style_review/` |
| Trio review | `filtering/family_batches11_14/trio_five_read_qc_20261110/` |

**Still to verify:** upstream VCF/VQSR/phenotype paths, local annotation resource versions, specific commands from run logs, and the output path for the 331-variant ClinVar shortlist. The shortlist script is available in the supplied analysis package, but was not found in the scanned server inventory.

Never commit private input paths containing participant information, participant-level tables, VCFs, BAMs, CRAMs, or clinical/pedigree workbooks.
