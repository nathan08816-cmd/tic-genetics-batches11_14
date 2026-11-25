# TIC Genetics: Batches 11–14 Variant Analysis

Rare loss-of-function (LoF) and coding-splicing variant analysis of joint-called exomes (GRCh38).

## Workflow

```text
Joint-called VCF + pedigrees
  → Site QC (bcftools) and functional annotation (ANNOVAR)
  → LoF/splicing selection, genotype QC, internal AF, family evidence
  → gnomAD v4.1 (echtvar), ClinVar and genomic-region flags (vcfanno)
  → Frequency/quality screening: PASS, REVIEW, EXCLUDE
  → ClinVar Pathogenic or Likely pathogenic variant and gene lists (PASS and REVIEW separate)
```

| Step | Script folder | Main output |
|---|---|---|
| 1. Site QC and ANNOVAR | `scripts/01_site_qc/` | Split-ALT annotated VCF |
| 2. LoF and family | `scripts/02_lof_and_family/` | Variant and family-event tables |
| 3. Local annotations | `scripts/03_local_annotation/` | gnomAD v4.1 and ClinVar evidence |
| 4. Candidate filtering | `scripts/04_local_v41_filter/` | PASS/REVIEW/EXCLUDE tables |
| 5. ClinVar shortlist | `scripts/05_clinvar_priority/` | Pathogenic/Likely pathogenic variants and genes |
| Additional | `scripts/06_trio_validation/`, `scripts/07_gene_analysis_historical/` | Trio and exploratory gene reviews |
| Diagnostics | `scripts/08_diagnostics/` | Aggregate frequency/ClinVar QC checks |

## Historical audit generator

`scripts/03_local_annotation/historical_audit/one_run_jointcall_finalizer.py` creates the comprehensive intermediate audit used by the later gnomAD v4.1 finalizer. This historical script uses earlier API/ANNOVAR frequency checks; **its candidate classifications are superseded by Stage 04**. Source genetic data and detailed candidate tables remain on the lab server.

## Main parameters

- Site `FILTER=PASS`; SNP VQSLOD ≥ −3.4427, indel VQSLOD ≥ 2.4465; capture interval ±20 bp.
- Frameshift insertion/deletion, stopgain, stoploss, or coding-splicing; exclude ncRNA-splicing-only.
- Genotype DP ≥ 10, GQ ≥ 20, ALT reads ≥ 4; heterozygous ALT fraction 0.25–0.75.
- Internal cohort AF ≤ 0.015; maximum observed gnomAD v4.1 AF < 0.001, gnomAD `FILTER=PASS`, AN > 0.
- Missing gnomAD observations, low-complexity/segmental-duplication regions, and non-autosomal cases go to `REVIEW`. `AC0` is not treated as verified rarity. Joint quality failures remain excluded under the original policy, with source-specific `AC0` cases reported separately for reassessment.
- ClinVar shortlist includes **exact** `Pathogenic` or `Likely_pathogenic` labels; combined classifications are separate.

## Results

| Stage | Result |
|---|---:|
| Individuals / families | 1,245 / 659 |
| Starting LoF/splicing alleles | 14,124 |
| PASS / REVIEW / EXCLUDE | 5,067 / 1,324 / 7,733 |
| PASS genes | 3,829 |
| ClinVar P/LP PASS shortlist | **331 variants / 304 genes** |
| Additional ClinVar P/LP REVIEW | **39 variants / 38 genes**, not PASS |

See [workflow commands](docs/WORKFLOW.md), [server paths](docs/PATHS.md), and [script provenance](docs/SOURCE_MAP.md).

```bash
cp config/paths.env.example config/paths.env
# Fill lab-private paths, then:
bash workflow.sh --check
bash workflow.sh --clinvar
```

The ClinVar step preserves the original PASS results and additionally writes a separate REVIEW table.
The historical counts above have **not** been recalculated after reference-based indel normalization. New QC runs require the exact GRCh38 FASTA (`REFERENCE_FASTA`) and downstream reannotation.

**Reproducibility limits:** The upstream generator of `ALL_14124_COMPLETE_AUDIT.tsv`, the lab-specific vcfanno configuration/gnomAD archive, and sample QC/pedigree generation are not in this repository. An end-to-end run still needs those lab-only dependencies. This is filtering and annotation, not a statistical association test or causal interpretation. ClinVar P/LP may refer to diseases unrelated to Tourette. No participant-level data is included. Obtain lab approval before sharing the repository outside authorized collaborators.
