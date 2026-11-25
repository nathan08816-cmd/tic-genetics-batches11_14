# Workflow and commands

Commands below describe the executed approach. Set all paths to their lab-local locations and run on new output directories; do not overwrite completed analyses.

| Stage | Inputs | Operation | Outputs |
|---|---|---|---|
| 1. Site QC | Joint VCF, VQSR scores, capture intervals | Site PASS, VQSLOD, split ALT | Site-QC VCF |
| 2. ANNOVAR | Site-QC VCF, RefGene | Gene and consequence annotation | Multianno VCF |
| 3. Family screen | Genotype VCF, phenotypes, sample QC | DP/GQ/allele balance; family segregation and internal AF | 14,124-allele audit and family events |
| 4. Local annotation | LoF-site VCF, gnomAD v4.1, ClinVar | echtvar and vcfanno; exact allele matching | Local evidence TSV |
| 5. Final filtering | QC audit, evidence, genotype VCF | AF, quality, region and carrier filters | 5,067 PASS / 1,324 REVIEW / 7,733 EXCLUDE |
| 6. ClinVar shortlist | Final audit and family events | Exact Pathogenic/ Likely_pathogenic labels | 331 variants / 304 genes |

The historical intermediate audit (`ALL_14124_COMPLETE_AUDIT.tsv`) was generated with
`scripts/03_local_annotation/historical_audit/one_run_jointcall_finalizer.py` using the merged variant table, ANNOVAR records, genotype VCF,
phenotype and sample QC files, family events, and gnomAD/local-flag evidence.
Its original filtering decisions are **not** the final v4.1 results.
The original source paths are recorded in `docs/PATHS.md`.

Use `python3 scripts/03_local_annotation/historical_audit/one_run_jointcall_finalizer.py --help` for input flags. The script requires
lab-local data and is not invoked by `workflow.sh`.

## Commands

```bash
# Set paths on the lab server; config/paths.env must not be committed.
export PROJECT_ROOT=/private/path/TS_ANNOVAR_analysis
export SOURCE_DATA_DIR=/private/path/upstream_joint_calls
export SAMPLE_QC=/private/path/sample_qc.tsv
export PHENOTYPE_XLSX=/private/path/phenotypes.xlsx
export ANNOVAR_DIR=/private/path/annovar
export HUMANDB=/private/path/humandb
export REFERENCE_FASTA=/private/path/exact_jointcall_GRCh38_reference.fa

# 1. Site QC and 2. RefGene annotation (write to new run directories only)
PROJECT="$PROJECT_ROOT" DATA="$SOURCE_DATA_DIR" bash scripts/01_site_qc/01_site_qc_and_annovar.sh qc
PROJECT="$PROJECT_ROOT" DATA="$SOURCE_DATA_DIR" bash scripts/01_site_qc/01_site_qc_and_annovar.sh annotate

# 3. LoF/coding-splicing and family screen
python3 scripts/02_lof_and_family/02_lof_family_screen.py \
  --vcf /private/path/splitALT.vcf.gz --annovar-vcf /private/path/multianno.vcf \
  --phenotype "$PHENOTYPE_XLSX" --data-dir "$SOURCE_DATA_DIR" --outdir /private/new/family_screen

# 4. Local annotations (echtvar + vcfanno)
bash scripts/03_local_annotation/annotate_local_gnomad_flags.sh \
  --vcf /private/path/lof_site.vcf.gz --outdir /private/new/local_annotations \
  --echtvar-archive /private/path/gnomad_v4.1.echtvar.zip \
  --vcfanno-config /private/path/vcfanno.toml
# Run gnomad_flag_evidence_from_vcf.py with the INFO mappings confirmed for the reference.

# 5. gnomAD v4.1 finalizer using saved source-audit and local annotations
SOURCE_RUN=/private/path/completed_prior_run \
  bash scripts/04_local_v41_filter/run_local_v41_final.sh

# 6. ClinVar shortlist; use filled, ignored config/paths.env
bash workflow.sh --check
bash workflow.sh --clinvar
```

## Checks before changing the original results

```bash
# Aggregated diagnostics only (no sample/variant rows are published):
python3 scripts/08_diagnostics/audit_annotation_summary.py --audit "$FINAL_V41_AUDIT"
# Unit and synthetic integration checks (no lab data required):
python3 -m unittest discover -s tests -v
```

Stage 01 now requires a reference FASTA and runs `bcftools norm -f` to
left-normalize indels. **This was not part of the historical output**: running
this updated step will change some allele keys and requires rebuilding downstream
annotations, family evidence and gnomAD/ClinVar exact-allele lookups. Never mix
pre- and post-normalization outputs.

The current ClinVar shortlist also writes a **separate** exact P/LP REVIEW
variant table. These records are NOT promoted to PASS. In the saved audit, 39
additional P/LP variants were REVIEW; the original 331 PASS variants and 304
genes remain unchanged.

Current gnomAD reference quality policy excludes joint-filter failures. The
source-specific `AC0` subset needs manual interpretation and is reported in
aggregate diagnostics; it must not be equated with a verified good AF=0 call.

## QC and interpretation

- Site: `FILTER=PASS`, capture ±20 bp, VQSLOD SNP ≥ −3.4427 / indel ≥ 2.4465.
- Genotype: DP ≥ 10, GQ ≥ 20; ALT reads ≥ 4; heterozygous ALT fraction 0.25–0.75; homozygous-ref ALT fraction ≤ 0.05; homozygous-alt ≥ 0.80. Failed calls become missing.
- Internal AF: ALT alleles / confidently called cohort alleles; exclude > 0.015.
- gnomAD: maximum available overall/exome/genome/ancestry AF < 0.001, allele number > 0, quality PASS. Unobserved variants remain REVIEW.
- LCR, segmental duplication and non-autosomal cases remain REVIEW. ClinVar is annotated independently of gnomAD.
- Segregation is supporting evidence, not a prerequisite for PASS. De novo-compatible genotype patterns require read-level confirmation.
- The separate ClinVar shortlist includes *exact* `Pathogenic` or `Likely_pathogenic`, regardless of the ClinVar disease association. ClinVar review stars and inheritance are not selection criteria.

## Other work

- Pedigree classification, family recurrence, and NDD gene-list overlap: historical scripts in the server inventory.
- Trio validation: pass a lab-only target TSV containing sample identifiers; never commit that table.
- Exploratory ToppGene/HumanBase analysis was performed for a separate expression dataset, not this ClinVar shortlist.

The script source was recovered from the server, but these portable copies have not been fully rerun end-to-end after removing environment-specific paths and private identifiers. The original full-audit generator, source-reference configuration, and unrelated/sample QC generation remain lab-private inputs; they must be documented and tested for end-to-end replication. Internal AF uses the full cohort (0.015), not the unrelated-only AF; this and any downstream enrichment/statistical testing require a separately agreed analysis plan.
