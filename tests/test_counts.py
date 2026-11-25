import csv
from pathlib import Path

p=Path(__file__).resolve().parents[1]/'results/AGGREGATE_COUNTS.tsv'
with p.open() as stream:
    rows={x['Metric']:int(x['Count']) for x in csv.DictReader(stream, delimiter='\t')}
assert rows['Local_v41_PASS_variants']+rows['Local_v41_REVIEW_variants']+rows['Local_v41_EXCLUDE_variants']==rows['LoF_or_coding_splicing_variants']
assert rows['ClinVar_exact_Pathogenic_PASS_variants']+rows['ClinVar_exact_Likely_pathogenic_PASS_variants']==rows['ClinVar_exact_P_or_LP_PASS_variants']
assert rows['Affected_individuals']+rows['Unaffected_individuals']+rows['Unknown_phenotype_individuals']==rows['Individuals']
print('Count checks passed')
