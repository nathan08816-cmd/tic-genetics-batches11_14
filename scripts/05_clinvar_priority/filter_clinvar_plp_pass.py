#!/usr/bin/env python3
"""Extract exact ClinVar Pathogenic/Likely_pathogenic from local-v4.1 PASS audit.

All external-frequency, QC and family filtering has already been done upstream.
No inference about Tourette causality. Does not change source files.
Python 3.9+, standard library only.
"""
import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

PRIMARY = {'Pathogenic', 'Likely_pathogenic'}
COMBINED = 'Pathogenic/Likely_pathogenic'
VARIANT_COLUMNS = [
    'Variant_ID', 'Gene', 'Function', 'ClinVar_Class', 'ClinVar',
    'ClinVar_Disease', 'ClinVar_Review_Status', 'ClinVar_Stars',
    'gnomAD_AF', 'Internal_AF', 'Affected_Carriers', 'Unaffected_Carriers',
    'Affected_Families', 'Strict_Families', 'Inconsistent_Families',
    'Segregation', 'Filter_Status'
]
GENE_COLUMNS = [
    'Gene', 'ClinVar_P_or_LP_Variants', 'Pathogenic_Variants',
    'Likely_pathogenic_Variants', 'Distinct_Affected_Families',
    'Affected_Variant_Carrier_Observations',
    'Unaffected_Variant_Carrier_Observations',
    'Strict_Compatible_Variants', 'Lowest_gnomAD_AF',
    'ClinVar_Diseases', 'Variant_IDs'
]


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        for row in reader:
            yield row


def split_genes(s):
    return sorted({g.strip() for g in s.split(';') if g.strip() not in ('', '.', 'NA')})


def label(row):
    # e.g. 'Pathogenic [4850530]' -> 'Pathogenic'.  Do NOT substring match.
    return (row['ClinVar'].split(' [', 1)[0]).strip()


def write_tsv(path, fields, rows):
    with open(path, 'w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True,
                        help='Completed ALL_14124_LOCAL_V41_REASSESSED_AUDIT.tsv')
    parser.add_argument('--events', type=Path, required=True,
                        help='Existing all_family_carrier_events.tsv.gz')
    parser.add_argument('--outdir', type=Path, required=True,
                        help='NEW empty output directory (never overwrites)')
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        parser.error('Output directory is nonempty; refusing to overwrite')
    if not args.audit.is_file() or not args.events.is_file():
        parser.error('Input audit/events file not found')

    # All selected variants originate from preexisting PASS calls.
    all_pass = []
    seen = set()
    for r in read_tsv(args.audit):
        if r['Filter_Status'] != 'PASS':
            continue
        vid = r['Variant_ID']
        if vid in seen:
            raise ValueError('Repeated PASS Variant_ID ' + vid)
        seen.add(vid)
        all_pass.append(r)
    if not all_pass:
        raise ValueError('No PASS rows in audit')
    primary = [r for r in all_pass if label(r) in PRIMARY]
    combined = [r for r in all_pass if label(r) == COMBINED]
    if not primary:
        raise ValueError('No exact ClinVar Pathogenic or Likely_pathogenic rows')

    selected_ids = {r['Variant_ID'] for r in primary}
    fam_by_vid = defaultdict(set)
    affected_by_vid = Counter()
    for e in read_tsv(args.events):
        vid = e['Variant_ID']
        if vid not in selected_ids:
            continue
        n = int(e['Affected_ALT'])
        if n > 0:
            if not e['Family_ID']:
                raise ValueError('Missing Family_ID for '+vid)
            fam_by_vid[vid].add(e['Family_ID'])
            affected_by_vid[vid] += n
    for r in primary:
        vid = r['Variant_ID']
        if affected_by_vid[vid] != int(r['Affected_Carriers']):
            raise ValueError('Affected-carrier mismatch for '+vid)
        if len(fam_by_vid[vid]) != int(r['Affected_Families']):
            raise ValueError('Affected-family mismatch for '+vid)
        if r.get('ClinVar_Source') != 'exact_allele_vcfanno':
            raise ValueError('Non-exact ClinVar evidence for '+vid)

    def var_record(r):
        return {**r, 'ClinVar_Class': label(r)}

    primary_out = [var_record(r) for r in primary]
    combined_out = [var_record(r) for r in combined]
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir/'ClinVar_P_or_LP_PASS_Variants.tsv', VARIANT_COLUMNS, primary_out)
    write_tsv(args.outdir/'ClinVar_Combined_PLP_PASS_Additional.tsv', VARIANT_COLUMNS, combined_out)

    # Do not mislabel sums of per-allele carrier counts as unique persons.
    genes = defaultdict(lambda: {'ids': set(), 'p': 0, 'lp': 0,
                                 'affected_families': set(), 'affected_obs': 0,
                                 'unaffected_obs': 0, 'strict': 0, 'freqs': [],
                                 'diseases': set()})
    for r in primary:
        for gene in split_genes(r['Gene']):
            x = genes[gene]
            if r['Variant_ID'] in x['ids']:
                continue
            x['ids'].add(r['Variant_ID'])
            x['p'] += label(r) == 'Pathogenic'
            x['lp'] += label(r) == 'Likely_pathogenic'
            x['affected_families'].update(fam_by_vid[r['Variant_ID']])
            x['affected_obs'] += int(r['Affected_Carriers'])
            x['unaffected_obs'] += int(r['Unaffected_Carriers'])
            x['strict'] += r['Segregation'] == 'Strict-compatible'
            if r['gnomAD_AF']:
                x['freqs'].append(float(r['gnomAD_AF']))
            disease = r.get('ClinVar_Disease', '')
            if disease and disease.lower() != 'not_provided':
                x['diseases'].add(disease)

    gene_out = []
    for gene, x in genes.items():
        gene_out.append({
            'Gene': gene,
            'ClinVar_P_or_LP_Variants': len(x['ids']),
            'Pathogenic_Variants': x['p'],
            'Likely_pathogenic_Variants': x['lp'],
            'Distinct_Affected_Families': len(x['affected_families']),
            'Affected_Variant_Carrier_Observations': x['affected_obs'],
            'Unaffected_Variant_Carrier_Observations': x['unaffected_obs'],
            'Strict_Compatible_Variants': x['strict'],
            'Lowest_gnomAD_AF': min(x['freqs']) if x['freqs'] else '',
            'ClinVar_Diseases': '; '.join(sorted(x['diseases'])),
            'Variant_IDs': '; '.join(sorted(x['ids']))
        })
    gene_out.sort(key=lambda r: (-r['Distinct_Affected_Families'],
                                -r['ClinVar_P_or_LP_Variants'],
                                -r['Strict_Compatible_Variants'], r['Gene']))
    write_tsv(args.outdir/'ClinVar_P_or_LP_PASS_Genes.tsv', GENE_COLUMNS, gene_out)
    combined_genes = {g for r in combined for g in split_genes(r['Gene'])}
    all_genes = set(genes) | combined_genes
    note = f'''TIC Genetics Batches 11-14 | ClinVar pathogenic-priority subset

Input: local gnomAD v4.1 final PASS audit; original family carrier events.
PASS variants in source: {len(all_pass)}
EXACT ClinVar Pathogenic: {sum(label(r)=='Pathogenic' for r in primary)}
EXACT ClinVar Likely_pathogenic: {sum(label(r)=='Likely_pathogenic' for r in primary)}
PRIMARY TOTAL variants: {len(primary)}
PRIMARY TOTAL distinct annotated genes: {len(genes)}
Distinct affected families/allele records validated against family-event file: {len(primary)}/{len(primary)}
ADDITIONAL ClinVar Pathogenic/Likely_pathogenic combined label: {len(combined)} variants in {len(combined_genes)} genes
UNION including combined label: {len(primary)+len(combined)} variants in {len(all_genes)} genes

PRIMARY selections require exact ClinVar classification of 'Pathogenic' or 'Likely_pathogenic', an exact-allele ClinVar vcfanno match, and Filter_Status=PASS.
Excluded from primary: Pathogenic/Likely_pathogenic combined label (supplied separately), conflicting interpretations, VUS, benign, not annotated, RECORD/no-significance, Pathogenic|Affects, and all REVIEW/EXCLUDE variants.

Interpretation: ClinVar disease classifications may concern disorders OTHER THAN Tourette. PASS means passing the upstream computational filters, not gene-disease validation. Gene tables split semicolon-delimited RefGene symbols, dedupe repeated symbols, and count unique affected FAMILIES across selected variants. 'Affected_Variant_Carrier_Observations' sums variant-level carrier counts; it is NOT the number of unique people per gene (requires joint VCF genotypes to compute).

No BAM/CRAM read review, reannotation, or new ClinVar query was performed during this extraction. The labels reflect the saved local ClinVar snapshot.
'''
    (args.outdir/'CLINVAR_SELECTION_SUMMARY.txt').write_text(note, encoding='utf-8')
    print(note)
    print('OUTPUT DIRECTORY:', args.outdir)

if __name__ == '__main__':
    main()
