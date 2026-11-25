#!/usr/bin/env python3
"""Aggregate diagnostics from a private finished audit; never prints variant/sample IDs."""
import argparse
import csv
from collections import Counter
from pathlib import Path


def classify_function(name):
    text=(name or '').lower()
    if 'frameshift' in text: return 'frameshift'
    if 'stopgain' in text: return 'stopgain'
    if 'stoploss' in text: return 'stoploss'
    if 'splic' in text: return 'splicing'
    return 'other'


def diagnostics(rows):
    count=Counter(); kinds={}
    for row in rows:
        count['total']+=1
        status=row['Filter_Status']
        count['status_'+status]+=1
        if row.get('Local_gnomAD_Filter')=='AC0':count['joint_AC0']+=1
        ex=row.get('Local_gnomAD_exomes_filters','')
        gn=row.get('Local_gnomAD_genomes_filters','')
        if 'AC0' in ex or 'AC0' in gn:
            count['source_AC0']+=1
            if int(row.get('Affected_Carriers') or 0)>0:
                count['source_AC0_affected']+=1
                if row.get('ClinVar','').split(' [',1)[0] in ('Pathogenic','Likely_pathogenic'):
                    count['source_AC0_affected_P_or_LP']+=1
        label=row.get('ClinVar','').split(' [',1)[0]
        if label in ('Pathogenic','Likely_pathogenic'):
            count[f'clinvar_plp_{status}']+=1
        name=classify_function(row.get('Function',''))
        k=kinds.setdefault(name,Counter())
        k['total']+=1
        try: an=int(row.get('Local_gnomAD_AN') or '-1')
        except ValueError:an=-1
        if an<=0:k['no_usable_gnomad_AN']+=1
    return count,kinds


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--audit',required=True,type=Path)
    args=ap.parse_args()
    with args.audit.open(newline='',encoding='utf-8') as f:
        rd=csv.DictReader(f,delimiter='\t')
        required={'Variant_ID','Filter_Status','ClinVar','Function','Local_gnomAD_AN',
                  'Local_gnomAD_exomes_filters','Local_gnomAD_genomes_filters'}
        if not required.issubset(rd.fieldnames or []):
            ap.error('Input missing final-audit fields')
        count,kinds=diagnostics(rd)
    print('Aggregate audit summary (no variant or sample identifiers)')
    for key in sorted(count):print(f'{key}\t{count[key]}')
    print('Functional category\tTotal\tNo usable gnomAD AN\tRate (%)')
    for kind,k in sorted(kinds.items()):
        print(f'{kind}\t{k["total"]}\t{k["no_usable_gnomad_AN"]}\t{k["no_usable_gnomad_AN"]/k["total"]*100:.2f}')
    print('Interpretation: Missing gnomAD data is not evidence of rarity; a higher indel missingness rate is not proof of an indel representation mismatch.')

if __name__=='__main__':main()
