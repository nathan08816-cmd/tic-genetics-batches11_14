#!/usr/bin/env python3
"""Read-only targeted VCF review for lab-supplied apparent de novo LoF/splicing sites.

Requires bcftools and the joint-call split-ALT QC-site VCF. The script does not
edit the source VCF or claim independent read/BAM validation.  It writes a
per-sample genotype/DP/GQ/AD report and one summary row per site. Targets supplied separately in a private TSV.
"""
import argparse
import csv
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_BCFTOOLS = 'bcftools'


def load_targets(path):
    """Read lab-private target list. Never commit this list to GitHub."""
    required = {'Variant_ID', 'Gene', 'Child_Sample', 'Father_Sample', 'Mother_Sample'}
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError('Expected TSV columns: ' + ', '.join(sorted(required)))
        rows = list(reader)
    if not rows:
        raise ValueError('No targets provided')
    for row in rows:
        parts = row['Variant_ID'].split('-', 3)
        if len(parts) != 4 or not parts[1].isdigit():
            raise ValueError('Invalid Variant_ID: ' + row['Variant_ID'])
        if not all(row.get(key) for key in required):
            raise ValueError('Missing required target fields')
    return rows


def tsv_write(path, fields, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        wr.writeheader()
        wr.writerows(rows)


def qc_sample_flags(path):
    if not path or not Path(path).is_file():
        return {}, 'Sample QC file unavailable; check separately'
    with open(path, encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    return {r['sample']:r for r in rows}, ''


def get_record(bcftools, vcf, chrom, pos, ref, alt, samples):
    errors = []
    for contig in (chrom, 'chr'+chrom):
        command = [bcftools, 'view', '-r', f'{contig}:{pos}-{pos}',
                   '-s', ','.join(samples), '-Ov', str(vcf)]
        proc = subprocess.run(command, text=True, capture_output=True, check=False)
        if proc.returncode:
            errors.append(f'{contig}: bcftools exit {proc.returncode}: {proc.stderr.strip()[:250]}')
            continue
        headers = [line.split('\t') for line in proc.stdout.splitlines() if line.startswith('#CHROM\t')]
        records = [line.split('\t') for line in proc.stdout.splitlines() if line and not line.startswith('#')]
        matches = [r for r in records if len(r) >= 10 and str(r[1]) == pos and
                   r[3].upper() == ref.upper() and r[4].upper() == alt.upper()]
        if not matches:
            errors.append(f'{contig}: no exact position/REF/ALT match ({len(records)} region records)')
            continue
        if len(matches) > 1:
            return None, None, f'Duplicate exact VCF site ({len(matches)})'
        if not headers:
            return None, None, 'No #CHROM header in VCF output'
        return headers[-1][9:], matches[0], ''
    return None, None, ' | '.join(errors)


def parse_sample(format_fields, raw):
    d = dict(zip(format_fields, raw.split(':')))
    gt = d.get('GT', '.')
    result = {'GT':gt,'DP':d.get('DP','.'),'GQ':d.get('GQ','.'),
              'AD':d.get('AD','.'),'PL':d.get('PL','.'),'AB':'.',
              'QC_call':'NO_CALL','QC_note':''}
    alleles = re.split(r'[|/]', gt)
    if not alleles or any(a not in ('0','1') for a in alleles):
        result['QC_note'] = 'Missing or nonsplit/multiallelic GT'
        return result
    try:
        dp, gq = int(result['DP']), float(result['GQ'])
    except (TypeError, ValueError):
        result['QC_note'] = 'Missing DP/GQ'
        return result
    if dp < 10 or gq < 20:
        result['QC_note'] = 'DP<10 or GQ<20'
        return result
    try:
        ads = [int(v) for v in result['AD'].split(',')]
    except ValueError:
        result['QC_note'] = 'Missing AD'
        return result
    if len(ads) < 2 or min(ads[:2]) < 0 or sum(ads[:2]) == 0:
        result['QC_note'] = 'Incomplete AD'
        return result
    ab = ads[1] / sum(ads[:2])
    result['AB'] = round(ab, 5)
    nalt = alleles.count('1')
    if nalt == 0:
        result['QC_call'] = 'REF' if ab <= 0.05 else 'NO_CALL'
        result['QC_note'] = 'Confident hom-ref' if ab <= 0.05 else 'Reference ALT-read fraction>0.05'
    elif nalt == len(alleles):
        good = ab >= 0.8 and ads[1] >= 4
        result['QC_call'] = 'ALT' if good else 'NO_CALL'
        result['QC_note'] = ('Haploid ALT: verify sex-specific ploidy' if good and len(alleles)==1
                             else 'Confident hom-alt' if good
                             else 'Hom-alt read balance/ALT reads failed')
    elif nalt == 1:
        good = ads[1] >= 4 and 0.25 <= ab <= 0.75
        result['QC_call'] = 'ALT' if good else 'NO_CALL'
        result['QC_note'] = 'Confident heterozygote' if good else 'Het read balance/ALT reads failed'
    else:
        result['QC_note'] = 'Unexpected genotype'
    return result


def parent_reference_warnings(role, call):
    """QC-REF is not enough to call a site de novo if parental ALT reads exist."""
    if role == 'affected_child' or call['QC_call'] != 'REF':
        return []
    try:
        alt_reads = int(str(call['AD']).split(',')[1])
    except (IndexError, TypeError, ValueError):
        return [f'{role}: reference AD unavailable; check original VCF']
    if alt_reads > 0:
        return [f'{role}: {alt_reads} ALT-supporting read(s) despite REF genotype; review potential mosaicism/artifact']
    return []


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--bcftools', default=DEFAULT_BCFTOOLS)
    ap.add_argument('--targets-tsv', required=True, help='Restricted table; never commit')
    ap.add_argument('--sample-qc', default='')
    ap.add_argument('--outdir', required=True)
    args=ap.parse_args()
    vcf=Path(args.vcf)
    bcftools = args.bcftools if Path(args.bcftools).is_file() else shutil.which(args.bcftools)
    if not vcf.is_file(): sys.exit('ERROR: VCF not found: '+str(vcf))
    if not bcftools: sys.exit('ERROR: bcftools not found: '+args.bcftools)
    dest=Path(args.outdir)
    if dest.exists() and any(dest.iterdir()):
        sys.exit('ERROR: output directory is nonempty; use a new --outdir to preserve prior analyses')
    dest.mkdir(parents=True, exist_ok=True)
    qc, qc_warning=qc_sample_flags(args.sample_qc)
    targets=load_targets(args.targets_tsv)
    per_sample=[]; summary=[]
    for target in targets:
        vid=target['Variant_ID']
        gene=target['Gene']
        fam=target.get('Family_Label','')
        ch,pos,ref,alt=vid.split('-',3)
        people=[target['Child_Sample'],target['Father_Sample'],target['Mother_Sample']]
        roles=['affected_child','father_unaffected','mother_unaffected']
        names, record, error = get_record(bcftools, vcf, ch, pos, ref, alt, people)
        out={'Variant_ID':vid,'Gene':gene,'Family_ID':fam,'Site_FILTER':'.',
             'Site_QUAL':'.','VQSLOD':'.','QD':'.','MQ':'.',
             'Child_GT':'.','Father_GT':'.','Mother_GT':'.',
             'Child_AD':'.','Father_AD':'.','Mother_AD':'.',
             'Child_DP':'.','Father_DP':'.','Mother_DP':'.',
             'Child_GQ':'.','Father_GQ':'.','Mother_GQ':'.',
             'Child_AB':'.','Father_AB':'.','Mother_AB':'.',
             'Genotype_QC_result':'UNVERIFIED','Issue':''}
        if error:
            out['Issue']=error
            summary.append(out)
            print(vid, 'UNVERIFIED', error, flush=True)
            continue
        site_info = dict((it.split('=',1)+[''])[:2] for it in record[7].split(';') if '=' in it)
        out.update({'Site_FILTER':record[6],'Site_QUAL':record[5],
                    'VQSLOD':site_info.get('VQSLOD','.'),'QD':site_info.get('QD','.'),
                    'MQ':site_info.get('MQ','.')})
        fmt=record[8].split(':')
        named=dict(zip(names,record[9:]))
        statuses=[]; warnings=[]
        for person,role in zip(people,roles):
            if person not in named:
                warnings.append(f'{role}: sample missing in output')
                statuses.append('MISSING')
                continue
            call=parse_sample(fmt,named[person]); statuses.append(call['QC_call'])
            qrow=qc.get(person,{})
            flags=qrow.get('qc_flags','') or ''
            if flags and flags.lower() not in ('nan','.', 'none'):
                warnings.append(f'{role}: sample QC flag={flags}')
            per_sample.append({'Variant_ID':vid,'Gene':gene,'Family_ID':fam,'Role':role,
                               'Sample_ID':person,**call,'Sample_QC_flags':flags or '.',
                               'Genetic_sex':qrow.get('genetic_sex','.')})
            prefix={'affected_child':'Child','father_unaffected':'Father','mother_unaffected':'Mother'}[role]
            for field in ['GT','AD','DP','GQ','AB']:
                out[f'{prefix}_{field}']=call[field]
            if call['QC_call']=='NO_CALL': warnings.append(f'{role}: {call["QC_note"]}')
            # A QC-REF parent may still have ALT reads under the 5% threshold.
            # Treat those as manual-review triggers for de novo analysis.
            if len(re.split(r'[|/]', str(call['GT']))) != 2 and call['QC_call'] != 'NO_CALL':
                warnings.append(f'{role}: haploid/non-diploid genotype; verify ploidy and sex chromosome')
            warnings.extend(parent_reference_warnings(role, call))
        if record[6]!='PASS':warnings.append('VCF FILTER is not PASS')
        if statuses!=['ALT','REF','REF']:
            warnings.append(f'Not ALT/REF/REF after genotype QC: {statuses}')
        if statuses==['ALT','REF','REF'] and record[6]=='PASS' and not warnings:
            out['Genotype_QC_result']='PATTERN_SUPPORTED_IN_SOURCE_VCF'
        else:
            out['Genotype_QC_result']='NEEDS_INVESTIGATION'
        out['Issue']='; '.join(warnings) if warnings else 'No genotype-QC warning; BAM/IGV review still required'
        summary.append(out)
        print(vid, out['Genotype_QC_result'], out['Issue'], flush=True)
    cols=['Variant_ID','Gene','Family_ID','Site_FILTER','Site_QUAL','VQSLOD','QD','MQ',
          'Child_GT','Father_GT','Mother_GT','Child_AD','Father_AD','Mother_AD',
          'Child_DP','Father_DP','Mother_DP','Child_GQ','Father_GQ','Mother_GQ',
          'Child_AB','Father_AB','Mother_AB','Genotype_QC_result','Issue']
    tsv_write(dest/'TRIO_TARGETED_VARIANT_VCF_CHECK.tsv',cols,summary)
    tsv_write(dest/'TRIO_TARGETED_SAMPLE_GENOTYPES.tsv',
              ['Variant_ID','Gene','Family_ID','Role','Sample_ID','GT','DP','GQ','AD','PL','AB','QC_call','QC_note','Sample_QC_flags','Genetic_sex'],per_sample)
    with open(dest/'TRIO_VCF_CHECK_SUMMARY.txt','w',encoding='utf-8') as w:
        w.write('Joint-call targeted trio genotype QC verification\n')
        w.write('Run UTC: '+datetime.datetime.now(datetime.timezone.utc).isoformat()+'\n')
        w.write('VCF: '+str(vcf)+'\n')
        w.write('BCFtools: '+str(bcftools)+'\n')
        w.write('Sample QC note: '+(qc_warning or 'Loaded sample QC file')+'\n')
        w.write('PASS genotype patterns: '+str(sum(r['Genotype_QC_result']=='PATTERN_SUPPORTED_IN_SOURCE_VCF' for r in summary))+'/'+str(len(targets))+'\n')
        for r in summary:
            w.write(f"{r['Variant_ID']} ({r['Gene']}): {r['Genotype_QC_result']} | {r['Issue']}\n")
        w.write('\nCAVEAT: A compatible genotype pattern and quality metrics DO NOT prove a biological de novo event.\n')
        w.write('Inspect parental/child aligned reads in BAM/CRAM with IGV or samtools: strand, read positions,\n')
        w.write('map/base quality, local indels/repeats, contamination, mosaic ALT reads, and sample identity.\n')
        w.write('ANNOVAR splicing/stopgain prediction also requires transcript-level interpretation.\n')
    with open(dest/'TRIO_RUN_METADATA.json','w',encoding='utf-8') as j:
        json.dump({'input_vcf':str(vcf),'output_directory':str(dest),
                   'results':{r['Variant_ID']:r['Genotype_QC_result'] for r in summary},
                   'read_level_manual_validation_completed':False},j,indent=2)
    print('OUTPUT:',dest)
    if any(r['Genotype_QC_result']!='PATTERN_SUPPORTED_IN_SOURCE_VCF' for r in summary):
        return 2
    return 0

if __name__ == '__main__':
    sys.exit(main())
