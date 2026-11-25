#!/usr/bin/env python3
"""Exact-allele local gnomAD v4.1, ClinVar, and region annotation extractor.
Uses INFO keys from the local reference configuration. No missing-record=0 assumptions.
"""
import argparse,csv,gzip,re,collections,math
from pathlib import Path

def opening(p):return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p,encoding='utf-8')
def value(info,key):return info.get(key,'').strip()
def positivefloat(s):
    if s in ('','.','-1','NA','MISSING','None','null'):return None
    try:
        x=float(s)
        if not math.isfinite(x) or not 0<=x<=1:raise ValueError
        return x
    except ValueError:raise RuntimeError(f'Unexpected gnomAD AF {s!r}')
def flag(info,k):return '1' if k in info and info[k].upper() not in ('0','FALSE','.') else '0'
def summary(r):
    ans=[r.get('ClinVar_Significance',''),r.get('ClinVar_ID','')]
    return any(x and x!='.' for x in ans)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf',required=True,type=Path);p.add_argument('--output',type=Path,default=Path('local_gnomad_flag_evidence.tsv'))
    p.add_argument('--expected-variants',type=int,default=14124)
    p.add_argument('--show-info-ids',action='store_true')
    p.add_argument('--filter-key',action='append',default=[],help='Must be gnomad_filter per the lab document')
    a=p.parse_args()
    if a.filter_key and a.filter_key != ['gnomad_filter']:
        raise RuntimeError("Document specifies sole joint status field 'gnomad_filter'; do not treat dataset-specific reason strings as separate PASS requirements")
    headers=set();rows=[];ids=set();counts=collections.Counter()
    for line in opening(a.vcf):
        if line.startswith('##INFO=<ID='):
            m=re.match(r'##INFO=<ID=([^,>]+)',line)
            if m:headers.add(m.group(1))
        if line.startswith('#'):continue
        t=line.rstrip('\n').split('\t')
        if len(t)<8:raise RuntimeError('Invalid VCF data row')
        chrom,pos,ref,alt=t[0],t[1],t[3],t[4]
        if ',' in alt:raise RuntimeError(f'Unsplit multiallelic VCF at {chrom}:{pos}')
        key=f'{chrom.upper().removeprefix("CHR")}-{int(pos)}-{ref.upper()}-{alt.upper()}'
        if key in ids:raise RuntimeError('Duplicate normalized key '+key)
        ids.add(key)
        info={}
        for chunk in t[7].split(';'):
            if '=' in chunk:
                k,v=chunk.split('=',1);info[k]=v
            elif chunk and chunk!='.':info[chunk]='1'
        f=value(info,'gnomad_filter')
        try:an=int(value(info,'gnomad_AN') or '-1')
        except ValueError:raise RuntimeError('Invalid gnomad_AN at '+key)
        keys=['gnomad_AF','gnomad_AF_grpmax','gnomad_AF_exomes','gnomad_AF_genomes']
        freqs={k:positivefloat(value(info,k)) for k in keys}
        valid=[x for x in freqs.values() if x is not None]
        af=max(valid) if valid else None
        if f=='AC0':
            status='UNRESOLVED';detail='AC0: no QC-passing ALT evidence; review'
        elif an<=0 or af is None:
            status='UNRESOLVED';detail='NO_GNOMAD_RECORD' if an<0 else 'NO_VALID_FREQUENCY'
        elif f=='PASS':
            status='EXCLUDE_AF' if af>=.001 else 'PASS';detail='gnomad_filter=PASS'
        elif f in ('','.','MISSING','UNKNOWN'):
            status='UNRESOLVED';detail='MISSING_GNOMAD_FILTER'
        else:
            status='FLAGGED';detail='gnomad_filter='+f
        counts[status]+=1
        row={'Variant_ID':key,'Flag_Assessment':status,'Flag_Detail':detail,'Flag_Source':'echtvar_gnomad_v4.1+vcfanno_exact_allele',
             'Local_gnomAD_AF':'' if af is None else format(af,'.9g'),
             'Local_gnomAD_AN':an,'Local_gnomAD_Filter':f,
             'Local_gnomAD_AF_joint':value(info,'gnomad_AF'),'Local_gnomAD_AF_grpmax':value(info,'gnomad_AF_grpmax'),
             'Local_gnomAD_AF_exomes':value(info,'gnomad_AF_exomes'),'Local_gnomAD_AF_genomes':value(info,'gnomad_AF_genomes'),
             'Local_gnomAD_fafmax95':value(info,'gnomad_fafmax95'),
             'Local_gnomAD_exomes_filters':value(info,'gnomad_exomes_filters'),
             'Local_gnomAD_genomes_filters':value(info,'gnomad_genomes_filters'),
             'Local_LCR':flag(info,'gnomad_lcr'),'Local_SegDup':flag(info,'gnomad_segdup'),
             'ClinVar_Significance':value(info,'clinvar_sig'), 'ClinVar_ID':value(info,'clinvar_id'),
             'ClinVar_Stars':value(info,'clinvar_stars'),'ClinVar_Review':value(info,'clinvar_revstat'),
             'ClinVar_Disease':value(info,'clinvar_disease')}
        rows.append(row)
    if a.show_info_ids:
        print('INFO header names:\n'+'\n'.join(sorted(headers)))
        print('Observed exact-allele variants:',len(rows))
        return
    needed={'gnomad_filter','gnomad_AN','gnomad_AF_grpmax','clinvar_sig','clinvar_id','gnomad_lcr','gnomad_segdup'}
    if not needed<=headers:raise RuntimeError('Missing expected gnomAD/ClinVar fields: '+str(sorted(needed-headers))+'; inspect vcfanno base path / lua / echtvar output')
    if len(rows)!=a.expected_variants:raise RuntimeError(f'Expected {a.expected_variants}, saw {len(rows)}')
    if not any(int(r['Local_gnomAD_AN'])>0 for r in rows):raise RuntimeError('No local gnomAD records matched! Check chromosome prefix and GRCh38 representation.')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with open(a.output,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    print('Wrote',a.output,'rows',len(rows),'outcomes',dict(counts),'matched gnomAD alleles',sum(int(r['Local_gnomAD_AN'])>0 for r in rows),'ClinVar matched',sum(summary(r) for r in rows))
    print('Unresolved or absent gnomAD record does NOT mean frequency zero. Source dataset is v4.1 and API results may be v4.0.')

if __name__=='__main__':main()
