#!/usr/bin/env python3
"""Reproduce the Batch 11 LoF/splicing gene-review categories for Batch 11--14.

ANNOTATIONS, NOT DIAGNOSES. Reads site/genotype ANNOVAR rollup, 24 gnomAD
results, exact-key older B11 gnomAD results, a read-only ClinVar GRCh38 VCF,
and a curated NDD gene list. Writes NEW files only to --outdir.

--reference-fasta optional: when supplied, normalize candidate and nearby
ClinVar alleles with bcftools norm before an exact REF/ALT comparison. Without
it, matching is EXPLICITLY coordinate+REF+ALT exact and indels are unresolved
when absent. This does NOT establish any variant's disease pathogenicity.

Requires bcftools on the Rutgers server, Python 3 standard library.
"""
import argparse
import collections
import csv
import datetime
import gzip
import json
import os
import re
import subprocess
from pathlib import Path

PROJECT=Path(os.environ.get('PROJECT_ROOT', '/path/to/project'))
RESULTS=PROJECT/'filtering/family_batches11_14/lof_joint_20261026/results'
BASE=Path(__file__).resolve().parent
CONSEQ=['stopgain','frameshift_deletion','frameshift_insertion','stoploss','splicing']
MISSING={'','.', 'NA', 'None', 'nan'}
CLIN_FIELDS=['ClinVar_match_status','ClinVar_exact_Variant_ID','ClinVar_CLNSIG','ClinVar_category','ClinVar_CLNREVSTAT','ClinVar_CLNDN','ClinVar_record_ID','ClinVar_allele_ID','ClinVar_source']
GN_FIELDS=['gnomAD_screen_source','gnomAD_screen_status','gnomAD_screen_max_AF','gnomAD_screen_decision','gnomAD_screen_reason']
NDD_FIELDS=['NDD_overlap','NDD_source','NDD_disease']

def date(): return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
def log(s): print('['+date()+'] '+s,flush=True)
def op(p,mode='rt'):
    return gzip.open(p,mode,encoding='utf-8',newline='') if str(p).endswith('.gz') else open(p,mode,encoding='utf-8',newline='')
def read_tsv(p):
    if not Path(p).is_file(): raise FileNotFoundError('Missing '+str(p))
    with op(p) as f:
        d=csv.DictReader(f,delimiter='\t')
        return list(d),list(d.fieldnames or [])
def dump_tsv(path,rows,fields):
    with op(path,'wt') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def norm_chrom(chrom):
    c=str(chrom).strip()
    if c.lower().startswith('chr'): c=c[3:]
    if c=='M': c='MT'
    return c.upper()

def norm_key(chrom,pos,ref,alt):
    return '-'.join([norm_chrom(chrom),str(pos),str(ref).upper(),str(alt).upper()])

def is_standard(row):
    return (bool(re.fullmatch('[ACGT]+',str(row['Ref']).upper()))
            and bool(re.fullmatch('[ACGT]+',str(row['Alt']).upper()))
            and int(row['Pos'])>=1)

def genes(row):
    v=str(row.get('Gene_refGene','')).replace(',', ';')
    return list(dict.fromkeys(g.strip() for g in v.split(';') if g.strip() not in MISSING))

def run(cmd, output=None):
    log('Running: '+' '.join(str(x) for x in cmd))
    if output:
        with open(output,'w',encoding='utf8') as o:
            subprocess.run([str(x) for x in cmd],stdout=o,check=True)
    else: subprocess.run([str(x) for x in cmd],check=True)

def clinvar_class(text):
    s=str(text).strip()
    if s in MISSING: return 'NOT_REPORTED'
    lo=s.lower().replace(' ','_')
    if any(x in lo for x in ('conflict','_vs_','/uncertain_significance','uncertain_significance/','pathogenic/benign','benign/pathogenic')):
        return 'CONFLICTING_OR_MIXED'
    if ('pathogenic' in lo) and ('benign' in lo):return 'CONFLICTING_OR_MIXED'
    if ('pathogenic' in lo) and ('uncertain' in lo):return 'CONFLICTING_OR_MIXED'
    if 'low_penetrance' in lo:return 'REVIEW_LOW_PENETRANCE'
    if 'pathogenic' in lo:return 'PATHOGENIC_OR_LIKELY_PATHOGENIC'
    if 'benign' in lo:return 'BENIGN_OR_LIKELY_BENIGN'
    if 'uncertain' in lo:return 'UNCERTAIN_SIGNIFICANCE'
    return 'OTHER_CLINICAL_SIGNIFICANCE'

def parse_info(s):
    d={}
    for part in s.split(';'):
        if '=' in part:
            k,v=part.split('=',1);d[k]=v
        elif part:d[part]='YES'
    return d

def vcf_records(p):
    with open(p,'rt',encoding='utf-8',errors='replace') as f:
        for line in f:
            if not line or line[0]=='#':continue
            a=line.rstrip('\n').split('\t')
            if len(a)<8:continue
            for alt in a[4].split(','):
                yield {'key':norm_key(a[0],a[1],a[3],alt),
                       'chrom':a[0],'pos':a[1],'ref':a[3],'alt':alt,'id':a[2],
                       'info':parse_info(a[7])}

def clinvar_lookup(rows,clinvar,out,ref=None,flank=250):
    """Use a bounded indexed ClinVar query; with a FASTA normalize both sides."""
    # Read ClinVar chrom prefixes via header, without reading full 191MB archive.
    try:
        header=subprocess.check_output(['bcftools','view','-h',str(clinvar)],text=True)
    except FileNotFoundError:
        raise RuntimeError('bcftools is not on PATH; use the same environment as the earlier variant pipeline.')
    cv_chrom={}
    for line in header.splitlines():
        if line.startswith('##contig=<ID='):
            v=line.split('ID=',1)[1].split(',',1)[0].split('>',1)[0]
            cv_chrom.setdefault(norm_chrom(v),v)
    if not cv_chrom: raise RuntimeError('Could not read ClinVar VCF contigs.')
    # Look within 250bp of each input anchor. Do not interpret absence as ClinVar benign.
    bed=[]
    for r in rows:
        chrom=cv_chrom.get(norm_chrom(r['Chr']))
        if chrom is None:continue
        p=int(r['Pos']);bed.append((chrom,max(0,p-1-flank),p+flank))
    chrom_order={v:i for i,v in enumerate(cv_chrom.values())}
    bed.sort(key=lambda x:(chrom_order.get(x[0],999),x[1],x[2]))
    merged=[]
    for x in bed:
        if merged and merged[-1][0]==x[0] and x[1]<=merged[-1][2]:
            merged[-1]=(x[0],merged[-1][1],max(x[2],merged[-1][2]))
        else:merged.append(x)
    with open(out/'candidate_intervals_plus250bp.bed','w') as f:
        for x in merged:f.write(f'{x[0]}\t{x[1]}\t{x[2]}\n')
    sites=out/'clinvar_near_candidate_positions.vcf'
    run(['bcftools','view','-R',str(out/'candidate_intervals_plus250bp.bed'),'-Ov','-o',str(sites),str(clinvar)])
    clinvar_path=sites
    original_map={r['Variant_ID']:r['Variant_ID'] for r in rows if is_standard(r)}
    if ref:
        # Require an existing FASTA index; never generate sidecars next to an original FASTA.
        idx=Path(str(ref)+'.fai')
        if not idx.is_file():raise RuntimeError('Reference needs existing .fai index: '+str(idx))
        fai=[line.split('\t')[0] for line in idx.read_text().splitlines() if line]
        ref_chrom={norm_chrom(v):v for v in fai}
        # All candidates must be build-compatible; don't silently change mismatched REF bases.
        for r in rows:
            if is_standard(r) and norm_chrom(r['Chr']) not in ref_chrom:
                raise RuntimeError('Reference FASTA lacks chromosome '+str(r['Chr']))
        rename_file=out/'clinvar_to_reference_chroms.tsv'
        with open(rename_file,'w') as f:
            for canon,old in cv_chrom.items():
                if canon in ref_chrom and old!=ref_chrom[canon]:f.write(old+'\t'+ref_chrom[canon]+'\n')
        mapped=out/'clinvar_nearby_renamed.vcf'
        if rename_file.stat().st_size:
            run(['bcftools','annotate','--rename-chrs',str(rename_file),'-Ov','-o',str(mapped),str(sites)])
        else:mapped=sites
        newcv=out/'clinvar_nearby_normalized.vcf'
        run(['bcftools','norm','-f',str(ref),'-m','-any','-Ov','-o',str(newcv),str(mapped)])
        clinvar_path=newcv
        query_vcf=out/'candidates_for_normalization.vcf'
        with open(query_vcf,'w') as f:
            f.write('##fileformat=VCFv4.2\n')
            for canon,v in ref_chrom.items():
                if canon in cv_chrom:f.write(f'##contig=<ID={v}>\n')
            f.write('##INFO=<ID=OID,Number=1,Type=String,Description="Original variant ID">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            for r in rows:
                if not is_standard(r):continue
                orig=r['Variant_ID']
                f.write(f"{ref_chrom[norm_chrom(r['Chr'])]}\t{r['Pos']}\t.\t{r['Ref']}\t{r['Alt']}\t.\tPASS\tOID={orig}\n")
        norm_candidates=out/'candidates_normalized.vcf'
        run(['bcftools','norm','-f',str(ref),'-m','-any','-Ov','-o',str(norm_candidates),str(query_vcf)])
        original_map={r['info'].get('OID'):r['key'] for r in vcf_records(norm_candidates) if r['info'].get('OID')}
        if len(original_map)!=sum(is_standard(r) for r in rows):
            raise RuntimeError('Candidate normalization lost an allele, inspect '+str(norm_candidates))
    match_by_key=collections.defaultdict(list)
    for item in vcf_records(clinvar_path):match_by_key[item['key']].append(item)
    log(f'Nearby ClinVar allele records: {sum(map(len,match_by_key.values())):,}, unique keys: {len(match_by_key):,}')
    result={}
    for row in rows:
        vid=row['Variant_ID'];is_indel=len(row['Ref'])!=1 or len(row['Alt'])!=1
        data={x:'.' for x in CLIN_FIELDS};data['ClinVar_source']=str(clinvar.name)
        if not is_standard(row):
            data['ClinVar_match_status']='UNSUPPORTED_ALT_SYMBOL';result[vid]=data;continue
        key=original_map.get(vid,vid)
        hits=match_by_key.get(key,[])
        data['ClinVar_exact_Variant_ID']=key
        if not hits:
            data['ClinVar_match_status']='NO_MATCH_AFTER_NORMALIZATION' if ref else ('NO_EXACT_MATCH_INDEL_REVIEW' if is_indel else 'NO_EXACT_MATCH')
            result[vid]=data;continue
        categories=[];vals={k:set() for k in ['CLNSIG','CLNDN','CLNREVSTAT','CLNALLELEID']};ids=set()
        for item in hits:
            ids.add(item['id'])
            for name in vals:
                v=item['info'].get(name,'.')
                if v not in MISSING: vals[name].add(v)
            categories.append(clinvar_class(item['info'].get('CLNSIG','.')))
        data['ClinVar_match_status']='EXACT_MATCH_NORMALIZED' if ref else 'EXACT_MATCH_UNNORMALIZED'
        data['ClinVar_CLNSIG']=';'.join(sorted(vals['CLNSIG'])) or '.'
        data['ClinVar_CLNDN']=';'.join(sorted(vals['CLNDN'])) or '.'
        data['ClinVar_CLNREVSTAT']=';'.join(sorted(vals['CLNREVSTAT'])) or '.'
        data['ClinVar_allele_ID']=';'.join(sorted(vals['CLNALLELEID'])) or '.'
        data['ClinVar_record_ID']=';'.join(sorted(ids)) or '.'
        cats=set(categories)
        if len(cats)>1:cat='MULTIPLE_CLINVAR_CLASSIFICATIONS_REVIEW'
        else:cat=next(iter(cats))
        data['ClinVar_category']=cat
        result[vid]=data
    return result

def load_ndd(p):
    records,_=read_tsv(p);ndd=collections.defaultdict(lambda:{'source':set(),'disease':set()})
    for r in records:
        symbol=str(r.get('Gene','')).strip()
        if not symbol or symbol in MISSING:continue
        for key,col in [('source','Source'),('disease','Disease')]:
            v=str(r.get(col,'')).strip()
            if v not in MISSING:ndd[symbol][key].add(v)
    return ndd

def load_gnomad(prior,current):
    old={}
    if Path(prior).exists():
        for r in read_tsv(prior)[0]:old[r['Variant_ID']]=r
    fresh={}
    if Path(current).exists():
        for r in read_tsv(current)[0]:fresh[r['Variant_ID']]=r
    return old,fresh

def gnomad_row(vid,old,fresh):
    d={k:'.' for k in GN_FIELDS}
    if vid in fresh:
        r=fresh[vid];d.update({
          'gnomAD_screen_source':'B11_14_NEW_24','gnomAD_screen_status':r.get('gnomAD_v4_status','.'),
          'gnomAD_screen_max_AF':r.get('gnomAD_v4_max_observed_AF','.'),
          'gnomAD_screen_decision':r.get('gnomAD_v4_decision','.'),
          'gnomAD_screen_reason':r.get('gnomAD_v4_reason','.')})
    elif vid in old:
        r=old[vid];d.update({
          'gnomAD_screen_source':'B11_EXACT_VARIANT_REUSED',
          'gnomAD_screen_status':r.get('gnomAD_v4_Status','.'),
          'gnomAD_screen_max_AF':r.get('gnomAD_v4_max_AF','.'),
          'gnomAD_screen_decision':r.get('gnomAD_v4_Decision','.'),
          'gnomAD_screen_reason':r.get('gnomAD_v4_exclusion_reason','.')})
    else:
        d['gnomAD_screen_source']='NOT_QUERIED_IN_AVAILABLE_RESULTS';d['gnomAD_screen_status']='PENDING';d['gnomAD_screen_decision']='PENDING'
    return d

def get_family_stats(p):
    groups=collections.defaultdict(dict)
    if not Path(p).is_file():return groups
    with op(p) as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key=row['Variant_ID'];fid=row['Family_ID'];s=row['Strict_screen_status']
            # one row per variant/family expected; combine any unexpected duplicates conservatively
            if fid in groups[key] and groups[key][fid]!=s:raise RuntimeError(f'Conflicting family statuses for {key}, {fid}')
            groups[key][fid]=s
    return groups

def build_gene_table(rows,fam,ndd):
    agg={}
    for r in rows:
        if int(r.get('Affected_QC_ALT_samples') or 0)<=0:continue
        for symbol in genes(r):
            a=agg.setdefault(symbol,{'ids':set(),'families':set(),'strict_families':set(),'inconsistent_families':set(),
                'by_consequence':collections.Counter(),'gn_retain':0,'gn_pending':0,'gn_exclude':0,
                'clin_plp':0,'clin_benign':0,'clin_no_match':0,'de_novo_events':0,'de_novo_variants':0})
            vid=r['Variant_ID'];a['ids'].add(vid)
            for fid,status in fam.get(vid,{}).items():
                if status=='COMPATIBLE_STRICT_SCREEN':a['strict_families'].add(fid);a['families'].add(fid)
                elif status=='INCONSISTENT_WITH_STRICT_RULE':a['inconsistent_families'].add(fid);a['families'].add(fid)
                elif status=='AFFECTED_CARRIER_ONLY':a['families'].add(fid)
            cons=r.get('ANNOVAR_function_reason','.')
            if cons in CONSEQ:a['by_consequence'][cons]+=1
            decision=r['gnomAD_screen_decision']
            if decision in {'RETAIN','RETAIN_FOR_REVIEW_NOT_FOUND'}:a['gn_retain']+=1
            elif decision=='EXCLUDE':a['gn_exclude']+=1
            else:a['gn_pending']+=1
            cc=r['ClinVar_category']
            if cc=='PATHOGENIC_OR_LIKELY_PATHOGENIC':a['clin_plp']+=1
            elif cc=='BENIGN_OR_LIKELY_BENIGN':a['clin_benign']+=1
            if str(r['ClinVar_match_status']).startswith('NO_'):a['clin_no_match']+=1
            dn=int(r.get('Possible_de_novo_children_review') or 0)
            a['de_novo_events']+=dn;a['de_novo_variants']+=int(dn>0)
    fields=['Gene','Unique_variants','Distinct_affected_families_B11_14','Strict_compatible_families',
            'Inconsistent_families','Potential_de_novo_variant_alleles','Potential_de_novo_child_variant_events',
            'Stopgain_variants','Frameshift_deletion_variants','Frameshift_insertion_variants',
            'Stoploss_variants','Splicing_annotated_variants','gnomAD_retained_or_absent_SNVS_for_review',
            'gnomAD_excluded','gnomAD_pending_or_representation_review',
            'ClinVar_PLP_exact_matches','ClinVar_benign_LB_exact_matches','ClinVar_no_match',
            'All_variant_IDs']+NDD_FIELDS
    result=[]
    for gene,a in agg.items():
        info=ndd.get(gene,{})
        d={'Gene':gene,'Unique_variants':len(a['ids']),
           'Distinct_affected_families_B11_14':len(a['families']),
           'Strict_compatible_families':len(a['strict_families']),
           'Inconsistent_families':len(a['inconsistent_families']),
           'Potential_de_novo_variant_alleles':a['de_novo_variants'],
           'Potential_de_novo_child_variant_events':a['de_novo_events'],
           'gnomAD_retained_or_absent_SNVS_for_review':a['gn_retain'],
           'gnomAD_excluded':a['gn_exclude'],
           'gnomAD_pending_or_representation_review':a['gn_pending'],
           'ClinVar_PLP_exact_matches':a['clin_plp'],
           'ClinVar_benign_LB_exact_matches':a['clin_benign'],
           'ClinVar_no_match':a['clin_no_match'],
           'All_variant_IDs':';'.join(sorted(a['ids'])),
           'NDD_overlap':'Yes' if info else 'No',
           'NDD_source':';'.join(sorted(info.get('source',()))) or '.',
           'NDD_disease':';'.join(sorted(info.get('disease',()))) or '.'}
        for consequence in CONSEQ:
            col={'stopgain':'Stopgain_variants','frameshift_deletion':'Frameshift_deletion_variants',
                 'frameshift_insertion':'Frameshift_insertion_variants','stoploss':'Stoploss_variants',
                 'splicing':'Splicing_annotated_variants'}[consequence]
            d[col]=a['by_consequence'][consequence]
        result.append(d)
    result.sort(key=lambda x:(-x['Potential_de_novo_variant_alleles'],-x['Strict_compatible_families'],
                     -x['Unique_variants'],x['Gene']))
    return result,fields

def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('--catalog',type=Path,default=RESULTS/'batch1114_lof_variant_rollup.tsv.gz')
    ap.add_argument('--family-screen',type=Path,default=RESULTS/'batch1114_lof_family_variant_screen.tsv.gz')
    ap.add_argument('--gnomad24',type=Path,default=RESULTS/'gnomad_first24_20261027/gnomad_r4_all_annotated.tsv.gz')
    ap.add_argument('--prior-gnomad',type=Path,default=BASE/'prior_batch11_gnomad_keyed_20261027.tsv.gz')
    ap.add_argument('--additional-gnomad',type=Path,action='append',default=[],
                    help='optional results from subsequent gnomAD runs (repeat to merge multiple batches)')
    ap.add_argument('--ndd',type=Path,default=BASE/'ndd_gene_sources_20261016.tsv.gz')
    ap.add_argument('--clinvar-vcf',type=Path,required=True)
    ap.add_argument('--reference-fasta',type=Path,default=None,
                    help='optional matching indexed GRCh38 FASTA: normalize both datasets before matching')
    ap.add_argument('--outdir',type=Path,default=RESULTS/'batch11_style_ClinVar_NDD_20261027')
    args=ap.parse_args()
    for file in (args.catalog,args.ndd,args.clinvar_vcf):
        if not file.is_file():ap.error('Input missing: '+str(file))
    if args.reference_fasta and not args.reference_fasta.is_file():ap.error('Reference FASTA missing')
    args.outdir.mkdir(parents=True,exist_ok=True)
    rows,fields=read_tsv(args.catalog)
    if len(rows)!=len(set(r['Variant_ID'] for r in rows)):
        raise RuntimeError('Duplicate alleles in input rollup; stop')
    log(f'Read {len(rows)} distinct joint-call annotated ALT alleles')
    ndd=load_ndd(args.ndd)
    old,fresh=load_gnomad(args.prior_gnomad,args.gnomad24)
    additional_ids=set()
    for add in args.additional_gnomad:
        extra,_=read_tsv(add)
        for item in extra:
            vid=item['Variant_ID']
            if item.get('gnomAD_v4_decision') not in ('', '.', 'PENDING', None):
                fresh[vid]=item;additional_ids.add(vid)
    log(f'Loaded {len(ndd)} NDD genes, {len(old)} prior gnomAD decisions and {len(fresh)} newer gnomAD entries')
    clin=clinvar_lookup(rows,args.clinvar_vcf,args.outdir,args.reference_fasta)
    family=get_family_stats(args.family_screen)
    for r in rows:
        vid=r['Variant_ID']
        r.update(clin[vid]);r.update(gnomad_row(vid,old,fresh))
        if vid in additional_ids:r['gnomAD_screen_source']='B11_14_NEW_ADDITIONAL'
        g=genes(r);matched=[ndd[s] for s in g if s in ndd]
        r.update({'NDD_overlap':'Yes' if matched else 'No',
                  'NDD_source':';'.join(sorted(set(v for x in matched for v in x['source']))) or '.',
                  'NDD_disease':';'.join(sorted(set(v for x in matched for v in x['disease']))) or '.'})
    extra=CLIN_FIELDS+GN_FIELDS+NDD_FIELDS
    full_fields=fields+extra
    dump_tsv(args.outdir/'batch1114_ALL_14124_LoF_splicing_review.tsv.gz',rows,full_fields)
    g,ggfields=build_gene_table(rows,family,ndd)
    dump_tsv(args.outdir/'batch1114_Gene_summary_with_NDD.tsv',g,ggfields)
    dump_tsv(args.outdir/'batch1114_NDD_gene_overlaps.tsv',[x for x in g if x['NDD_overlap']=='Yes'],ggfields)
    dump_tsv(args.outdir/'batch1114_exact_ClinVar_PLP_affected_carriers.tsv',[
        r for r in rows if int(r.get('Affected_QC_ALT_samples') or 0)>0 and r['ClinVar_category']=='PATHOGENIC_OR_LIKELY_PATHOGENIC'],full_fields)
    dump_tsv(args.outdir/'batch1114_24_gnomad_ClinVar_NDD_review.tsv',[
        r for r in rows if r['Variant_ID'] in fresh],full_fields)
    stats=collections.Counter((r['ClinVar_category'] for r in rows))
    decisions=collections.Counter(r['gnomAD_screen_decision'] for r in rows)
    clinstats=collections.Counter(r['ClinVar_match_status'] for r in rows)
    with open(args.outdir/'SUMMARY_README.txt','w') as f:
        f.write('Batches 11--14 | Batch 11-style LoF/splicing gene review | '+date()+'\n')
        f.write(f'ClinVar source: {args.clinvar_vcf}; GRCh38\n')
        f.write('Reference normalization: '+(str(args.reference_fasta) if args.reference_fasta else 'NONE - exact raw representation only')+'\n')
        f.write(f'Unique ANNOVAR LoF/splicing ALT alleles: {len(rows)}\n')
        f.write(f'Genes with a confident affected carrier: {len(g)}; NDD overlaps: {sum(x["NDD_overlap"]=="Yes" for x in g)}\n')
        f.write(f'ClinVar exact matches: {sum(v for k,v in clinstats.items() if k.startswith("EXACT"))}\n')
        f.write('ClinVar match breakdown: '+json.dumps(dict(clinstats),sort_keys=True)+'\n')
        f.write('ClinVar significance breakdown (includes not-found): '+json.dumps(dict(stats),sort_keys=True)+'\n')
        f.write('gnomAD decision breakdown (do not treat PENDING as rare): '+json.dumps(dict(decisions),sort_keys=True)+'\n')
        f.write(f'B11 exact-key reuse loaded: {len(old)}, NEW entries loaded: {len(fresh)}, additional: {len(additional_ids)}\n')
        f.write('No ClinVar match does NOT mean benign or unreported after perfect normalization.\n')
        f.write('P/LP classification applies to ClinVar-reported condition(s), NOT necessarily Tourette disorder.\n')
        f.write('A P/LP record with low review status or conflicting submissions merits review, not a pathogenic claim.\n')
        f.write('Family status is an exploratory strict screen; potential de novo is NOT confirmed.\n')
        f.write('gnomAD previously checked exact matches are reused, all other variants are PENDING.\n')
        f.write('All original Yale calls, old analyses, and previous gnomAD checkpoints remain unchanged.\n')
    log('COMPLETE: '+str(args.outdir));print((args.outdir/'SUMMARY_README.txt').read_text())

if __name__=='__main__': main()
