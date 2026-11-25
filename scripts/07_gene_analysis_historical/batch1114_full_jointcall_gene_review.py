#!/usr/bin/env python3
"""Batch 11-14: affected/unaffected, multi-family genes, and gnomAD v4 review.

Python 3.9+, standard library only. Intended for a lab Linux server.
READS existing TSV(.gz) tables; NEVER modifies VCFs, source files or old checkpoints.
Creates separate outputs in --outdir. Only the `query` action contacts gnomAD,
transmitting GRCh38 chromosome/position/REF/ALT, NOT sample IDs or phenotypes.
Use only with approval to query the external resource.

Actions:
  inspect  : display input column names before analyzing
  prepare  : validate input schema, build complete tables + gene-level recurrence
  query    : rate-limited, resumable gnomAD v4 queries for PRIORITY missing alleles
  export   : merge checkpointed gnomAD results, rebuild gene-level reports

Required --variants: 14,124-row annotated LoF/splicing review or equivalent.
Required --events  : FAMILY-VARIANT screen with exact allele, family ID and explicit
                     affected/unaffected ALT genotype counts (or affected IDs).
Optional --rollup : per-allele cohort counts if missing from --variants.
Optional --extra  : e.g., completed 24-variant gnomAD/ClinVar annotations.

Interpretation:
* >=2 distinct ALT variants across >=2 distinct affected-carrier families is
  an EXPLORATORY recurrence heuristic, NOT proof of segregation/association.
* Sum of carriers across variants = carrier-variant observations, not unique
  affected individuals per gene.
* gnomAD NOT_FOUND / no AF / pending / indel mismatch != AF zero or 'rare'.
* gnomAD `flags` describe annotation/interpretation, not necessarily site QC;
  only exome/genome `filters` are used for site-quality review.
* Existing historical EXCLUDE decisions are reported separately, since they
  may include annotation flags (e.g. lc_lof) rather than site-quality failures.
* Unnormalized indels require separate GRCh38 normalization and manual review.
"""
import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import math
import os
import re
import sqlite3
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib import request, error

PROJECT = Path(os.environ.get('PROJECT_ROOT', '/path/to/project'))
DEF_VARIANTS = PROJECT/'filtering/family_batches11_14/lof_joint_20261026/results/batch11_style_ClinVar_NDD_20261027/batch1114_ALL_14124_LoF_splicing_review.tsv.gz'
DEF_EVENTS = PROJECT/'filtering/family_batches11_14/lof_joint_20261026/results/batch1114_lof_family_variant_screen.tsv.gz'
DEF_ROLLUP = PROJECT/'filtering/family_batches11_14/lof_joint_20261026/results/batch1114_lof_variant_rollup.tsv.gz'
DEF_EXTRA = PROJECT/'filtering/family_batches11_14/lof_joint_20261026/results/batch11_style_ClinVar_NDD_20261027/batch1114_24_gnomad_ClinVar_NDD_review.tsv'
DEF_OUT = PROJECT/'filtering/family_batches11_14/full_multifamily_review_20261027'
API = 'https://gnomad.broadinstitute.org/api'
AF_CUTOFF = 0.001
MISSING = {'', '.', 'na', 'nan', 'none', 'null', 'n/a'}
QC_FIELDS = ['Affected_QC_ALT_samples','Unaffected_QC_ALT_samples']
QUERY = '''query VariantReview($id: String!, $ds: DatasetId!) {
  variant(variantId: $id, dataset: $ds) {
    variant_id flags
    exome { af filters populations { id ac an } }
    genome { af filters populations { id ac an } }
  }
}'''


def timestamp():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def msg(s):
    print(f'[{timestamp()}] {s}',flush=True)


def norm(s):
    return re.sub('[^a-z0-9]','',str(s).lower())


def open_text(p,mode='rt'):
    return gzip.open(p,mode,encoding='utf-8-sig' if 'r' in mode else 'utf-8',newline='') if str(p).endswith('.gz') else open(p,mode,encoding='utf-8-sig' if 'r' in mode else 'utf-8',newline='')


def header(p):
    if not p.is_file():
        raise ValueError('Missing file: '+str(p))
    with open_text(p) as f:
        h=next(csv.reader(f,delimiter='\t'),[])
    if len(h)<2: raise ValueError('Not a tab-delimited header: '+str(p))
    return h


def rows(p):
    with open_text(p) as f:
        rdr=csv.DictReader(f,delimiter='\t')
        for i,r in enumerate(rdr,2):
            if None in r: raise ValueError(f'Extra TSV fields in {p} line {i}')
            yield i,r


def col(fields,*aliases,required=False):
    d={norm(v):v for v in fields}
    for a in aliases:
        if norm(a) in d:return d[norm(a)]
    if required:
        raise ValueError('Missing required field '+str(aliases)+'; available: '+str(fields))
    return None


def get(r,k):
    return str((r.get(k) if k else '') or '').strip()


def missing(s):
    return str(s).strip().lower() in MISSING


def intval(s,where):
    if missing(s): return None
    try:
        f=float(s)
        if not f.is_integer() or f<0: raise ValueError()
        return int(f)
    except (ValueError,TypeError):raise ValueError(f'Invalid nonnegative integer at {where}: {s!r}')


def afval(s):
    if missing(s):return None
    try:
        v=float(s)
        return v if math.isfinite(v) and 0<=v<=1 else None
    except (ValueError,TypeError):return None


def variant_id(r,h,where):
    k=col(h,'Variant_ID','VariantID','variant_key','Allele_ID')
    vid=get(r,k)
    if vid and not missing(vid):
        m=re.fullmatch(r'(?:chr)?([0-9]+|X|Y|M|MT)[-:]([0-9]+)[-:]([A-Za-z*<>.]+)[-:]([A-Za-z*<>.]+)',vid,re.I)
        if not m:raise ValueError(f'Invalid allele ID at {where}: {vid[:120]}')
        c,p,ref,alt=m.groups()
    else:
        c=get(r,col(h,'Chr','CHROM','#CHROM','chromosome',required=True))
        p=get(r,col(h,'Pos','POS','Start','position',required=True))
        ref=get(r,col(h,'Ref','REF',required=True))
        alt=get(r,col(h,'Alt','ALT',required=True))
    c=c.upper().removeprefix('CHR') if sys.version_info>=(3,9) else re.sub('^CHR','',c.upper())
    if not p.isdecimal() or int(p)<1 or not c or not ref or not alt:
        raise ValueError(f'Invalid coordinates at {where}')
    return f'{c}-{int(p)}-{ref.upper()}-{alt.upper()}'


def gene_symbols(s):
    if missing(s):return set()
    return {x.strip() for x in re.split('[,;|]',str(s)) if not missing(x)}


def choose_count(r,h,affected,where):
    if affected:
        ck=col(h,'Affected_QC_ALT_samples','Affected_ALT_samples','Affected_ALT','n_affected_carriers','Affected_carrier_count','Affected_Carrier_Count','Affected_carriers_count')
        ids=col(h,'Affected_Carrier_Samples','Affected_QC_ALT_sample_ids','Affected_ALT_sample_ids','Affected_Carrier_IDs')
    else:
        ck=col(h,'Unaffected_QC_ALT_samples','Unaffected_ALT_samples','Unaffected_ALT','n_unaffected_carriers','Unaffected_carrier_count','Unaffected_Carrier_Count','Unaffected_carriers_count')
        ids=col(h,'Unaffected_Carrier_Samples','Unaffected_QC_ALT_sample_ids','Unaffected_ALT_sample_ids','Unaffected_Carrier_IDs')
    if ck:
        return intval(get(r,ck),where+' '+ck)
    if ids:
        value=get(r,ids)
        return len([x for x in re.split('[,;|]',value) if not missing(x)]) if not missing(value) else 0
    return None


def parse_previous(r,h):
    status=get(r,col(h,'gnomAD_screen_status','gnomAD_v4_status','gnomAD_status')).upper()
    decision=get(r,col(h,'gnomAD_screen_decision','gnomAD_v4_decision','gnomAD_decision')).upper()
    a=get(r,col(h,'gnomAD_screen_max_AF','gnomAD_v4_max_observed_AF','gnomAD_max_AF'))
    af=afval(a)
    reason=get(r,col(h,'gnomAD_screen_reason','gnomAD_v4_reason','Screen_note'))
    if decision=='EXCLUDE':
        return {'status':'PRIOR_EXCLUDE_REVIEW','af':af,'reason':'Historical exclusion: '+reason}
    if decision in ('REVIEW_REPRESENTATION','REVIEW_NO_FREQUENCY'):
        return {'status':decision,'af':af,'reason':reason}
    if status in ('NOT_FOUND','NOTFOUND') or decision=='RETAIN_FOR_REVIEW_NOT_FOUND':
        return {'status':'NOT_FOUND_UNRESOLVED','af':None,'reason':'Exact-key gnomAD lookup did not find allele; AF unknown'}
    if decision=='RETAIN' and af is not None and af<AF_CUTOFF:
        return {'status':'PRIOR_RARE_RETAIN','af':af,'reason':'Earlier frequency/quality screen retained'}
    if af is not None and af>=AF_CUTOFF:
        return {'status':'PRIOR_AF_GE_0.001','af':af,'reason':'Earlier observed AF >= cutoff'}
    if not missing(status) and status in ('FOUND',) and af is not None:
        return {'status':'PRIOR_AF_UNVERIFIED','af':af,'reason':'Observed AF but no verified quality/retain decision'}
    return None


def score_prior(r):
    return {'PRIOR_RARE_RETAIN':8,'PRIOR_AF_GE_0.001':7,'PRIOR_EXCLUDE_REVIEW':6,'NOT_FOUND_UNRESOLVED':5,'REVIEW_REPRESENTATION':5,'REVIEW_NO_FREQUENCY':4,'PRIOR_AF_UNVERIFIED':3}.get(r['status'],0)


def digest(p):
    # read-once hash to protect against accidentally mixing checkpoint inputs
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def write_tsv(path,fieldnames,data,compressed=False):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.working')
    if compressed:
        f=gzip.open(tmp,'wt',encoding='utf-8',newline='')
    else:f=open(tmp,'w',encoding='utf-8',newline='')
    with f:
        w=csv.DictWriter(f,fieldnames=fieldnames,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader()
        for item in data:w.writerow(item)
    os.replace(tmp,path)
    msg(f'WROTE {path}')


def inspect(a):
    for label,path in [('variants',a.variants),('events',a.events),('rollup',a.rollup),('extra',a.extra)]:
        if path and path.exists():
            print(f'\n[{label}] {path}\n'+'\n'.join(header(path)))
        elif label in ('variants','events'):
            print(f'\n[{label}] MISSING: {path}')


def process(a):
    vh=header(a.variants);eh=header(a.events)
    gene_col=col(vh,'Gene_refGene','Gene.refGene','Gene','Genes','Gene_Name','gene_symbol',required=True)
    event_family=col(eh,'Family_ID','Family','FamilyID','family_key','family_id',required=True)
    # Explicit family-level carrier counts or carrier lists are essential.
    if not (col(eh,'Affected_QC_ALT_samples','Affected_ALT_samples','Affected_ALT','Affected_carrier_count','n_affected_carriers','Affected_Carrier_Count','Affected_carriers_count') or col(eh,'Affected_Carrier_Samples','Affected_QC_ALT_sample_ids','Affected_ALT_sample_ids','Affected_Carrier_IDs')):
        raise ValueError('Events require an explicit affected-carrier count or ID column. Run inspect and share headers; cannot infer carriers from row existence.')
    if not (col(eh,'Unaffected_QC_ALT_samples','Unaffected_ALT_samples','Unaffected_ALT','Unaffected_carrier_count','n_unaffected_carriers','Unaffected_Carrier_Count','Unaffected_carriers_count') or col(eh,'Unaffected_Carrier_Samples','Unaffected_QC_ALT_sample_ids','Unaffected_ALT_sample_ids','Unaffected_Carrier_IDs')):
        raise ValueError('Events require an explicit unaffected-carrier count or ID column. Run inspect; not safe to assume absent field means zero.')

    V={}; families=defaultdict(dict)
    for n,r in rows(a.variants):
        k=variant_id(r,vh,f'{a.variants}:{n}')
        genes=gene_symbols(get(r,gene_col))
        if not genes:continue
        if k in V:
            raise ValueError(f'Duplicate allele in variant catalog {k}; inspect row {n} before deduplicating')
        V[k]={'vid':k,'genes':genes,'affected':choose_count(r,vh,True,f'{a.variants}:{n}'),
              'unaffected':choose_count(r,vh,False,f'{a.variants}:{n}'),
              'consequence':get(r,col(vh,'ExonicFunc_refGene','ExonicFunc.refGene','Consequence','ANNOVAR_function_reason')) or get(r,col(vh,'Func_refGene','Func.refGene')),
              'de_novo':intval(get(r,col(vh,'Possible_de_novo_children_review')),f'possible de novo {k}') or 0,
              'clinvar':get(r,col(vh,'ClinVar_CLNSIG','CLNSIG','ClinVar_significance')),
              'ndd':get(r,col(vh,'NDD_overlap','NDD_Gene_Overlap')),
              'prior':parse_previous(r,vh)}
    if not V:raise ValueError('No usable variant annotations in '+str(a.variants))
    msg(f'Loaded {len(V)} annotated distinct alleles')

    if a.rollup and not a.rollup.exists():
        msg(f'Optional rollup file not found: {a.rollup}; using catalog/family-event carrier counts')
    if a.rollup and a.rollup.exists():
        rh=header(a.rollup)
        matched=0
        for n,r in rows(a.rollup):
            k=variant_id(r,rh,f'{a.rollup}:{n}')
            if k not in V:continue
            matched+=1
            for field,affected in [('affected',True),('unaffected',False)]:
                x=choose_count(r,rh,affected,f'{a.rollup}:{n}')
                if V[k][field] is None and x is not None:V[k][field]=x
            x=intval(get(r,col(rh,'Possible_de_novo_children_review')),f'possible de novo {k}')
            if x is not None:V[k]['de_novo']=max(V[k]['de_novo'],x)
        msg(f'Rollup exact allele matches: {matched}')

    if a.extra and not a.extra.exists():
        msg(f'Optional annotation file not found: {a.extra}; previous 24-variant gnomAD results cannot be imported from it')
    if a.extra and a.extra.exists():
        xh=header(a.extra);matched=0
        for n,r in rows(a.extra):
            k=variant_id(r,xh,f'{a.extra}:{n}')
            if k not in V:continue
            matched+=1
            prior=parse_previous(r,xh)
            if prior and (not V[k]['prior'] or score_prior(prior)>score_prior(V[k]['prior'])):V[k]['prior']=prior
            dn=intval(get(r,col(xh,'Possible_de_novo_children_review')),f'possible de novo {k}')
            if dn is not None:V[k]['de_novo']=max(V[k]['de_novo'],dn)
            if missing(V[k]['clinvar']):V[k]['clinvar']=get(r,col(xh,'ClinVar_CLNSIG','CLNSIG'))
            if missing(V[k]['ndd']):V[k]['ndd']=get(r,col(xh,'NDD_overlap'))
        msg(f'Additional exact-allele annotations matched: {matched}')

    matched=0;unmatched=0; event_rows=0
    for n,r in rows(a.events):
        event_rows+=1
        k=variant_id(r,eh,f'{a.events}:{n}')
        fam=get(r,event_family)
        if missing(fam):raise ValueError(f'Missing family ID for variant {k} at event row {n}')
        affected=choose_count(r,eh,True,f'{a.events}:{n}')
        unaffected=choose_count(r,eh,False,f'{a.events}:{n}')
        if affected is None or unaffected is None:
            raise ValueError(f'Missing carrier count in family-event row {n}; cannot assume zero')
        if k not in V:
            unmatched+=1
            continue
        matched+=1
        pair=(affected,unaffected)
        old=families[k].get(fam)
        if old is not None and old!=pair:
            raise ValueError(f'Conflicting duplicate family-allele event counts {k} family {fam}; verify source')
        families[k][fam]=pair
    if matched==0:raise ValueError('No event alleles matched the variant catalog; verify genome build and ALT splitting')
    msg(f'Matched family-event rows: {matched}; unmatched: {unmatched}; total: {event_rows}')

    inconsistencies=[]
    for k,v in V.items():
        fa=families.get(k,{})
        af=sum(x[0] for x in fa.values());un=sum(x[1] for x in fa.values())
        if v['affected'] is None:v['affected']=af if fa else None
        elif fa and v['affected']!=af:inconsistencies.append((k,'affected',v['affected'],af))
        if v['unaffected'] is None:v['unaffected']=un if fa else None
        elif fa and v['unaffected']!=un:inconsistencies.append((k,'unaffected',v['unaffected'],un))
        v['affected_families']={f for f,(a,u) in fa.items() if a>0}
        v['unaffected_families']={f for f,(a,u) in fa.items() if u>0}
        v['family_events_observed']=bool(fa)
        v['counts_complete']=(v['affected'] is not None and v['unaffected'] is not None)
    if inconsistencies:
        # Preserve source cohort-level carrier counts; family-screening genotype subsets may differ.
        msg(f'WARNING: {len(inconsistencies)} cohort-vs-family event count discrepancies; see discrepancy report; do not treat mismatched counts as validated.')
        write_tsv(a.outdir/'carrier_count_discrepancies.tsv',['Variant_ID','Group','Cohort_count','Family_event_sum'],
                  ({'Variant_ID':k,'Group':c,'Cohort_count':x,'Family_event_sum':y} for k,c,x,y in inconsistencies))
    return V


def gene_summary(V):
    G=defaultdict(lambda:{'variants':set(),'families':set(),'unaff_families':set(),'rare':set(),'unknown':set(),'excluded':set(),'dn':set(),'plp':set(),'ndd':set(),'carrier_events_aff':0,'carrier_events_unaff':0})
    for k,v in V.items():
        if v['affected'] is None or v['unaffected'] is None or v['affected']<1:continue
        for gene in v['genes']:
            o=G[gene];o['variants'].add(k)
            o['families'].update(v['affected_families']);o['unaff_families'].update(v['unaffected_families'])
            o['carrier_events_aff']+=v['affected'];o['carrier_events_unaff']+=v['unaffected']
            status=v['gnomad']['status']
            if status in ('PRIOR_RARE_RETAIN','GNOMAD_RARE_PASS'):o['rare'].add(k)
            elif status in ('PRIOR_AF_GE_0.001','PRIOR_EXCLUDE_REVIEW','GNOMAD_AF_GE_0.001','GNOMAD_SITE_FILTER_REVIEW'):o['excluded'].add(k)
            else:o['unknown'].add(k)
            if v['de_novo']>0:o['dn'].add(k)
            clin=v['clinvar'].lower()
            if 'pathogenic' in clin and not any(z in clin for z in ('benign','conflict','uncertain')):o['plp'].add(k)
            if v['ndd'].lower() in ('yes','y','true','1'):o['ndd'].add(k)
    return G


def gn_result_from_api(payload):
    v=(payload.get('data') or {}).get('variant')
    if v is None:return {'status':'NOT_FOUND_UNRESOLVED','af':None,'reason':'No exact gnomAD r4 allele match; AF unknown','filters':'.','variant_flags':'.'}
    fs=[];afs=[];popmax=None;ex='';gen=''
    for tag in ('exome','genome'):
        call=v.get(tag) or {}
        af=afval(call.get('af'))
        if tag=='exome':ex='' if af is None else str(af)
        else:gen='' if af is None else str(af)
        if af is not None:afs.append(af)
        flags=call.get('filters') or []
        if isinstance(flags,str):flags=[flags]
        fs.extend(tag+':'+str(x) for x in flags if str(x).strip() and str(x).upper() not in ('PASS','.'))
        for pop in call.get('populations') or []:
            try:
                an=float(pop.get('an'));num=float(pop.get('ac'))
                ratio=num/an if an>0 and num>=0 else None
            except (TypeError,ValueError,ZeroDivisionError):ratio=None
            if ratio is not None and 0<=ratio<=1:
                popmax=ratio if popmax is None else max(popmax,ratio)
    if popmax is not None:afs.append(popmax)
    overall=max(afs) if afs else None
    interpretive=v.get('flags') or []
    if isinstance(interpretive,str):interpretive=[interpretive]
    if fs:status='GNOMAD_SITE_FILTER_REVIEW';reason='Non-PASS exome/genome filter(s)'
    elif overall is None:status='GNOMAD_NO_USABLE_AF';reason='Found but frequency not available'
    elif overall>=AF_CUTOFF:status='GNOMAD_AF_GE_0.001';reason='Maximum AF >= 0.001'
    else:status='GNOMAD_RARE_PASS';reason='Maximum observed AF <0.001 and no non-PASS site filters'
    return {'status':status,'af':overall,'reason':reason,'filters':';'.join(fs) or '.',
            'variant_flags':';'.join(str(z) for z in interpretive) or '.',
            'exome_af':ex,'genome_af':gen,'population_max_af':str(popmax) if popmax is not None else ''}


def api_query(vid):
    body=json.dumps({'query':QUERY,'variables':{'id':vid,'ds':'gnomad_r4'}}).encode('utf-8')
    req=request.Request(API,data=body,headers={'Content-Type':'application/json','User-Agent':'WES-candidate-review/1.0'},method='POST')
    try:
        with request.urlopen(req,timeout=35) as response:
            data=json.load(response)
    except error.HTTPError as e:
        raise RuntimeError(f'HTTP {e.code}: {e.read(500).decode("utf-8","replace")}')
    if data.get('errors'):
        raise RuntimeError('GraphQL errors: '+str(data['errors'])[:400])
    if 'data' not in data or 'variant' not in (data.get('data') or {}):
        raise RuntimeError('Unexpected gnomAD API structure; do not checkpoint')
    result=gn_result_from_api(data)
    ref,alt=vid.split('-',3)[2:]
    if result['status']=='NOT_FOUND_UNRESOLVED' and (len(ref)!=1 or len(alt)!=1):
        result['status']='REVIEW_INDEL_REPRESENTATION'
        result['reason']='Indel has no exact gnomAD match; normalize against GRCh38 before interpreting population frequency'
    return result


def api_ready(vid):
    p=vid.split('-',3)
    return len(p)==4 and p[0] in ({str(i) for i in range(1,23)}|{'X','Y'}) and bool(re.fullmatch('[ACGT]+',p[2])) and bool(re.fullmatch('[ACGT]+',p[3]))


def merge_checkpoint(V,out):
    db=out/'gnomad_priority.sqlite3'
    saved={}
    if db.exists():
        with sqlite3.connect(db) as c:
            saved={k:json.loads(t) for k,t in c.execute('SELECT vid, result_json FROM results')}
    for k,v in V.items():
        prior=v['prior'] or {'status':'PENDING','af':None,'reason':'Not queried or no verified lookup'}
        if k in saved:
            v['gnomad']={**saved[k],'source':'NEW_GNOMAD_R4'}
        else:v['gnomad']={**prior,'source':'EXISTING_ANNOTATION' if v['prior'] else 'MISSING'}
        if not api_ready(k) and v['gnomad']['status']=='PENDING':
            v['gnomad']={'status':'UNSUPPORTED_ALLELE','af':None,'reason':'Symbolic/spanning deletion or nonstandard allele; no direct GraphQL lookup','source':'NOT_QUERIED'}


def tables(a,V):
    merge_checkpoint(V,a.outdir)
    G=gene_summary(V)
    recurrent={gene for gene,g in G.items() if len(g['variants'])>=a.min_variants and len(g['families'])>=a.min_families}
    outvariants=[]
    for k,v in V.items():
        z=v['gnomad'];af=z.get('af')
        outvariants.append({'Variant_ID':k,'Gene':';'.join(sorted(v['genes'])),'Consequence':v['consequence'],
          'Affected_QC_ALT_samples':v['affected'] if v['affected'] is not None else '',
          'Unaffected_QC_ALT_samples':v['unaffected'] if v['unaffected'] is not None else '',
          'Carrier_counts_complete':int(v['counts_complete']),
          'Distinct_affected_families':len(v['affected_families']),'Distinct_unaffected_families':len(v['unaffected_families']),
          'Possible_de_novo_events_PRELIMINARY':v['de_novo'],'ClinVar_CLNSIG':v['clinvar'],'NDD_overlap':v['ndd'],
          'Member_of_multifamily_gene':int(bool(v['genes'] & recurrent)),
          'gnomAD_v4_status':z['status'],'gnomAD_v4_max_observed_AF':str(af) if af is not None else '',
          'gnomAD_v4_source':z.get('source',''),'gnomAD_v4_reason':z.get('reason',''),
          'gnomAD_v4_filters':z.get('filters',''),'gnomAD_v4_variant_interpretive_flags':z.get('variant_flags','')})
    vf=list(outvariants[0])
    write_tsv(a.outdir/'jointcall_all_variants_annotated.tsv.gz',vf,outvariants,True)
    genelist=[]
    for gene,g in G.items():
        genelist.append({'Gene':gene,'Unique_variants_with_affected_carriers':len(g['variants']),
           'Distinct_families_with_affected_carriers':len(g['families']),
           'Distinct_families_with_unaffected_carriers':len(g['unaff_families']),
           'Affected_carrier_variant_observations_NOT_UNIQUE_PEOPLE':g['carrier_events_aff'],
           'Unaffected_carrier_variant_observations_NOT_UNIQUE_PEOPLE':g['carrier_events_unaff'],
           'Rarity_confirmed_variants':len(g['rare']),'Population_pending_or_unresolved_variants':len(g['unknown']),
           'Previously_excluded_or_common_variants':len(g['excluded']),
           'Preliminary_possible_de_novo_variant_alleles':len(g['dn']),
           'ClinVar_PLP_annotated_variant_alleles':len(g['plp']),
           'NDD_overlapping_variant_alleles':len(g['ndd']),
           'Meets_2_variants_2_families':int(gene in recurrent)})
    if not genelist:raise ValueError('No affected-carrier genes with complete counts were found; inspect field mappings and data')
    genelist.sort(key=lambda r:(-r['Meets_2_variants_2_families'],-r['Distinct_families_with_affected_carriers'],
                   -r['Unique_variants_with_affected_carriers'],r['Gene']))
    gf=list(genelist[0])
    write_tsv(a.outdir/'jointcall_all_genes_summary.tsv',gf,genelist)
    write_tsv(a.outdir/'jointcall_multifamily_genes.tsv',gf,(g for g in genelist if g['Meets_2_variants_2_families']))
    query_priority=[]
    # Do not assume a missing frequency implies the variant is rare.
    for k,v in V.items():
        if v['affected'] is None or v['unaffected'] is None or v['affected']<1 or v['gnomad']['status'] not in ('PENDING','PRIOR_AF_UNVERIFIED'):continue
        if not api_ready(k):continue
        recurrent_members=v['genes'] & recurrent
        if a.query_scope=='priority' and not recurrent_members and v['de_novo']==0:continue
        # de novo is preliminary and doesn't confer pathogenicity.
        ranking=(100 if v['de_novo'] else 0)+(50 if recurrent_members else 0)+min(30,5*len(v['affected_families']))
        ranking+=10 if v['unaffected']==0 else 0
        query_priority.append({'Variant_ID':k,'Gene':';'.join(sorted(v['genes'])),
          'Priority_points_NOT_CLINICAL':ranking,'Reason':('POSSIBLE_DENOVO;' if v['de_novo'] else '')+('MULTIFAMILY_GENE' if recurrent_members else ''),
          'Affected_carriers':v['affected'],'Unaffected_carriers':v['unaffected'],
          'Distinct_affected_families':len(v['affected_families'])})
    query_priority.sort(key=lambda r:(-r['Priority_points_NOT_CLINICAL'],-r['Affected_carriers'],r['Variant_ID']))
    qh=list(query_priority[0]) if query_priority else ['Variant_ID','Gene','Priority_points_NOT_CLINICAL','Reason','Affected_carriers','Unaffected_carriers','Distinct_affected_families']
    # A fixed queue is essential for a resumable query checkpoint.
    queuefile=a.outdir/'gnomad_priority_queue.tsv'
    if not queuefile.exists():
        write_tsv(queuefile,qh,query_priority)
    else:
        write_tsv(a.outdir/'gnomad_still_pending.tsv',qh,query_priority)
    status=Counter(v['gnomad']['status'] for v in V.values())
    report=[f'Generated: {timestamp()}',f'Input distinct alleles: {len(V)}',f'Genes with >=1 confidently affected-carried ALT: {len(G)}',
      f'Genes with >={a.min_variants} distinct variants in >={a.min_families} affected-carrier families: {len(recurrent)}',
      f'Current query scope: {a.query_scope}; queryable variants without verified gnomAD screening: {len(query_priority)}',
      f'Variants with unavailable cohort/family carrier count evidence: {sum(not v["counts_complete"] for v in V.values())}',
      'Note: The candidate queue is a focused subset. It does NOT imply other variants were fully gnomAD screened.',
      'Note: Cohort-level affected/unaffected allele carriers are observational; not a formal case-control test.',
      'Note: Unknown/not-found population frequency != AF zero; indels require normalization review.',
      'Note: Gene-level carrier observations are not counts of distinct people.',
      'Note: Recurrent/ClinVar/de novo annotations are preliminary, not evidence of TD pathogenicity.']
    report += [f'gnomAD status {s}: {n}' for s,n in sorted(status.items())]
    (a.outdir/'SUMMARY.txt').write_text('\n'.join(report)+'\n',encoding='utf-8')
    for line in report:print(line)
    print('\nTop 15 recurrent genes:')
    for g in [z for z in genelist if z['Meets_2_variants_2_families']][:15]:
        print(g['Gene'],g['Unique_variants_with_affected_carriers'],g['Distinct_families_with_affected_carriers'],sep='\t')


def write_manifest(a):
    inputs={s:{'path':str(p.resolve()),'sha256':digest(p)} for s,p in
       [('variants',a.variants),('events',a.events),('rollup',a.rollup),('extra',a.extra)] if p is not None and p.is_file()}
    data={'inputs':inputs,'min_variants':a.min_variants,'min_families':a.min_families,'query_scope':a.query_scope,'reference':'GRCh38','dataset':'gnomad_r4'}
    manifest=a.outdir/'source_manifest.json'
    if manifest.exists():
        original=json.loads(manifest.read_text())
        if original!=data:
            raise ValueError('Existing output directory belongs to different source/checksum/criteria. Choose a NEW --outdir; old checkpoint protected.')
    else:
        a.outdir.mkdir(parents=True,exist_ok=True)
        manifest.write_text(json.dumps(data,indent=2)+'\n')


def load_prioritized(a):
    f=a.outdir/'gnomad_priority_queue.tsv'
    if not f.exists():raise ValueError('Run prepare first; missing '+str(f))
    return [r['Variant_ID'] for _,r in rows(f)]


def query_gnomad(a):
    if a.sleep<5 or a.max_new<1 or a.hours<=0:
        raise ValueError('Require --sleep >=5, --max-new >=1 and --hours >0')
    queue=load_prioritized(a)
    if not queue:
        msg('Queue empty; nothing to query');return
    db=a.outdir/'gnomad_priority.sqlite3';lock=a.outdir/'gnomad_query.lock'
    # Avoid multiple simultaneous processes using the same checkpoint.
    with lock.open('a+') as lf:
        try:fcntl.flock(lf,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Another query process is already running')
        with sqlite3.connect(db) as c:
            c.execute('CREATE TABLE IF NOT EXISTS metadata (k TEXT PRIMARY KEY, v TEXT)')
            c.execute('CREATE TABLE IF NOT EXISTS results (vid TEXT PRIMARY KEY, result_json TEXT, queried_at TEXT)')
            queue_hash=digest(a.outdir/'gnomad_priority_queue.tsv')
            existing=c.execute('SELECT v FROM metadata WHERE k="queue_sha256"').fetchone()
            if existing and existing[0]!=queue_hash:
                raise ValueError('Queue changed since checkpoint creation. Use new --outdir.')
            if not existing:c.execute('INSERT INTO metadata VALUES (?,?)',('queue_sha256',queue_hash))
            c.commit()
            done={v for (v,) in c.execute('SELECT vid FROM results')}
        todo=[v for v in queue if v not in done][:a.max_new]
        msg(f'Start gnomAD priority queries: {len(done)} checkpointed; {len(todo)} new requested; {len(queue)} in queue')
        start=time.monotonic();successes=0
        for i,vid in enumerate(todo):
            if time.monotonic()-start>a.hours*3600:
                msg('Time limit reached; checkpoint saved');break
            try:
                result=api_query(vid)
            except Exception as exc:
                msg(f'STOPPED on {vid}: {exc}; completed lookups are safely checkpointed')
                break
            with sqlite3.connect(db) as c:
                c.execute('INSERT OR IGNORE INTO results VALUES (?,?,?)',(vid,json.dumps(result),timestamp()))
                c.commit()
            successes+=1
            msg(f'gnomAD {successes}/{len(todo)} {vid} -> {result["status"]}')
            if i<len(todo)-1:time.sleep(a.sleep)
        msg(f'New completed lookups: {successes}. Run export to update reports; resume query anytime.')


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('action',choices=['inspect','prepare','query','export'])
    p.add_argument('--variants',type=Path,default=DEF_VARIANTS)
    p.add_argument('--events',type=Path,default=DEF_EVENTS)
    p.add_argument('--rollup',type=Path,default=DEF_ROLLUP)
    p.add_argument('--extra',type=Path,default=DEF_EXTRA)
    p.add_argument('--outdir',type=Path,default=DEF_OUT)
    p.add_argument('--min-variants',type=int,default=2)
    p.add_argument('--min-families',type=int,default=2)
    p.add_argument('--query-scope',choices=['priority','all-affected'],default='priority',
                   help='Priority: de novo and recurrent genes; all-affected: all affected-carrier variants missing verified gnomAD.')
    p.add_argument('--max-new',type=int,default=25)
    p.add_argument('--sleep',type=float,default=5.)
    p.add_argument('--hours',type=float,default=2.)
    a=p.parse_args()
    if a.action=='inspect':inspect(a);return
    if a.action=='query':
        if not (a.outdir/'source_manifest.json').exists():raise ValueError('Run prepare first')
        query_gnomad(a);return
    if a.min_variants<2 or a.min_families<2:raise ValueError('Use --min-variants >=2 and --min-families >=2')
    for f in [a.variants,a.events]:
        if not f.exists():raise ValueError('Missing required input '+str(f)+'; use inspect and update --variants / --events')
    write_manifest(a)
    V=process(a)
    tables(a,V)

if __name__=='__main__':
    try:main()
    except Exception as e:
        print('STOPPED SAFELY: '+str(e),file=sys.stderr,flush=True)
        sys.exit(2)
