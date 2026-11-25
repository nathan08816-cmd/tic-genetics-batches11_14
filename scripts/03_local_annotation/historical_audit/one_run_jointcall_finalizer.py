#!/usr/bin/env python3
"""Historical QC/annotation/internal-AF/segregation audit generator.

No external API calls. No changes to input data or existing directories.
Produces concise 9-column gene and 15-column allele TSVs, a full QC audit,
family summary, and a methods/results report.
Existing gnomAD 'prior rare' results remain REVIEW until independently verified.

Requires a 14,124-allele exact ANNOVAR TSV, the original filtered split-ALT VCF,
original 02_lof_family_screen.py, October 14 phenotype file, sample QC and
unrelated subset, and the merged gnomAD/family table. Requires lab-only inputs.
"""
from __future__ import annotations
import argparse, collections, csv, gzip, hashlib, importlib.util, json, math, os, re, sys, textwrap
from datetime import datetime, timezone
from pathlib import Path

# Historical intermediate-audit stage. Its original gnomAD/API and ANNOVAR
# frequency decisions are superseded by the local v4.1 filtering stage.
# Every path can be overridden with a command-line argument.
ROOT = Path(os.environ.get('PROJECT_ROOT', '/path/to/TS_ANNOVAR_analysis'))
SUB = ROOT / 'filtering/family_batches11_14'
DATA = Path(os.environ.get('SOURCE_DATA_DIR', '/path/to/jointcall_inputs'))
DEF_MERGED = ROOT / 'JointCall_Configured_Local_gnomAD_20261110/all_jointcall_variants_merged.tsv.gz'
DEF_ANNOT = SUB / 'lof_joint_20261026/results/batch11_style_ClinVar_NDD_20261027/batch1114_ALL_14124_LoF_splicing_review.tsv.gz'
DEF_VCF = SUB / 'lof_joint_20261026/results/batch1114_lof_QCsite_candidates.with_genotypes.vcf.gz'
DEF_QC = DATA / 'sample_qc/sample_qc.tsv'
DEF_UNREL = DATA / 'sample_qc/unrelated_samples.txt'
DEF_PHENO = ROOT / 'family_structure/inputs/phenotypes.xlsx'
DEF_ORIG = SUB / 'scripts/02_lof_family_screen.py'
DEF_GNOMAD = ROOT / 'JointCall_Configured_Local_gnomAD_20261110/gnomad_v4_results.tsv'
DEF_EVENTS = ROOT / 'JointCall_Configured_Local_gnomAD_20261110/all_family_carrier_events.tsv.gz'
MISSING = {'', '.', 'na', 'n/a', 'nan', 'none', 'null', 'not queried'}
GENES = ['Gene','Candidate_Variants','Variant_Types','Affected_Carriers','Unaffected_Carriers',
         'Affected_Families','Strict_Segregation_Variants','Lowest_gnomAD_AF','ClinVar_Summary']
ALLELES = ['Variant_ID','Gene','Function','gnomAD_AF','gnomAD_QC','Internal_AF',
           'Affected_Carriers','Unaffected_Carriers','Affected_Families','Strict_Families',
           'Inconsistent_Families','Segregation','ClinVar','Filter_Status','Filter_Reason']


def log(s): print(f'[{datetime.now().strftime("%H:%M:%S")}] {s}',flush=True)
def require(path):
    if not path.is_file() or path.stat().st_size == 0:raise RuntimeError(f'Missing/empty required file: {path}')

def opener(p,mode='rt'):
    return (gzip.open(p,mode,encoding='utf-8-sig',newline='') if str(p).endswith('.gz') else
            open(p,mode,encoding='utf-8-sig',newline=''))

def tsv_rows(path):
    with opener(path) as f:
        rd=csv.DictReader(f,delimiter='\t')
        if not rd.fieldnames:raise RuntimeError(f'No TSV header: {path}')
        for i,r in enumerate(rd,2):
            if None in r:raise RuntimeError(f'Invalid TSV fields: {path}:{i}')
            yield r

def tsv_header(path):
    with opener(path) as f:return next(csv.reader(f,delimiter='\t'))

def write_tsv(path, fields, rows):
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n');w.writeheader()
        w.writerows(rows)
    log(f'Wrote {path.name}')

def normid(k):
    x=str(k).strip();a=x.split('-',3)
    if len(a)!=4:raise RuntimeError(f'Unparseable Variant_ID {x!r}')
    chrom=a[0].upper().removeprefix('CHR')
    return f'{chrom}-{int(a[1])}-{a[2].upper()}-{a[3].upper()}'

def vkey(chrom,pos,ref,alt):
    return normid(f'{chrom}-{pos}-{ref}-{alt}')

def getid(r):
    k=r.get('Variant_ID','').strip()
    if k:return normid(k)
    return vkey(r.get('Chr') or r.get('CHROM'),r.get('Pos') or r.get('POS'),r.get('Ref') or r.get('REF'),r.get('Alt') or r.get('ALT'))

def pick(r,*names):
    for n in names:
        if n in r and str(r[n]).strip().lower() not in MISSING:return str(r[n]).strip()
    return ''

def field_present(r,*names):return any(n in r for n in names)

def af_val(v,where):
    if str(v if v is not None else '').strip().lower() in MISSING:return None
    try:
        x=float(v)
        if 0<=x<=1 and math.isfinite(x):return x
    except (TypeError,ValueError):pass
    raise RuntimeError(f'Invalid allele frequency {v!r} at {where}')

def outaf(v):return '' if v is None else f'{v:.9g}'

def nint(v,k):
    try:return int(v)
    except (ValueError,TypeError):raise RuntimeError(f'Invalid integer {k}: {v!r}')

def import_original(p):
    spec=importlib.util.spec_from_file_location('original_jointcall_02',p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    for name in ('parse_gt','load_phenotypes','load_sample_warnings','functional_reason'):
        if not hasattr(m,name):raise RuntimeError(f'Original 02 module lacks {name}')
    return m

def genes_from(raw):
    if not raw or raw.lower() in MISSING:return []
    return sorted({s.strip() for s in re.split(r'[;,|]',raw) if s.strip() and s.strip().lower() not in MISSING})

def clinical(r):
    # Source column presence is different from a genuinely missing annotation.
    labels=('CLNSIG','ClinVar_CLNSIG','ClinVar_significance','ClinVar_Significance',
            'Clinvar_Significance','ClinVar_ClinicalSignificance','ClinicalSignificance')
    records=('ClinVar_record','ClinVar_Record','ClinVar_record_ID','ClinVar_ID','ClinVar_accession','Clinvar_record')
    sig=pick(r,*labels);record=pick(r,*records)
    if sig:
        return sig+(f' [{record}]' if record and record not in sig else ''),'present'
    if record:return f'Record: {record} (significance unavailable)','record_only'
    if field_present(r,*labels,*records):return 'Not annotated','absent'
    return 'Not supplied','not_supplied'

def clinical_bucket(s):
    x=s.lower()
    if 'pathogenic' in x and not x.startswith('uncertain'):
        # Conflicting interpretations retained separately.
        if 'conflict' in x:return 'Conflicting'
        if 'benign' in x:return 'Conflicting'
        return 'P/LP'
    if 'uncertain' in x or re.search(r'\bvus\b',x):return 'VUS'
    if 'benign' in x:return 'Benign/LB'
    if s in ('Not annotated','Not supplied'):return s
    return 'Other'

def first_annotation_field(r,*names):
    for name in names:
        if name in r:return r[name],name
    return '',''

def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--merged',type=Path,default=DEF_MERGED)
    ap.add_argument('--annotations',type=Path,default=DEF_ANNOT)
    ap.add_argument('--vcf',type=Path,default=DEF_VCF)
    ap.add_argument('--original-qc-code',type=Path,default=DEF_ORIG)
    ap.add_argument('--sample-qc',type=Path,default=DEF_QC)
    ap.add_argument('--unrelated',type=Path,default=DEF_UNREL)
    ap.add_argument('--phenotype',type=Path,default=DEF_PHENO)
    ap.add_argument('--gnomad',type=Path,default=DEF_GNOMAD)
    ap.add_argument('--family-events',type=Path,default=DEF_EVENTS)
    ap.add_argument('--refreshed-gnomad',type=Path,help='Optional independently re-queried prior variants, same schema as gnomad_v4_results.tsv')
    ap.add_argument('--local-flags',type=Path,help='Output of gnomad_flag_evidence_from_vcf.py based on lab-local echtvar/vcfanno reference mappings')
    ap.add_argument('--allow-no-local-flags',action='store_true',help='Emit provisional lists without required local flag verification; default keeps otherwise passing variants in REVIEW')
    ap.add_argument('--outdir',type=Path,default=SUB/'final_candidate_analysis_20261110')
    ap.add_argument('--internal-af-cutoff',type=float,default=.015)
    ap.add_argument('--external-af-cutoff',type=float,default=.001)
    ap.add_argument('--expect-variants',type=int,default=14124)
    ap.add_argument('--expect-samples',type=int,default=1245)
    ap.add_argument('--expect-families',type=int,default=659)
    ap.add_argument('--expect-unrelated',type=int,default=953)
    ap.add_argument('--expect-masked',type=int,default=9)
    a=ap.parse_args()
    for p in [a.merged,a.annotations,a.vcf,a.original_qc_code,a.sample_qc,a.unrelated,a.phenotype,a.gnomad,a.family_events]:require(p)
    if a.refreshed_gnomad:require(a.refreshed_gnomad)
    if a.local_flags:require(a.local_flags)
    if a.outdir.exists() and any(a.outdir.iterdir()):raise RuntimeError(f'Refusing to overwrite nonempty output directory {a.outdir}')
    # Require exact allele annotation coverage before creating outputs.
    orig=import_original(a.original_qc_code)
    merged={}
    for r in tsv_rows(a.merged):
        k=getid(r)
        if k in merged:raise RuntimeError('Duplicate merged allele '+k)
        merged[k]=r
    if len(merged)!=a.expect_variants:raise RuntimeError(f'Expected {a.expect_variants} variants; found {len(merged)}')
    log(f'Loaded {len(merged):,} merged alleles')
    ann={}
    for r in tsv_rows(a.annotations):
        k=getid(r)
        if k not in merged:continue
        if k in ann:raise RuntimeError('Duplicate annotated allele '+k)
        gene=pick(r,'Gene_refGene','Gene.refGene','Gene')
        func=pick(r,'Func_refGene','Func.refGene')
        exon=pick(r,'ExonicFunc_refGene','ExonicFunc.refGene')
        consequence=orig.functional_reason(func,exon)
        if not gene or not consequence:raise RuntimeError(f'Missing gene or non-LoF annotation at {k}: gene={gene!r} function={func!r}, exonic={exon!r}')
        ann[k]={**r,'_gene':gene,'_func':func,'_exon':exon,'_consequence':consequence}
    if len(ann)!=len(merged):
        missing=list(merged.keys()-ann.keys())[:8]
        raise RuntimeError(f'Annotation matching incomplete {len(ann)}/{len(merged)}. First missing: {missing}')
    log('Exact ANNOVAR allele matching and LoF/splicing function checks passed')
    # gnomAD new query data overrides and verifies merged decision status.
    query={}
    for p in [a.gnomad]+([a.refreshed_gnomad] if a.refreshed_gnomad else []):
        for r in tsv_rows(p):
            k=getid(r)
            if k not in merged:raise RuntimeError('gnomAD result allele missing from merged file '+k)
            if k in query and p == a.gnomad:raise RuntimeError('Duplicate gnomAD result '+k)
            query[k]=r
    log(f'Loaded {len(query):,} available gnomAD query results')
    # Additional local gnomAD FILTER evidence from echtvar/vcfanno. We do not infer
    # source FILTER semantics; this file must have been generated from exact field
    # names verified against lab documentation. API/local conflicts => REVIEW.
    local={}
    if a.local_flags:
        for r in tsv_rows(a.local_flags):
            k=getid(r)
            if k not in merged:raise RuntimeError('Local flag allele absent in merged file '+k)
            if k in local:raise RuntimeError('Duplicate local flag allele '+k)
            status=pick(r,'Flag_Assessment').upper()
            if status not in {'PASS','FLAGGED','UNRESOLVED','EXCLUDE_AF'}:raise RuntimeError('Unknown local flag status '+k+': '+status)
            if not r.get('Flag_Source','').strip():raise RuntimeError('Missing provenance for local gnomAD flag record '+k)
            for required_col in ('Local_gnomAD_AF','Local_gnomAD_AN','ClinVar_Significance','ClinVar_ID'):
                if required_col not in r:raise RuntimeError('Local annotation column missing: '+required_col+'; use CONFIGURED extractor, not the old flag-only script')
            local[k]=r
        if len(local)!=len(merged):raise RuntimeError(f'Incomplete local flag file: {len(local)}/{len(merged)}; do not treat missing allele as PASS')
        log('Loaded complete local gnomAD flag cross-check table')
    elif not a.allow_no_local_flags:
        log('Local echtvar/vcfanno flag evidence absent: otherwise passing variants will remain REVIEW')
    expected_new={'RARE_RETAIN':'NEW_GNOMAD_RARE_PASS','EXCLUDE_AF':'NEW_GNOMAD_EXCLUDE_AF',
                  'EXCLUDE_QUALITY_FLAGS':'NEW_GNOMAD_EXCLUDE_QUALITY',
                  'NOT_FOUND_UNRESOLVED':'NEW_GNOMAD_NOT_FOUND_REVIEW'}
    for r in tsv_rows(a.gnomad):
        k=getid(r);previous=merged[k].get('gnomAD_final_decision','')
        if previous!=expected_new[r['Gnomad_Decision']]:
            raise RuntimeError(f'Merged gnomAD decision differs from original query at {k}: {previous} vs {r["Gnomad_Decision"]}')
    known_decisions={'RARE_RETAIN','EXCLUDE_AF','EXCLUDE_QUALITY_FLAGS','NOT_FOUND_UNRESOLVED'}
    for k,r in query.items():
        d=r.get('Gnomad_Decision','')
        if d not in known_decisions:raise RuntimeError('Unknown gnomAD decision '+k+': '+d)
        v=af_val(r.get('v4_max_AF_any'),k+' gnomAD AF')
        if d=='RARE_RETAIN':
            if v is None or v>=a.external_af_cutoff:raise RuntimeError(f'Invalid rare pass: {k} AF={v}')
            if r.get('Any_gnomAD_callset_filtered','').strip()=='1':raise RuntimeError('gnomAD rare pass has quality filter '+k)
        elif d=='EXCLUDE_AF' and v is not None and v<a.external_af_cutoff:
            raise RuntimeError(f'AF exclusion contradictory to verified max AF: {k} AF={v}')
    # Read original phenotype/QC only. Validate expected sample membership first from VCF header.
    contamination=orig.load_sample_warnings(a.sample_qc)
    with open(a.unrelated) as f:
        unrelated={line.strip().split()[0] for line in f if line.strip() and not line.startswith('#')}
    if len(contamination)!=a.expect_masked or len(unrelated)!=a.expect_unrelated:
        raise RuntimeError(f'Unexpected flagged/unrelated set sizes {len(contamination)}/{len(unrelated)}')
    # Family event lookup: all old per-variant/per-family rows must be corroborated.
    event={}
    event_by_variant=collections.defaultdict(set)
    event_class_by_variant=collections.defaultdict(collections.Counter)
    for r in tsv_rows(a.family_events):
        k=getid(r)
        if k not in merged:raise RuntimeError('Family event allele absent from merged: '+k)
        fam=r['Family_ID']
        key=(k,fam)
        if key in event:raise RuntimeError('Duplicate family event '+str(key))
        event[key]=tuple(nint(r[f],key) for f in ['Affected_ALT','Unaffected_ALT','Unknown_ALT','Affected_REF','Unaffected_REF','Affected_NO_CALL','Unaffected_NO_CALL'])
        event_by_variant[k].add(key)
        event_class_by_variant[k][r.get('Family_status','')]+=1
    log(f'Loaded {len(event):,} prior family event rows for count validation')
    # Confirm previously reported per-variant segregation class totals agree with event table.
    statuses=[('Families_strict_compatible','COMPATIBLE_STRICT_SCREEN'),
              ('Families_affected_carrier_only','AFFECTED_CARRIER_ONLY'),
              ('Families_inconsistent','INCONSISTENT_WITH_STRICT_RULE'),
              ('Families_incomplete','INCOMPLETE_GENOTYPES_REVIEW')]
    for k,r in merged.items():
        for col,prior_status in statuses:
            if col in r and nint(r[col],k+' '+col)!=event_class_by_variant[k][prior_status]:
                raise RuntimeError(f'Segregation count mismatch {k} {col}: merged={r[col]}, events={event_class_by_variant[k][prior_status]}')
    log('Segregation class counts agree with original event classifications')
    # The VCF is already filtered and split into one alternate allele per row.
    stats={};site_seen=set();sample_ids=[];phenos={};family_by_sample={}
    family_data=collections.defaultdict(lambda:collections.Counter())
    phenocount=collections.Counter()
    qc_mismatches=[];event_mismatches=[]
    with opener(a.vcf) as f:
        for line in f:
            if line.startswith('##'):continue
            if line.startswith('#CHROM'):
                if sample_ids:raise RuntimeError('Duplicate VCF header')
                sample_ids=line.rstrip('\n').split('\t')[9:]
                if len(sample_ids)!=a.expect_samples or len(set(sample_ids))!=len(sample_ids):raise RuntimeError(f'VCF has {len(sample_ids)} nonunique/wrong number of samples')
                if not contamination<=set(sample_ids) or not unrelated<=set(sample_ids):raise RuntimeError('QC or unrelated set not a subset of VCF')
                phenos=orig.load_phenotypes(a.phenotype,set(sample_ids))
                for sid in sample_ids:
                    family_by_sample[sid]=sid.rsplit('-',1)[0]
                    td=phenos[sid]['td']
                    if td not in {'affected','unaffected','unknown'}:raise RuntimeError('Unknown TD label for '+sid)
                    phenocount[td]+=1;family_data[family_by_sample[sid]][td]+=1
                if len(family_data)!=a.expect_families:
                    raise RuntimeError(f'Unexpected family count {len(family_data)} expected {a.expect_families}')
                log(f'VCF phenotype check passed: {len(sample_ids)} samples, {len(family_data)} families, {dict(phenocount)}')
                continue
            if line.startswith('#'):continue
            if not sample_ids:raise RuntimeError('No VCF header found')
            p=line.rstrip('\n').split('\t')
            if len(p)!=9+len(sample_ids):raise RuntimeError('VCF sample column mismatch')
            if ',' in p[4]:raise RuntimeError('Expected split ALT VCF; found multiple ALTs at '+p[0]+':'+p[1])
            k=vkey(p[0],p[1],p[3],p[4])
            if k not in merged:raise RuntimeError('VCF variant not present in annotated 14,124 table: '+k)
            if k in site_seen:raise RuntimeError('Duplicate VCF allele '+k)
            site_seen.add(k)
            if p[6]!='PASS':raise RuntimeError('Filtered VCF contains FILTER !=PASS at '+k+': '+p[6])
            fmt=p[8].split(':')
            counters={x:[0,0] for x in ['ALL','UNRELATED','UNAFFECTED','UNAFFECTED_UNRELATED']}
            counts=collections.Counter();families=collections.defaultdict(lambda:collections.Counter())
            a_carrier=set();u_carrier=set();affected_families=set()
            for i,sid in enumerate(sample_ids):
                td=phenos[sid]['td'];family=family_by_sample[sid]
                if sid in contamination:status='NO_CALL';dosage=0;ploidy=0
                else:
                    status,_=orig.parse_gt(fmt,p[i+9],10,20,.05,.25,.75,4)
                    dosage=0;ploidy=0
                    if status!='NO_CALL':
                        gt=p[i+9].split(':')[fmt.index('GT')]
                        alleles=re.split('[/|]',gt)
                        if not alleles or any(x not in ('0','1') for x in alleles):raise RuntimeError('Unexpected GT on passing sample '+k)
                        ploidy=len(alleles);dosage=alleles.count('1')
                        if (dosage>0)!=(status=='ALT'):raise RuntimeError('parse_gt ALT disagreement at '+k+' / '+sid)
                counts[(td,status)]+=1
                if status=='ALT':
                    families[family][f'{td}_ALT']+=1
                    if td=='affected':a_carrier.add(sid);affected_families.add(family)
                    elif td=='unaffected':u_carrier.add(sid)
                elif status=='REF' and td!='unknown':families[family][f'{td}_REF']+=1
                elif status=='NO_CALL' and td!='unknown':families[family][f'{td}_NO_CALL']+=1
                if status!='NO_CALL':
                    for group,include in [('ALL',True),('UNRELATED',sid in unrelated),('UNAFFECTED',td=='unaffected'),('UNAFFECTED_UNRELATED',td=='unaffected' and sid in unrelated)]:
                        if include:counters[group][0]+=dosage;counters[group][1]+=ploidy
            for td,prefix in [('affected','Affected'),('unaffected','Unaffected'),('unknown','Unknown')]:
                for status,label in [('ALT','QC_ALT'),('REF','REF'),('NO_CALL','NO_CALL')]:
                    col=f'{prefix}_{label}_samples'
                    if col not in merged[k]:continue
                    before=nint(merged[k][col],k+' '+col)
                    after=counts[(td,status)]
                    if before!=after:qc_mismatches.append((k,col,before,after))
            # Validate every historical event and verify no new event appears.
            current={}
            for fam,c in families.items():
                if c['affected_ALT'] or c['unaffected_ALT'] or c['unknown_ALT']:
                    tup=(c['affected_ALT'],c['unaffected_ALT'],c['unknown_ALT'],c['affected_REF'],c['unaffected_REF'],c['affected_NO_CALL'],c['unaffected_NO_CALL'])
                    current[(k,fam)]=tup
                    if event.get((k,fam))!=tup:event_mismatches.append(((k,fam),event.get((k,fam)),tup))
            # Important: the old events can include a family with no ALT? If so stop on mismatch.
            event_rows=event_by_variant.get(k,set())
            if set(current)!=event_rows:
                event_mismatches.append((k,'event key difference',len(set(current)^event_rows)))
            stats[k]={'group':counters,'affected_carriers':a_carrier,'unaffected_carriers':u_carrier,
                      'affected_families':affected_families,'counts':counts}
            if len(site_seen)%2000==0:log(f'Checked {len(site_seen):,} of {len(merged):,} variant records')
    if site_seen!=set(merged):raise RuntimeError(f'VCF allele mismatch, found {len(site_seen)} expected {len(merged)}; first unmatched {list(set(merged)-site_seen)[:3]}')
    if qc_mismatches or event_mismatches:
        a.outdir.mkdir(parents=True,exist_ok=True)
        write_tsv(a.outdir/'STOP_COUNT_MISMATCH.tsv',['Variant_ID','Field','Prior','Current'],
                  [{'Variant_ID':k,'Field':c,'Prior':old,'Current':new} for k,c,old,new in qc_mismatches[:1000]])
        write_tsv(a.outdir/'STOP_EVENT_MISMATCH.tsv',['Issue'],[{'Issue':str(x)} for x in event_mismatches[:1000]])
        raise RuntimeError(f'COUNT VALIDATION FAILED: {len(qc_mismatches)} genotype-count and {len(event_mismatches)} family-event discrepancies. No final lists produced.')
    log('All variant genotype counts and family carrier events reproduced')
    # Derive exact per-allele evidence/decisions; historical gnomAD passes are never fully VERIFIED.
    result=[];full=[];audit_count=collections.Counter();review_reason_count=collections.Counter();eligible_ids=set()
    allele_gene_map={}
    for k,r in merged.items():
        an=ann[k];st=stats[k];cv,cv_status=clinical(an)
        # vcfanno exact REF/ALT ClinVar is authoritative over ANNOVAR/VEP overlaps.
        # The local CSV is generated even when no ClinVar entry matched.
        lr=local.get(k)
        if lr and 'ClinVar_Significance' in lr:
            sig=pick(lr,'ClinVar_Significance')
            cv_id=pick(lr,'ClinVar_ID')
            cv=(sig+(f' [{cv_id}]' if cv_id else '')) if sig else ('Record: '+cv_id if cv_id else 'Not annotated')
            cv_status='exact_allele_vcfanno' if sig or cv_id else 'absent_exact_allele_vcfanno'
        gene=an['_gene'];consequence=an['_consequence']
        afs={g:(ac/an_ if an_ else None) for g,(ac,an_) in st['group'].items()}
        q=query.get(k)
        previous=r.get('gnomAD_final_decision','')
        has_prior=previous=='PRIOR_RARE_PASS_NOT_REQUERIED'
        local_row=local.get(k)
        local_status=local_row.get('Flag_Assessment','') if local_row else 'NOT_RUN'
        local_detail=local_row.get('Flag_Detail','') if local_row else ''
        local_af=af_val(local_row.get('Local_gnomAD_AF',''),k+' local gnomAD AF') if local_row else None
        local_good=local_status=='PASS' and local_af is not None and local_af<a.external_af_cutoff
        local_discordant=local_status in ('FLAGGED','EXCLUDE_AF') and (q and q.get('Gnomad_Decision')=='RARE_RETAIN' or has_prior)

        if q:
            qdec=q['Gnomad_Decision'];gnaf=af_val(q.get('v4_max_AF_any'),k)
            qc={'RARE_RETAIN':'PASS','EXCLUDE_AF':'FAIL_AF','EXCLUDE_QUALITY_FLAGS':'FAIL_QUALITY','NOT_FOUND_UNRESOLVED':'UNRESOLVED'}[qdec]
        else:
            qdec='';gnaf=af_val(r.get('gnomAD_max_AF'),k)
            qc='PRIOR_REVIEW' if has_prior else ('PRIOR_EXCLUSION' if previous=='PRIOR_EXCLUDE_NOT_REQUERIED' else 'UNRESOLVED')
        # Keep the more conservative documented AF when v4.0 and local v4.1 both report values.
        if local_af is not None and (gnaf is None or local_af>gnaf):gnaf=local_af
        old_nn,old_nn_field=first_annotation_field(an,'non_neuro_AF_popmax','Old_non_neuro_AF_popmax','gnomAD_non_neuro_AF_popmax')
        old_all,old_all_field=first_annotation_field(an,'AF_popmax','Old_AF_popmax')
        old1=af_val(old_nn,k+' old non-neuro AF');old2=af_val(old_all,k+' old AF popmax')
        oldfields_present=bool(old_nn_field and old_all_field)
        # Preserve comparison where input contains historical frequency fields.
        reasons=[]
        if st['counts'][('affected','ALT')]==0:reasons.append('No affected carrier')
        if afs['ALL'] is None:reasons.append('No callable internal AF')
        elif afs['ALL']>a.internal_af_cutoff:reasons.append('Internal AF > 0.015')
        if old1 is not None and old1>=a.external_af_cutoff:reasons.append('ANNOVAR non_neuro_AF_popmax >= 0.001')
        if old2 is not None and old2>=a.external_af_cutoff:reasons.append('ANNOVAR AF_popmax >= 0.001')
        if qdec=='EXCLUDE_AF':reasons.append('gnomAD AF >= 0.001')
        elif qdec=='EXCLUDE_QUALITY_FLAGS':reasons.append('gnomAD quality flags')
        elif previous=='PRIOR_EXCLUDE_NOT_REQUERIED' and not q:reasons.append('Previously excluded gnomAD record (not revalidated)')
        # Even previously retained AF must be coherent with the written cutoff.
        if has_prior and not q and not local_good and (gnaf is None or gnaf >= a.external_af_cutoff):
            reasons.append('Prior gnomAD rarity record lacks valid AF < 0.001')
        positive_gnomad=(qdec=='RARE_RETAIN') or (not q and has_prior) or (local_good and qdec not in ('EXCLUDE_AF','EXCLUDE_QUALITY_FLAGS') and previous!='PRIOR_EXCLUDE_NOT_REQUERIED')
        if reasons:
            disposition='EXCLUDE'
        else:
            doubts=[]
            if not q:
                if has_prior and not local_good:doubts.append('Prior gnomAD pass not independently rechecked')
                elif not local_good:doubts.append('gnomAD unresolved or not queried')
            elif qdec!='RARE_RETAIN':doubts.append('gnomAD not-found or unresolved in API v4.0; local v4.1 may differ')
            if local_status in ('FLAGGED','EXCLUDE_AF'):
                doubts.append('Local gnomAD v4.1 frequency/FILTER conflicts with API v4.0 or previous decision')
            elif local_status=='UNRESOLVED':
                doubts.append('Local gnomAD FILTER annotation unresolved')
            elif local_status=='NOT_RUN' and not a.allow_no_local_flags:
                doubts.append('PI-requested echtvar/vcfanno gnomAD flags not yet validated')
            if not oldfields_present:doubts.append('Original ANNOVAR frequency columns not supplied')
            if cv_status=='not_supplied':doubts.append('ClinVar columns not supplied')
            if r.get('Chromosome_scope')!='AUTOSOMAL':doubts.append('Sex/other chromosome requires ploidy/QC review')
            disposition='REVIEW' if doubts else 'PASS'
            reasons=doubts or ['All checked filters passed; segregation is an evidence tier']
        # A local FLAGGED record is not silently equated with an API FILTER failure:
        # discordant site/source mappings need human review before exclusion.
        if local_status in ('FLAGGED','EXCLUDE_AF') and qc=='PASS':qc='REVIEW_FLAG_CONFLICT'
        elif local_status in ('FLAGGED','EXCLUDE_AF') and qc=='PRIOR_REVIEW':qc='REVIEW_LOCAL_FLAG'
        if positive_gnomad and disposition in ('PASS','REVIEW'):
            eligible_ids.add(k)
        bucket=r.get('Final_review_bucket','')
        family_label={'CANDIDATE_STRICT_COMPATIBLE':'Strict-compatible','CANDIDATE_AFFECTED_CARRIER_ONLY':'Affected-carrier-only','CANDIDATE_FAMILY_REVIEW':'Family review'}.get(bucket,'No affected-carrier candidate')
        # For every 14k row, do not apply the family tier as an exclusion beyond NO AFFECTED.
        summary={
            'Variant_ID':k,'Gene':gene,'Function':consequence,'gnomAD_AF':outaf(gnaf),'gnomAD_QC':qc,
            'Internal_AF':outaf(afs['ALL']),'Affected_Carriers':len(st['affected_carriers']),
            'Unaffected_Carriers':len(st['unaffected_carriers']),'Affected_Families':len(st['affected_families']),
            'Strict_Families':nint(r['Families_strict_compatible'],k+' strict families'),
            'Inconsistent_Families':nint(r['Families_inconsistent'],k+' inconsistent families'),
            'Segregation':family_label,'ClinVar':cv,'Filter_Status':disposition,'Filter_Reason':'; '.join(reasons)
        }
        if len(st['affected_families'])!=nint(r['Families_with_affected_ALT'],k+' families affected'):
            raise RuntimeError('Family count mismatch '+k)
        audit={**summary,'Internal_ALL_AC':st['group']['ALL'][0],'Internal_ALL_AN':st['group']['ALL'][1],
               'Internal_Unrelated_AF':outaf(afs['UNRELATED']),
               'Internal_Unaffected_AF':outaf(afs['UNAFFECTED']),
               'Internal_Unaffected_Unrelated_AF':outaf(afs['UNAFFECTED_UNRELATED']),
               'Old_non_neuro_AF_popmax':outaf(old1),'Old_AF_popmax':outaf(old2),
               'Old_non_neuro_source':old_nn_field or 'NOT_SUPPLIED',
               'Old_AF_popmax_source':old_all_field or 'NOT_SUPPLIED',
               'Gnomad_Decision':qdec or previous,'Gnomad_Source':'QUERY_RESULT' if q else 'PREVIOUS_OR_MISSING',
               'Local_Flag_Status':local_status,'Local_Flag_Details':local_detail,
               'Local_Flag_Source':(local_row.get('Flag_Source','') if local_row else 'NOT_RUN'),
               'ClinVar_Source':cv_status,
               'ClinVar_Disease':pick(local_row,'ClinVar_Disease') if local_row else pick(an,'ClinVar_CLNDN','CLNDN','ClinVar_disease'),
               'ClinVar_Review_Status':pick(local_row,'ClinVar_Review') if local_row else pick(an,'ClinVar_CLNREVSTAT','CLNREVSTAT'),
               'ClinVar_Stars':pick(local_row,'ClinVar_Stars') if local_row else '',
               'Local_AF_v4.1':outaf(local_af),'Local_LCR':local_row.get('Local_LCR','') if local_row else '',
               'Local_SegDup':local_row.get('Local_SegDup','') if local_row else '',
               'Chromosome_scope':r.get('Chromosome_scope',''),
               'Original_Family_Evidence':r.get('Family_evidence',''),
               'Eligible_for_Prioritized_Candidate_Table':int(k in eligible_ids)}
        result.append(summary);full.append(audit);audit_count[disposition]+=1
        if disposition=='REVIEW':
            for d in reasons:review_reason_count[d]+=1
        allele_gene_map[k]=genes_from(gene)
    # Parent/offspring inheritance is not checked by this step; it preserves validated prior segregation category.
    candidates=[r for r in result if r['Variant_ID'] in eligible_ids]
    passing=[r for r in candidates if r['Filter_Status']=='PASS']
    reviewing=[r for r in candidates if r['Filter_Status']=='REVIEW']
    population_unresolved=[r for r in result if r['Filter_Status']=='REVIEW' and r['Variant_ID'] not in eligible_ids]
    def build_genes(allelelist):
        ix=collections.defaultdict(lambda:{'variant_ids':set(),'affected':set(),'unaffected':set(),'families':set(),
                                            'types':set(),'strict':0,'af':[],'cv':collections.Counter()})
        for r in allelelist:
            k=r['Variant_ID'];st=stats[k]
            for g in allele_gene_map[k]:
                z=ix[g];z['variant_ids'].add(k);z['affected'].update(st['affected_carriers']);z['unaffected'].update(st['unaffected_carriers']);z['families'].update(st['affected_families'])
                z['types'].update(r['Function'].split(';'))
                z['strict']+=int(r['Strict_Families']>0)
                if r['gnomAD_AF']!='':z['af'].append(float(r['gnomAD_AF']))
                z['cv'][clinical_bucket(r['ClinVar'])]+=1
        rows=[]
        for g,z in ix.items():
            c=z['cv'];cs='; '.join(f'{label}: {c[label]}' for label in ['P/LP','Conflicting','VUS','Benign/LB','Other','Not annotated','Not supplied'] if c[label]) or 'Not annotated'
            rows.append({'Gene':g,'Candidate_Variants':len(z['variant_ids']),'Variant_Types':'; '.join(sorted(z['types'])),
                         'Affected_Carriers':len(z['affected']),'Unaffected_Carriers':len(z['unaffected']),
                         'Affected_Families':len(z['families']),'Strict_Segregation_Variants':z['strict'],
                         'Lowest_gnomAD_AF':outaf(min(z['af']) if z['af'] else None),'ClinVar_Summary':cs})
        return sorted(rows,key=lambda x:(-int(x['Strict_Segregation_Variants']),-int(x['Candidate_Variants']),x['Gene']))
    genes_all=build_genes(candidates);genes_pass=build_genes(passing);genes_review=build_genes(reviewing)
    tiers=collections.Counter(r['Segregation'] for r in candidates)
    sizes=collections.Counter('1' if sum(c.values())==1 else '2' if sum(c.values())==2 else '3' if sum(c.values())==3 else '4+' for c in family_data.values())
    a.outdir.mkdir(parents=True,exist_ok=True)
    write_tsv(a.outdir/'jointcall_candidate_variants_for_dr_xing.tsv',ALLELES,candidates)
    write_tsv(a.outdir/'jointcall_candidate_genes_for_dr_xing.tsv',GENES,genes_all)
    write_tsv(a.outdir/'jointcall_candidate_variants_final_verified.tsv',ALLELES,passing)
    write_tsv(a.outdir/'jointcall_candidate_genes_final_verified.tsv',GENES,genes_pass)
    # Matching the agreed 9-column gene and 15-column variant final-filter file names.
    # These contain only PASS entries, not REVIEW entries.
    write_tsv(a.outdir/'jointcall_candidate_variants_final.tsv',ALLELES,passing)
    write_tsv(a.outdir/'jointcall_candidate_genes_final.tsv',GENES,genes_pass)
    write_tsv(a.outdir/'jointcall_candidate_variants_need_review.tsv',ALLELES,reviewing)
    write_tsv(a.outdir/'jointcall_candidate_genes_need_review.tsv',GENES,genes_review)
    write_tsv(a.outdir/'population_unresolved_not_yet_candidates.tsv',ALLELES,population_unresolved)
    write_tsv(a.outdir/'ALL_14124_COMPLETE_AUDIT.tsv',list(full[0].keys()),full)
    write_tsv(a.outdir/'family_summary.tsv',['Family_ID','Sequenced_Individuals','Affected','Unaffected','Unknown'],[
        {'Family_ID':fam,'Sequenced_Individuals':sum(c.values()),'Affected':c['affected'],'Unaffected':c['unaffected'],'Unknown':c['unknown']}
        for fam,c in sorted(family_data.items())])
    # Generate concise methods/results and email using this run's real counts only.
    n_prior=sum(1 for k,r in merged.items() if r.get('gnomAD_final_decision')=='PRIOR_RARE_PASS_NOT_REQUERIED' and k not in query)
    newcounts=collections.Counter(x['Gnomad_Decision'] for x in tsv_rows(a.gnomad))
    date=datetime.now().strftime('%B %d, %Y')
    summary=textwrap.dedent(f'''\
    # TIC Genetics — joint-called Batches 11–14: methods and candidate summary
    Generated: {date}

    ## Dataset
    Sequenced individuals: {len(sample_ids):,} (TD affected {phenocount['affected']:,}; unaffected {phenocount['unaffected']:,}; unknown {phenocount['unknown']:,}).
    Families represented: {len(family_data):,}.
    Family size breakdown (1 / 2 / 3 / >=4 sequenced): {dict(sizes)}.
    These sizes are not identical to the historical clinical family-type labels.
    Contamination-proxy samples masked in genotype evidence: {len(contamination)}.
    Kinship-unrelated reference samples: {len(unrelated)}.
    Site-filtered split-ALT LoF/coding-splicing variants: {len(merged):,}.

    ## Analysis and filters
    1. Source is prefiltered, split-ALT joint-called GRCh38 VCF, checked FILTER=PASS.
       Yale's original QC and extra VQSR/capture filters were applied upstream; this stage verifies
       the saved passing site set, not the entire upstream raw calling history.
    2. Matched exact chromosome-position-REF-ALT keys to {len(ann):,} ANNOVAR records;
       verified retained stopgain, stoploss, frameshift or coding-splicing consequences.
       This analysis does not include the separately studied nonsynonymous category.
    3. Reused original genotype QC: DP>=10, GQ>=20, ALT reads>=4, heterozygous
       ALT balance 0.25–0.75, reference ALT fraction<=0.05, homozygous ALT balance>=0.8.
       Contamination-proxy calls masked. Reconciled every saved affected/unaffected/unknown
       carrier count and each original family event against genotype data.
    4. Internal AF calculated as ALT copies / confidently called allele copies at each site;
       excluded when full joint-call AF > {a.internal_af_cutoff}.
       Also reported unrelated, unaffected, and unaffected-unrelated AF in complete audit.
       The original 0.015 threshold was exploratory; unaffected-only primary filtering
       was not applied.
    5. External gnomAD v4 max population/callset AF threshold < {a.external_af_cutoff};
       filter low-quality gnomAD sites; checked both legacy ANNOVAR popmax fields when supplied.
       Missing fields and gnomAD not-found cases are not interpreted as zero frequency.
       Previously retained population calls not independently re-queried remain REVIEW.
       Echtvar/vcfanno local FILTER-source cross-check: {'COMPLETE' if local else ('NOT RUN; all otherwise passing records held in REVIEW' if not a.allow_no_local_flags else 'NOT RUN; user chose provisional mode')}.
       Any conflict between local flag evidence and API results remains REVIEW pending a
       source-field and allele-coordinate audit rather than being silently excluded.
    6. Family segregation is an evidence tier, not a mandatory genetic model; strict-compatible,
       affected-carrier-only, and conflicting/incomplete are retained separately.
       No required gene recurrence, NDD overlap, de novo finding, or ClinVar P/LP label.
       ClinVar records are annotations, not proof of Tourette association.

    ## Results (variant level)
    New gnomAD v4 query outcomes from source results: {dict(sorted(newcounts.items()))}
    Prior rare-retained decisions still un-requeried: {n_prior:,}.
    PASS all currently checkable criteria with needed source provenance: {len(passing):,} alleles / {len(genes_pass):,} genes.
    REVIEW requiring data/provenance or scientific confirmation: {len(reviewing):,} alleles / {len(genes_review):,} genes.
    EXCLUDE: {audit_count['EXCLUDE']:,} alleles.
    Population unresolved without positive rare screening: {len(population_unresolved):,} alleles (separate, NOT counted as prioritized candidates).
    Main prioritization list (PASS+REVIEW with positive rarity evidence): {len(candidates):,} alleles / {len(genes_all):,} genes.
    Candidate segregation evidence distribution: {dict(tiers)}.
    Local gnomAD flag status distribution: {dict(collections.Counter(x['Flag_Assessment'] for x in local.values())) if local else 'NOT AVAILABLE'}.
    Most common review reasons: {dict(review_reason_count.most_common(12))}

    ## Interpretation and limitations
    Counts above are variant-level or unique-person gene-level as labeled. Gene carrier and
    affected-family counts are DISTINCT across all selected alleles in that gene.
    Candidate status does not establish pathogenicity, penetrance, or disease association.
    Exact historical ANNOVAR AF and ClinVar field completeness determines REVIEW status.
    Sex/other chromosome alleles remain REVIEW pending ploidy/QC interpretation.
    No external database calls were performed by this finalizer.
    This package does not independently rerun Yale's VQSR, annotation, or full VCF creation.
    Review prior gnomAD decisions and confirm internal AF policy before calling
    the complete PASS+REVIEW table a fully final disease-gene list.
    ''')
    (a.outdir/'METHODS_AND_FAMILY_RESULTS.md').write_text(summary,encoding='utf-8')
    manifest={'date_utc':datetime.now(timezone.utc).isoformat(),
              'inputs':{name:str(getattr(a,name)) for name in ['merged','annotations','vcf','sample_qc','unrelated','phenotype','original_qc_code','gnomad','family_events']},
              'local_flags_input':str(a.local_flags) if a.local_flags else None,
              'local_flags_required_for_PASS':not a.allow_no_local_flags,
              'thresholds':{'internal_full_af_exclude_gt':a.internal_af_cutoff,'population_af_pass_lt':a.external_af_cutoff},
              'counts':{'variants':len(merged),'families':len(family_data),'sample_total':len(sample_ids),'sample_affected':phenocount['affected'],
                        'sample_unaffected':phenocount['unaffected'],'sample_unknown':phenocount['unknown'],
                        'PASS':len(passing),'REVIEW':len(reviewing),'EXCLUDE':audit_count['EXCLUDE'],
                        'candidate_variants':len(candidates),'candidate_genes':len(genes_all),
                        'candidate_segregation_tiers':dict(tiers),'family_sizes':dict(sizes),
                        'population_unresolved_not_candidates':len(population_unresolved)},
              'warnings':list(review_reason_count)}
    (a.outdir/'RUN_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('\n' + '='*68)
    log(f'COMPLETE: PASS {len(passing):,}, REVIEW {len(reviewing):,}, EXCLUDE {audit_count["EXCLUDE"]:,}; '+
        f'priority candidates {len(candidates):,} variants in {len(genes_all):,} genes')
    log('Historical audit complete; final candidate status requires the local v4.1 reclassification stage.')

if __name__=='__main__':
    try:main()
    except (RuntimeError,ValueError,KeyError) as e:
        print('ERROR: '+str(e),file=sys.stderr);sys.exit(2)
