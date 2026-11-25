#!/usr/bin/env python3
"""Refinalize Batches 11-14 using only locally annotated gnomAD v4.1.

Inputs are existing validated 14,124-row QC/internal-AF/ANNOVAR/ClinVar audit and
exact-allele echtvar/vcfanno table. Previous gnomAD API results and historical
ANNOVAR AF fields are NOT criteria. No web/API queries, no data mutation.

With --vcf (required for complete gene lists), recomputes per-gene UNIQUE
individual carriers from the existing split-ALT genotype VCF with the original
QC function; checks each variant's affected/unaffected carrier counts and
internal AC/AN against prior audit before writing final gene lists.
"""
import argparse, collections, csv, gzip, importlib.util, json, math, re, sys
from pathlib import Path
from datetime import datetime, timezone

V_FIELDS=['Variant_ID','Gene','Function','gnomAD_AF','gnomAD_QC','Internal_AF',
          'Affected_Carriers','Unaffected_Carriers','Affected_Families',
          'Strict_Families','Inconsistent_Families','Segregation','ClinVar',
          'Filter_Status','Filter_Reason']
G_FIELDS=['Gene','Candidate_Variants','Variant_Types','Affected_Carriers',
          'Unaffected_Carriers','Affected_Families','Strict_Segregation_Variants',
          'Lowest_gnomAD_AF','ClinVar_Summary']
TARGET_VARIANTS=14124
CUTOFF=.001
INTERNAL=.015


def read_tsv(path):
    with open(path,encoding='utf-8-sig',newline='') as f:
        rd=csv.DictReader(f,delimiter='\t')
        if not rd.fieldnames: raise ValueError(f'No TSV header: {path}')
        for row in rd:
            if None in row: raise ValueError(f'Extra TSV fields: {path}')
            yield row


def write_tsv(path, fields, records):
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader();w.writerows(records)


def fval(v, name, absent_ok=True):
    if v is None or str(v).strip() in ('','.','NA','N/A','MISSING','-1'):
        return None if absent_ok else (_ for _ in ()).throw(ValueError(f'{name}: missing frequency'))
    try: x=float(v)
    except (TypeError,ValueError): raise ValueError(f'{name}: invalid float {v!r}')
    if not (math.isfinite(x) and 0 <= x <= 1):raise ValueError(f'{name}: invalid AF {v!r}')
    return x


def integer(v,name):
    try: return int(v)
    except (TypeError,ValueError): raise ValueError(f'{name}: invalid integer {v!r}')


def fmt(x): return '' if x is None else f'{x:.9g}'


def genes_of(g):
    return sorted({x.strip() for x in re.split('[;,|]',g) if x.strip() and x.strip()!='.'})


def key_of(chrom,pos,ref,alt):
    chrom=chrom.upper().removeprefix('CHR')
    return f'{chrom}-{int(pos)}-{ref.upper()}-{alt.upper()}'


def clinical_bucket(text):
    s=text.lower()
    if s=='not annotated':return 'Not annotated'
    if 'conflict' in s:return 'Conflicting'
    if 'benign' in s and 'pathogenic' in s:return 'Conflicting'
    if 'pathogenic' in s:return 'P/LP'
    if 'uncertain' in s or re.search(r'\bvus\b',s):return 'VUS'
    if 'benign' in s:return 'Benign/LB'
    return 'Other'


def local_assessment(r,allele):
    if r.get('Flag_Source','')!='echtvar_gnomad_v4.1+vcfanno_exact_allele':
        raise ValueError(f'{allele}: local gnomAD source not verified')
    afs=[fval(r.get(k,''),allele+':'+k) for k in (
       'Local_gnomAD_AF_joint','Local_gnomAD_AF_grpmax',
       'Local_gnomAD_AF_exomes','Local_gnomAD_AF_genomes')]
    af=max((v for v in afs if v is not None),default=None)
    recorded=fval(r.get('Local_gnomAD_AF',''),allele+':derived max AF')
    if (af is None)!=(recorded is None) or (af is not None and abs(af-recorded)>max(1e-8,af*2e-5)):
        raise ValueError(f'{allele}: max AF does not match extracted summary: {af} vs {recorded}')
    an=integer(r.get('Local_gnomAD_AN'),allele+':AN')
    filt=r.get('Local_gnomAD_Filter','').strip()
    if filt=='AC0':outcome='AC0_REVIEW'
    elif an<=0 or af is None: outcome='UNOBSERVED'
    elif filt in ('','.','MISSING','UNKNOWN'):outcome='UNRESOLVED_FILTER'
    elif filt!='PASS':outcome='FLAGGED'
    elif af>=CUTOFF:outcome='EXCLUDE_AF'
    else:outcome='PASS'
    mapped={'UNOBSERVED':'UNRESOLVED','AC0_REVIEW':'UNRESOLVED','UNRESOLVED_FILTER':'UNRESOLVED',
            'FLAGGED':'FLAGGED','EXCLUDE_AF':'EXCLUDE_AF','PASS':'PASS'}[outcome]
    if mapped!=r.get('Flag_Assessment'):
        raise ValueError(f'{allele}: contradictory saved Flag_Assessment {r.get("Flag_Assessment")!r} vs computed {outcome}')
    return outcome,af,filt,an


def tier(r):
    if integer(r['Affected_Carriers'],'Affected_Carriers')<=0:return 'No affected carrier'
    ev=r['Original_Family_Evidence']
    if ev=='STRICT_SUPPORTED_IN_AT_LEAST_ONE_FAMILY':return 'Strict-compatible'
    if ev=='AFFECTED_CARRIER_ONLY':return 'Affected-carrier-only'
    return 'Family review'


def reconstruct(audit, local):
    if len(audit)!=TARGET_VARIANTS or len(local)!=TARGET_VARIANTS or set(audit)!=set(local):
        raise ValueError(f'Variant coverage mismatch: audit={len(audit)}, local={len(local)}, need identical {TARGET_VARIANTS} keys')
    out=[];full=[];counts=collections.Counter()
    for key,r in audit.items():
        if r.get('Variant_ID')!=key:raise ValueError('Bad key '+key)
        lr=local[key]
        if not r.get('Gene') or not r.get('Function'):raise ValueError(key+': ANNOVAR annotation missing')
        if r.get('ClinVar_Source') not in ('exact_allele_vcfanno','absent_exact_allele_vcfanno'):
            raise ValueError(key+': ClinVar not confirmed as exact-allele matched')
        outcome,af,filter_text,an=local_assessment(lr,key)
        saved_af=fval(r.get('Local_AF_v4.1'),key+':audit local AF')
        if (af is None)!=(saved_af is None) or (af is not None and abs(af-saved_af)>max(1e-8,af*2e-5)):
            raise ValueError(key+': audit/local AF mismatch')
        if r.get('Local_Flag_Status')!=lr.get('Flag_Assessment'):
            raise ValueError(key+': audit/local filter mismatch')
        internal=fval(r.get('Internal_AF'),key+':internal')
        affected=integer(r.get('Affected_Carriers'),key+':affected')
        reasons=[]
        if affected<=0:reasons.append('No affected carrier')
        if internal is None:reasons.append('No callable internal AF')
        elif internal>INTERNAL:reasons.append('Internal AF > 0.015')
        if outcome=='EXCLUDE_AF':reasons.append('gnomAD v4.1 max AF >= 0.001')
        elif outcome=='FLAGGED':reasons.append('gnomAD v4.1 FILTER='+filter_text)
        elif outcome=='UNOBSERVED':reasons.append('Absent/unmatched in local gnomAD v4.1; no positive frequency evidence')
        elif outcome=='UNRESOLVED_FILTER':reasons.append('Missing local gnomAD quality filter')
        elif outcome=='AC0_REVIEW':reasons.append('gnomAD AC0: no QC-passing ALT observations; review before inclusion')
        if affected<=0 or internal is None or (internal is not None and internal>INTERNAL) or outcome in ('EXCLUDE_AF','FLAGGED'):
            status='EXCLUDE'
        else:
            if r.get('Chromosome_scope')!='AUTOSOMAL':reasons.append('Sex/other chromosome: ploidy review')
            if lr.get('Local_LCR')=='1':reasons.append('Low-complexity region: mapping review')
            if lr.get('Local_SegDup')=='1':reasons.append('Segmental duplication: mapping review')
            status='PASS' if not reasons else 'REVIEW'
        if status=='PASS':reasons=['Local gnomAD v4.1, internal AF and QC passed']
        line={
           'Variant_ID':key,'Gene':r['Gene'],'Function':r['Function'],
           'gnomAD_AF':fmt(af),'gnomAD_QC':('PASS' if outcome in ('PASS','EXCLUDE_AF') else filter_text if outcome in ('FLAGGED','AC0_REVIEW') else 'UNRESOLVED'),
           'Internal_AF':fmt(internal),'Affected_Carriers':affected,
           'Unaffected_Carriers':integer(r['Unaffected_Carriers'],'Unaffected_Carriers'),
           'Affected_Families':integer(r['Affected_Families'],'Affected_Families'),
           'Strict_Families':integer(r['Strict_Families'],'Strict_Families'),
           'Inconsistent_Families':integer(r['Inconsistent_Families'],'Inconsistent_Families'),
           'Segregation':tier(r),'ClinVar':r['ClinVar'],
           'Filter_Status':status,'Filter_Reason':'; '.join(reasons)
        }
        out.append(line)
        full.append({**line,**{col:r[col] for col in ('Internal_ALL_AC','Internal_ALL_AN',
                   'Internal_Unrelated_AF','Internal_Unaffected_AF','Internal_Unaffected_Unrelated_AF',
                   'ClinVar_Source','ClinVar_Disease','ClinVar_Review_Status','ClinVar_Stars',
                   'Local_LCR','Local_SegDup','Chromosome_scope','Original_Family_Evidence')},
            'Local_gnomAD_AN':an,'Local_gnomAD_Filter':filter_text,
            'Local_gnomAD_AF_joint':lr.get('Local_gnomAD_AF_joint',''),
            'Local_gnomAD_AF_grpmax':lr.get('Local_gnomAD_AF_grpmax',''),
            'Local_gnomAD_AF_exomes':lr.get('Local_gnomAD_AF_exomes',''),
            'Local_gnomAD_AF_genomes':lr.get('Local_gnomAD_AF_genomes',''),
            'Local_gnomAD_fafmax95':lr.get('Local_gnomAD_fafmax95',''),
            'Local_gnomAD_exomes_filters':lr.get('Local_gnomAD_exomes_filters',''),
            'Local_gnomAD_genomes_filters':lr.get('Local_gnomAD_genomes_filters',''),
            'Old_API_Used_for_Decision':'NO','Old_ANNOVAR_AF_Used_for_Decision':'NO'})
        counts[status]+=1;counts['local_'+outcome]+=1
        if affected>0:counts['affected_carrier_total']+=1
    return out,full,counts


def aggregate_gene_rows(variants, carriers):
    idx=collections.defaultdict(lambda:{'alleles':set(),'types':set(),'a':set(),'u':set(),'f':set(),'strict':0,'afs':[], 'cv':collections.Counter()})
    for v in variants:
        k=v['Variant_ID']
        if k not in carriers:raise ValueError('Missing VCF carrier sample evidence for '+k)
        a,u,f=carriers[k]
        for gene in genes_of(v['Gene']):
            d=idx[gene];d['alleles'].add(k);d['types'].update(v['Function'].split(';'))
            d['a'].update(a);d['u'].update(u);d['f'].update(f)
            d['strict']+=int(v['Strict_Families']>0)
            if v['gnomAD_AF']:d['afs'].append(float(v['gnomAD_AF']))
            d['cv'][clinical_bucket(v['ClinVar'])]+=1
    answer=[]
    for gene,d in idx.items():
        cv='; '.join(f'{k}: {d["cv"][k]}' for k in ('P/LP','Conflicting','VUS','Benign/LB','Other','Not annotated') if d['cv'][k])
        answer.append({'Gene':gene,'Candidate_Variants':len(d['alleles']),
         'Variant_Types':'; '.join(sorted(d['types'])),'Affected_Carriers':len(d['a']),
         'Unaffected_Carriers':len(d['u']),'Affected_Families':len(d['f']),
         'Strict_Segregation_Variants':d['strict'],
         'Lowest_gnomAD_AF':fmt(min(d['afs']) if d['afs'] else None),
         'ClinVar_Summary':cv or 'Not annotated'})
    return sorted(answer,key=lambda x:(-x['Strict_Segregation_Variants'],-x['Candidate_Variants'],x['Gene']))


def recompute_carriers(vcf,audit,ids,original_code,sample_qc,phenotype):
    spec=importlib.util.spec_from_file_location('old_jointcall_family_qc',original_code)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    masks=mod.load_sample_warnings(sample_qc)
    if len(masks)!=9:raise ValueError(f'Expected 9 contamination-proxy samples, found {len(masks)}')
    carrier={};seen=set()
    op=gzip.open if str(vcf).endswith('.gz') else open
    with op(vcf,'rt') as f:
        samples=None;families=None;tds=None
        for line in f:
            if line.startswith('##'):continue
            if line.startswith('#CHROM'):
                samples=line.rstrip('\n').split('\t')[9:]
                if len(samples)!=1245 or len(set(samples))!=len(samples):raise ValueError('Unexpected VCF sample list')
                td=mod.load_phenotypes(phenotype,set(samples))
                tds=[td[s]['td'] for s in samples]
                families=[s.rsplit('-',1)[0] for s in samples]
                if collections.Counter(tds)!=collections.Counter({'affected':699,'unaffected':534,'unknown':12}):
                    raise ValueError('Phenotype/sample totals differ from previously audited source')
                if len(set(families))!=659:raise ValueError('Expected 659 families')
                continue
            if line.startswith('#'):continue
            if samples is None:raise ValueError('VCF lacks #CHROM sample header')
            t=line.rstrip('\n').split('\t');key=key_of(t[0],t[1],t[3],t[4])
            if key not in audit:raise ValueError('Unmatched allele in VCF '+key)
            if key in seen:raise ValueError('Duplicate allele in VCF '+key)
            seen.add(key)
            if t[6]!='PASS':raise ValueError('VCF FILTER not PASS at '+key)
            if key not in ids:continue
            fmt=t[8].split(':');gt_i=fmt.index('GT')
            af=set();uf=set();fams=set();ac=0;an=0
            for i,sid in enumerate(samples):
                if sid in masks:continue
                call,_=mod.parse_gt(fmt,t[i+9],10,20,.05,.25,.75,4)
                if call=='NO_CALL':continue
                geno=t[i+9].split(':')[gt_i]
                alleles=re.split(r'[/|]',geno)
                if not alleles or any(z not in ('0','1') for z in alleles):raise ValueError('Unexpected genotype '+key)
                dose=alleles.count('1');ac+=dose;an+=len(alleles)
                if (dose>0)!=(call=='ALT'):raise ValueError('ALT genotype QC conflict '+key)
                if call=='ALT':
                    if tds[i]=='affected':af.add(sid);fams.add(families[i])
                    elif tds[i]=='unaffected':uf.add(sid)
            r=audit[key]
            for k,n in (('Affected_Carriers',len(af)),('Unaffected_Carriers',len(uf)),('Affected_Families',len(fams)),('Internal_ALL_AC',ac),('Internal_ALL_AN',an)):
                if integer(r[k],key+':'+k)!=n:raise ValueError(f'VALIDATION STOP {key} {k}: audit {r[k]} vs VCF {n}')
            carrier[key]=(af,uf,fams)
    if seen!=set(audit):raise ValueError(f'Full VCF coverage mismatch {len(seen)}/{len(audit)}')
    if set(carrier)!=set(ids):raise ValueError('Candidate sample coverage incomplete')
    return carrier


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--local-evidence',type=Path,required=True)
    p.add_argument('--vcf',type=Path,help='Required to produce true deduplicated gene lists')
    p.add_argument('--qc-code',type=Path,help='Original 02_lof_family_screen.py')
    p.add_argument('--sample-qc',type=Path)
    p.add_argument('--phenotype',type=Path)
    p.add_argument('--family-summary',type=Path)
    p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--preview-variants-only',action='store_true',help='Produce allele tables only; never claim gene list is finalized')
    args=p.parse_args()
    if args.outdir.exists():raise ValueError('Output dir exists; refusing to overwrite '+str(args.outdir))
    if not args.preview_variants_only:
        if any(not v for v in (args.vcf,args.qc_code,args.sample_qc,args.phenotype)):
            raise ValueError('Full mode requires --vcf --qc-code --sample-qc --phenotype')
    aud={};local={}
    for r in read_tsv(args.audit):
        k=r['Variant_ID'];
        if k in aud:raise ValueError('Duplicate audit variant '+k)
        aud[k]=r
    for r in read_tsv(args.local_evidence):
        k=r['Variant_ID'];
        if k in local:raise ValueError('Duplicate local variant '+k)
        local[k]=r
    variants,full,counts=reconstruct(aud,local)
    passing=[r for r in variants if r['Filter_Status']=='PASS']
    review=[r for r in variants if r['Filter_Status']=='REVIEW']
    selected=passing+review
    genes_pass=[];genes_review=[];genes_all=[]
    if not args.preview_variants_only:
        print('Validating VCF sample calls and deduplicating gene carriers...',flush=True)
        carrier=recompute_carriers(args.vcf,aud,{r['Variant_ID'] for r in selected},args.qc_code,args.sample_qc,args.phenotype)
        genes_pass=aggregate_gene_rows(passing,carrier)
        genes_review=aggregate_gene_rows(review,carrier)
        genes_all=aggregate_gene_rows(selected,carrier)
    args.outdir.mkdir(parents=True)
    write_tsv(args.outdir/'jointcall_candidate_variants_final.tsv',V_FIELDS,passing)
    write_tsv(args.outdir/'jointcall_candidate_variants_need_review.tsv',V_FIELDS,review)
    write_tsv(args.outdir/'jointcall_candidate_variants_pass_and_review.tsv',V_FIELDS,selected)
    write_tsv(args.outdir/'gnomAD_unobserved_variants_review.tsv',V_FIELDS,[r for r in review if 'Absent/unmatched in local gnomAD' in r['Filter_Reason']])
    write_tsv(args.outdir/'region_or_ploidy_variants_review.tsv',V_FIELDS,[r for r in review if ('region:' in r['Filter_Reason'] or 'chromosome:' in r['Filter_Reason'])])
    if not args.preview_variants_only:
        write_tsv(args.outdir/'jointcall_candidate_genes_final.tsv',G_FIELDS,genes_pass)
        write_tsv(args.outdir/'jointcall_candidate_genes_need_review.tsv',G_FIELDS,genes_review)
        write_tsv(args.outdir/'jointcall_candidate_genes_pass_and_review.tsv',G_FIELDS,genes_all)
    write_tsv(args.outdir/'ALL_14124_LOCAL_V41_REASSESSED_AUDIT.tsv',list(full[0]),full)
    if args.family_summary:
        rows=list(read_tsv(args.family_summary))
        if sum(integer(r['Sequenced_Individuals'],'family count') for r in rows)!=1245 or len(rows)!=659:
            raise ValueError('Family summary mismatch')
        write_tsv(args.outdir/'family_summary.tsv',list(rows[0]),rows)
    notes=[
      '# TIC Genetics Batches 11–14 — final local gnomAD v4.1 screen',
      '',f'Generated {datetime.now().strftime("%Y-%m-%d %H:%M")}',
      '', '## Input and policy',
      'Previous QC/LoF annotation, exact-allele ClinVar, internal-AF and family evidence were retained from the audited output.',
      'Only local echtvar gnomAD v4.1 INFO annotations determine external AF/quality.',
      'AF is MAX of local overall, exome, genome and ancestry-group maximum; PASS requires gnomAD_AN>0, FILTER=PASS and maximum AF<0.001.',
      'Internal cohort AF: exclude if >0.015 (cohort-level screening cutoff).',
      'Unmatched/absent gnomAD alleles are REVIEW, never assigned an observed AF of 0.',
      'LCR/segmental-duplication and non-autosomal variants require additional review, not automatic exclusion.',
      'Family segregation is an evidence tier; no NDD, recurrence, de novo, or ClinVar pathogenicity requirement.',
      'Historical API/v4.0 decisions and old ANNOVAR AF fields are NOT used for inclusion or exclusion.',
      'ANNOVAR consequences and gene names are retained; ClinVar is exact-allele matched by vcfanno.',
      '', '## Counts',
      f'Total LoF/splicing variants: {len(variants):,}',
      f'PASS: {len(passing):,}',f'REVIEW: {len(review):,}',f'EXCLUDE: {counts["EXCLUDE"]:,}',
      f'Local rarity+QC pass before chromosome and mapping review: {counts["local_PASS"]:,} (includes variants lacking an affected carrier and internal AF exclusions)',
      f'PASS genes: {len(genes_pass):,}' if genes_pass else 'PASS genes: not calculated in variant-only preview',
      f'PASS+REVIEW genes: {len(genes_all):,}' if genes_all else 'PASS+REVIEW genes: not calculated in variant-only preview',
      '', '## Important limitations',
      'PASS is filtering evidence, not a causal Tourette diagnosis or validated disease-gene association.',
      'This step does not re-run Yale raw calling/VQSR or ClinVar retrieval, and uses the lab-managed gnomAD v4.1 reference snapshot.',
      'For final gene carrier counts the server must validate the existing joint-call genotype VCF (full mode).',
      'Variant-only preview does not generate a gene list because multiple variants can share carriers.',
    ]
    (args.outdir/'METHODS_AND_RESULTS_LOCAL_V41.md').write_text('\n'.join(notes)+'\n',encoding='utf-8')
    (args.outdir/'RUN_MANIFEST_LOCAL_V41.json').write_text(json.dumps({
      'generated_utc':datetime.now(timezone.utc).isoformat(),
      'inputs':{'audit':str(args.audit),'local_evidence':str(args.local_evidence),'vcf':str(args.vcf) if args.vcf else None},
      'mode':'variant_only_preview' if args.preview_variants_only else 'full_vcf_verified',
      'thresholds':{'gnomad_local_v4.1_max_AF_lt':CUTOFF,'gnomad_filter':'PASS',
                    'internal_full_AF_exclude_gt':INTERNAL},
      'counts':{**counts,'PASS_genes':len(genes_pass) if not args.preview_variants_only else None,
                'candidate_genes':len(genes_all) if not args.preview_variants_only else None},
      'old_api_used':False,'old_annovar_frequency_used':False
    },indent=2),encoding='utf-8')
    print(f'DONE — {len(passing)} PASS, {len(review)} REVIEW, {counts["EXCLUDE"]} EXCLUDE',flush=True)
    if not args.preview_variants_only:print('GENES',len(genes_pass),'PASS and',len(genes_all),'PASS+REVIEW',flush=True)
    print('Folder:',args.outdir,flush=True)

if __name__=='__main__':
    try:main()
    except (ValueError,KeyError,RuntimeError) as e:
        print('ERROR:',e,file=sys.stderr);sys.exit(2)
