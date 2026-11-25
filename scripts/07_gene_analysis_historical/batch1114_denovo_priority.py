#!/usr/bin/env python3
"""Create de novo-centered review sheets for Batch 11-14 (GRCh38).

Requires a PRELIMINARY DE NOVO events TSV/TSV.gz as input. Never infers de novo
status from NDD/ClinVar annotations. Leaves all source files unchanged.

prepare: deduplicate by exact ALT allele and create gnomAD-ready unique list.
rank: join optional exact-allele annotation files incl. gnomAD export and
      output review-ranked events/alleles. Scores are TRIAGE, not pathogenicity.
"""
import argparse, csv, gzip, re
from collections import Counter, defaultdict
from pathlib import Path

MISSING = {'', '.', 'NA', 'N/A', 'NONE', 'NULL', 'NAN'}
REQUIRED = ['Variant_ID','Chr','Pos','Ref','Alt']
AF_LIMIT = 0.001

def read(path):
    path=Path(path)
    op=gzip.open if path.suffix=='.gz' else open
    with op(path,'rt',newline='',encoding='utf-8-sig') as f:
        rd=csv.DictReader(f,delimiter='\t')
        if not rd.fieldnames: raise ValueError(f'No header: {path}')
        return list(rd.fieldnames),list(rd)

def first(row, names):
    lower={k.lower():v for k,v in row.items()}
    for name in names:
        v=lower.get(name.lower())
        if v is not None and str(v).strip().upper() not in MISSING:
            return str(v).strip()
    return ''

def allele(row, fname):
    vid=first(row,['Variant_ID','variant_id','allele_id'])
    parts=None
    if vid:
        m=re.fullmatch(r'(?:chr)?([0-9]+|X|Y|MT|M)[-:]([0-9]+)[-:]([ACGT]+)[-:]([ACGT]+)',vid,re.I)
        if not m: raise ValueError(f'{fname}: invalid Variant_ID {vid!r}')
        parts=list(m.groups())
    else:
        parts=[first(row,['Chr','CHROM','#CHROM','chromosome']),first(row,['Pos','POS','position','Start']),first(row,['Ref','REF']),first(row,['Alt','ALT'])]
        if not all(parts): raise ValueError(f'{fname}: expected Variant_ID or Chr,Pos,Ref,Alt')
    c,p,r,a=parts
    c=c.removeprefix('chr').upper();r=r.upper();a=a.upper()
    if not p.isdecimal() or not re.fullmatch('[ACGT]+',r) or not re.fullmatch('[ACGT]+',a):
        raise ValueError(f'{fname}: invalid GRCh38 SNV/indel allele {parts}')
    return f'{c}-{int(p)}-{r}-{a}',c,str(int(p)),r,a

def write(path, fields, rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader();w.writerows(rows)
    print(f'Wrote {path} ({len(rows)} rows)')

def maybe_num(s):
    try:
        v=float(s)
        return v if 0<=v<=1 else None
    except (ValueError,TypeError):return None

def metadata(row):
    return {'Gene':first(row,['Gene','Gene_refGene','Gene.refGene','Gene_Name','Symbol','gene_symbol']),
            'Consequence':first(row,['Consequence','ExonicFunc_refGene','ExonicFunc.refGene','Func_refGene','Func.refGene','Effect','Variant_Classification','Annotation']),
            'ClinVar':first(row,['ClinVar','CLNSIG','ClinVar_CLNSIG','ClinVar_significance']),
            'NDD_overlap':first(row,['NDD_overlap','NDD','NDD_gene','is_NDD_gene'])}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['prepare','rank'])
    ap.add_argument('--denovo',required=True,help='Event-level preliminary de novo TSV(.gz)')
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--annotations',nargs='*',default=[],help='Optional allele-matched tables including gnomAD export')
    args=ap.parse_args();out=Path(args.outdir)
    fields,source=read(args.denovo)
    if not source:raise ValueError('Empty de novo event input')
    events=[];groups=defaultdict(list)
    for i,r in enumerate(source,2):
        v,c,p,ref,alt=allele(r,f'{args.denovo}:{i}')
        new=dict(r);new.update({'Variant_ID':v,'Chr':c,'Pos':p,'Ref':ref,'Alt':alt})
        events.append(new);groups[v].append(new)
    unique=[]
    for vid,rr in groups.items():
        item={k:rr[0][k] for k in REQUIRED}
        item['Preliminary_event_count']=len(rr)
        unique.append(item)
    print(f'PRELIMINARY events={len(events)}, exact unique alleles={len(unique)}')
    if args.action=='prepare':
        write(out/'denovo_unique_for_gnomad.tsv',REQUIRED+['Preliminary_event_count'],unique)
        write(out/'denovo_all_preliminary_events.tsv',list(dict.fromkeys(REQUIRED+fields)),events)
        return
    lookup=defaultdict(dict)
    for path in args.annotations:
        af,ar=read(path)
        matches=0
        for j,r in enumerate(ar,2):
            try:vid,*_=allele(r,f'{path}:{j}')
            except ValueError:continue
            if vid not in groups:continue
            matches+=1
            for key,val in r.items():
                if key in REQUIRED or val is None or str(val).strip().upper() in MISSING:continue
                # Never silently overwrite discrepant annotation values.
                existing=lookup[vid].get(key)
                if existing and existing!=val:
                    raise ValueError(f'Conflicting annotation for {vid}, field {key}: {existing!r} vs {val!r}')
                lookup[vid][key]=val
        print(f'Exact allele matches from {path}: {matches}')
    ranked=[]
    for base in unique:
        vid=base['Variant_ID'];ann=lookup[vid]
        # Fall back to event's own annotation fields when missing in optional tables.
        joined={**groups[vid][0],**ann}
        for k,v in metadata(joined).items():
            if not first(joined,[k]):joined[k]=v
        consequence=first(joined,['Consequence','ExonicFunc_refGene','ExonicFunc.refGene','Func_refGene','Func.refGene','Effect']).lower()
        cln=first(joined,['ClinVar','CLNSIG','ClinVar_CLNSIG']).lower()
        ndd=first(joined,['NDD_overlap','NDD','NDD_gene','is_NDD_gene']).lower()
        decision=first(joined,['gnomAD_v4_decision']).upper()
        status=first(joined,['gnomAD_v4_status']).upper()
        af=maybe_num(first(joined,['gnomAD_v4_max_observed_AF']))
        # Use maximum of components if aggregate absent.
        afs=[maybe_num(first(joined,[k])) for k in ('gnomAD_v4_exome_AF','gnomAD_v4_genome_AF','gnomAD_v4_population_max_AF')]
        if af is None and any(x is not None for x in afs):af=max(x for x in afs if x is not None)
        high_impact=any(s in consequence for s in ['frameshift','stopgain','stop_gained','splice_acceptor','splice_donor','splicing','stoploss','stop_lost']) and 'ncrna' not in consequence
        pathogenic=bool(re.search(r'(?<![a-z])pathogenic|likely.pathogenic',cln)) and not any(s in cln for s in ['conflict','uncertain','benign'])
        ndd_yes=ndd in ('yes','true','1','y')
        # Status remains unresolved even when an allele is absent from gnomAD.
        if decision=='EXCLUDE' or af is not None and af>=AF_LIMIT:
            population='EXCLUDE_OR_REVIEW_COMMON'
        elif decision=='RETAIN' and af is not None and af<AF_LIMIT:
            population='PASS_AF_AND_FLAGS'
        elif decision=='RETAIN_FOR_REVIEW_NOT_FOUND' or status=='NOT_FOUND':
            population='NOT_FOUND_AF_UNKNOWN'
        elif decision=='REVIEW_REPRESENTATION':population='CHECK_NORMALIZATION'
        elif decision=='REVIEW_NO_FREQUENCY':population='NO_USABLE_FREQUENCY'
        else:population='NOT_QUERIED_OR_UNRESOLVED'
        points=(4 if high_impact else 0)+(2 if pathogenic else 0)+(1 if ndd_yes else 0)+(2 if population=='PASS_AF_AND_FLAGS' else 0)
        review=('POPULATION_EXCLUSION_REVIEW' if population=='EXCLUDE_OR_REVIEW_COMMON' else
                'PRIORITY_REVIEW' if high_impact and population=='PASS_AF_AND_FLAGS' else
                'VERIFY_GNOMAD_AND_TRIO' if high_impact else 'SECONDARY_REVIEW')
        ranked.append({**base,**joined,**ann,'Preliminary_event_count':base['Preliminary_event_count'],
                       'max_gnomAD_AF_used':af if af is not None else '',
                       'Population_review_status':population,'High_impact_annotation':high_impact,
                       'NDD_evidence':ndd_yes,'ClinVar_PLP_label':pathogenic,
                       'Triage_score_not_pathogenicity':points,'Review_category':review,
                       'De_novo_status':'PRELIMINARY_UNVALIDATED'})
    ranked.sort(key=lambda r:(r['Review_category']=='POPULATION_EXCLUSION_REVIEW',
                              -r['Triage_score_not_pathogenicity'],r['Variant_ID']))
    finalfields=list(dict.fromkeys(REQUIRED+['Preliminary_event_count','Gene','Consequence','ClinVar','NDD_overlap']+
        [k for r in ranked for k in r]+['max_gnomAD_AF_used','Population_review_status',
        'High_impact_annotation','NDD_evidence','ClinVar_PLP_label','Triage_score_not_pathogenicity',
        'Review_category','De_novo_status']))
    write(out/'denovo_unique_prioritized.tsv',finalfields,ranked)
    by_id={r['Variant_ID']:r for r in ranked}
    eventfields=list(dict.fromkeys(finalfields+fields))
    write(out/'denovo_events_prioritized.tsv',eventfields,[{**r,**by_id[r['Variant_ID']]} for r in events])
    print('Review categories:',dict(Counter(r['Review_category'] for r in ranked)))
    print('CAUTION: de novo calls not validated; do not interpret missing gnomAD as AF=0.')

if __name__=='__main__':main()
