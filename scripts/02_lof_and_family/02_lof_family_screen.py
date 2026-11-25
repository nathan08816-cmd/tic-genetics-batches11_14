#!/usr/bin/env python3
"""Non-destructive joint B11-14 LoF/coding-splicing catalog and genotype-QC family screen.

Inputs: already site-QC'd split-ALT VCF and single-sample ANNOVAR annotated VCF.
Never labels de novo 'confirmed'; these are only read-quality-aware candidates.
No gnomAD population lookup is performed here, and missing external AF is NOT zero.
"""
import argparse
import collections
import csv
import gzip
import os
from pathlib import Path
import re
import subprocess
import sys
import shutil
from datetime import datetime, timezone
from openpyxl import load_workbook

ROOT = Path(os.environ.get('PROJECT_ROOT', '/path/to/project'))
DATA = Path(os.environ.get('SOURCE_DATA_DIR', '/path/to/source'))
WORK = ROOT / 'filtering/family_batches11_14/lof_joint_20261026'
VCF = WORK / '03_vqsr997_950.originalPASS.splitALT.vcf.gz'
ANNOT = WORK / '04_refgene.hg38_multianno.vcf'
PHENO = ROOT / 'private/phenotypes.xlsx'
OUT = WORK / 'results'
ALLOWED_EXON = {'frameshift_insertion', 'frameshift_deletion', 'stopgain', 'stoploss'}
MISSING = {'', '.', 'NA', 'nan', 'None', 'null', 'NULL'}


def fail(msg):
    raise RuntimeError(msg)


def key(chrom, pos, ref, alt):
    return ((chrom[3:] if chrom.lower().startswith('chr') else chrom).upper(), int(pos), ref.upper(), alt.upper())


def display(k):
    return '-'.join(map(str, k))


def annovarchar(s):
    return (s or '').replace(r'\x3b', ';').replace(r'\x3B', ';')


def parse_info(raw):
    result = {}
    for token in raw.split(';'):
        if '=' in token:
            name, value = token.split('=', 1)
            result[name] = annovarchar(value)
        elif token:
            result[token] = ''
    return result


def annovar_val(m, name):
    return m.get(name + '.refGene', m.get(name + '_refGene', '.'))


def functional_reason(func, exon):
    ex = annovarchar(exon).strip().lower()
    f = annovarchar(func).strip().lower()
    func_parts = {x.strip() for x in f.split(';')}
    splice = 'splicing' in func_parts and not f.startswith('ncrna')
    reasons = []
    if ex in ALLOWED_EXON:
        reasons.append(ex)
    if splice:
        reasons.append('splicing')
    return ';'.join(reasons)


def check_exists(p):
    if not p.is_file() or p.stat().st_size == 0:
        fail(f'Missing or empty: {p}')


def load_annovar(path):
    records = {}
    seen = 0
    header_ok = False
    classes = collections.Counter()
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('##INFO=<ID=Func.refGene') or line.startswith('##INFO=<ID=Func_refGene'):
                header_ok = True
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t', 8)
            if len(p) < 8:
                fail('Malformed ANNOVAR VCF record')
            if ',' in p[4]:
                fail('ANNOVAR input should contain split-ALT records; found multiallelic record')
            seen += 1
            info = parse_info(p[7])
            func = annovar_val(info, 'Func')
            exon = annovar_val(info, 'ExonicFunc')
            gene = annovar_val(info, 'Gene')
            reason = functional_reason(func, exon)
            if not reason:
                continue
            k = key(p[0], p[1], p[3], p[4])
            if k in records:
                fail('Duplicate annotated variant: ' + display(k))
            records[k] = {'_source_contig':p[0], 'Gene_refGene':gene, 'Func_refGene':func,
                          'ExonicFunc_refGene':exon, 'ANNOVAR_function_reason':reason}
            classes[reason] += 1
    if not header_ok:
        fail('No Func.refGene VCF INFO header found. Is this a valid ANNOVAR hg38 annotated VCF?')
    if not records:
        fail('No coding LoF/splicing records found. Check ANNOVAR annotation format and refGene DB.')
    print('ANNOVAR records read:', seen, 'qualifying LoF/splicing ALT records:', len(records), flush=True)
    print('Qualifying functional categories:', dict(classes), flush=True)
    return records, seen


def normalize(s):
    return re.sub('[^a-z0-9]', '', str(s).lower())


def td(v):
    v = str(v).strip()
    if v in ('1','1.0'):
        return 'affected'
    if v in ('0','0.0'):
        return 'unaffected'
    return 'unknown'


def load_phenotypes(path, samples):
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    it = ws.iter_rows(values_only=True)
    hdr = next(it)
    cols = {normalize(n): i for i,n in enumerate(hdr) if n is not None}
    for name in ('subjiddisplay', 'tourettedisorderid','fatherid','motherid'):
        if name not in cols:
            fail('Phenotype workbook missing: ' + name)
    records = {}
    for r in it:
        value = r[cols['subjiddisplay']]
        if value is None:
            continue
        sid = str(value).strip()
        if sid not in samples:
            continue
        obj = {'td':td(r[cols['tourettedisorderid']]),
               'father':r[cols['fatherid']], 'mother':r[cols['motherid']]}
        if sid in records and records[sid] != obj:
            fail('Conflicting duplicate phenotype row for sample: '+sid)
        records[sid] = obj
    wb.close()
    missing = samples-set(records)
    if missing:
        fail(f'Missing phenotype rows for {len(missing)} VCF samples; examples: {list(sorted(missing))[:5]}')
    return records


def resolve_parent(raw, child, samples):
    if raw is None:
        return None
    s = str(raw).strip()
    if s in MISSING or s == '0':
        return None
    if s in samples:
        return s
    if s.endswith('.0') and s[:-2].isdigit():
        s = s[:-2]
    fam = child.rsplit('-',1)[0]
    sid = fam + '-' + s.rsplit('-',1)[-1]
    return sid if sid in samples else None


def load_relationships(path):
    relation = {}
    with open(path, newline='') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            a,b = r['sample1'].strip(),r['sample2'].strip()
            k = frozenset([a,b])
            rel = r['relationship'].strip().lower()
            if k in relation and relation[k] != rel:
                fail('Conflicting relatedness labels')
            relation[k] = rel
    return relation


def load_sample_warnings(path):
    flagged = set()
    with open(path,newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if 'contamination_proxy' in (r.get('qc_flags') or ''):
                flagged.add(r['sample'].strip())
    return flagged


def is_po(label):
    return label in ('parent_offspring','parent-offspring','parent_child','po')


def trios_from_pedigree(records,samples,relations):
    trios = []
    counts = collections.Counter()
    for child, r in records.items():
        fa = resolve_parent(r['father'],child,samples)
        mo = resolve_parent(r['mother'],child,samples)
        if not fa or not mo or len({child,fa,mo}) < 3:
            continue
        counts['all_pedigree_complete'] += 1
        if not (is_po(relations.get(frozenset([child,fa]),'')) and is_po(relations.get(frozenset([child,mo]),''))):
            counts['relationship_flagged'] += 1
            continue
        counts['genetically_supported'] += 1
        if r['td']=='affected':
            counts['affected_child'] += 1
            if records[fa]['td']=='unaffected' and records[mo]['td']=='unaffected':
                counts['affected_child_two_unaffected_parents'] += 1
                trios.append((child,fa,mo))
    return trios,counts


def parse_gt(fmt, value, min_dp, min_gq, ref_max_ab, het_min_ab, het_max_ab, alt_min_reads):
    d = dict(zip(fmt,value.split(':')))
    gt = d.get('GT','.')
    if '.' in gt or not gt or gt in MISSING:
        return ('NO_CALL','missing_GT')
    alleles = re.split('[/|]',gt)
    if not alleles or any(x not in ('0','1') for x in alleles):
        return ('NO_CALL','unresolved_GT')
    try:
        dp = int(d.get('DP','.'))
        gq = float(d.get('GQ','.'))
    except (ValueError, TypeError):
        return ('NO_CALL','missing_DP_GQ')
    if dp < min_dp or gq < min_gq:
        return ('NO_CALL','low_DP_or_GQ')
    ad = d.get('AD','.')
    try:
        abases = [int(x) for x in ad.split(',')]
    except ValueError:
        return ('NO_CALL','missing_AD')
    if len(abases)<2 or min(abases[:2])<0 or sum(abases[:2])==0:
        return ('NO_CALL','missing_AD')
    frac = abases[1] / sum(abases[:2])
    nalt=alleles.count('1')
    if nalt==0:
        if frac > ref_max_ab:
            return ('NO_CALL','ref_alt_reads_too_high')
        return ('REF','')
    if nalt==len(alleles):
        if frac < 0.80 or abases[1]<alt_min_reads:
            return ('NO_CALL','alt_homozygote_balance')
        return ('ALT','')
    if nalt==1:
        if abases[1] < alt_min_reads or not(het_min_ab<=frac<=het_max_ab):
            return ('NO_CALL','het_balance_or_AD')
        return ('ALT','')
    return ('NO_CALL','unexpected_GT')


def opener_gz(path):
    return gzip.open(path,'wt',newline='',encoding='utf-8')


def write_vcf_bgzip(sourcevcf,positions,annotations,out_path,records,seqnames,trios,contam,args,paths):
    p_bed=paths['bed']
    contigs={(k[0], k[1]):v['_source_contig'] for k,v in annotations.items()}
    with open(p_bed,'w') as out:
        for chrom,pos in sorted(positions,key=lambda x: (int(x[0]) if x[0].isdigit() else {'X':23,'Y':24,'MT':25}.get(x[0],99),x[1])):
            # bcftools -T BED semantics: 0-based half-open; retain actual VCF contig spelling.
            actual_chrom=contigs[(chrom,pos)]
            out.write(f'{actual_chrom}\t{pos-1}\t{pos}\n')
    q = subprocess.Popen(['bcftools','view','-T',str(p_bed),'--targets-overlap','pos','-Ov',str(sourcevcf)],
                         stdout=subprocess.PIPE, stderr=None, text=True, bufsize=1024*1024)
    outbinary = open(out_path,'wb')
    bgzip = subprocess.Popen(['bcftools','view','-Oz','-'], stdin=subprocess.PIPE,stdout=outbinary)
    heads=[]
    counts=collections.Counter()
    fnames=['Variant_ID','Chr','Pos','Ref','Alt','Gene_refGene','Func_refGene','ExonicFunc_refGene','ANNOVAR_function_reason',
            'Original_FILTER','Affected_QC_ALT_samples','Unaffected_QC_ALT_samples','Affected_QC_REF_samples',
            'Family_count_affected_carrier','Family_count_strict_compatible','Family_count_strict_inconsistent',
            'Possible_de_novo_children_review','QC_no_call_sample_count','VCF_ALT_frequency_NOTE']
    familyfields=['Variant_ID','Gene_refGene','Family_ID','Affected_ALT','Unaffected_ALT','Affected_REF','Unaffected_REF',
                  'Affected_missing_or_low_quality','Unaffected_missing_or_low_quality','Strict_screen_status']
    dnvfields=['Variant_ID','Gene_refGene','Child_ID','Father_ID','Mother_ID','Evidence_status','QC_NOTE']
    with opener_gz(paths['catalog']) as catalog, opener_gz(paths['rollup']) as rollup, opener_gz(paths['family']) as famout, opener_gz(paths['dnv']) as denovo:
        cw=csv.DictWriter(catalog,delimiter='\t',fieldnames=fnames[:9]);cw.writeheader()
        rw=csv.DictWriter(rollup,delimiter='\t',fieldnames=fnames);rw.writeheader()
        fw=csv.DictWriter(famout,delimiter='\t',fieldnames=familyfields);fw.writeheader()
        dw=csv.DictWriter(denovo,delimiter='\t',fieldnames=dnvfields);dw.writeheader()
        try:
            assert q.stdout is not None and bgzip.stdin is not None
            for line in q.stdout:
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        observed=line.rstrip('\n').split('\t')[9:]
                        if observed != seqnames:
                            fail('VCF sample order changed unexpectedly')
                    bgzip.stdin.write(line.encode());continue
                p=line.rstrip('\n').split('\t')
                if len(p)!=9+len(seqnames):
                    fail('Unexpected genotype count at '+p[0]+':'+p[1])
                if ',' in p[4]:
                    fail('SplitALT VCF still includes multiallelic site at '+p[0]+':'+p[1])
                k=key(p[0],p[1],p[3],p[4])
                a=annotations.get(k)
                if not a:
                    continue
                bgzip.stdin.write(line.encode())
                counts['matched_alleles']+=1
                rec={'Variant_ID':display(k),'Chr':k[0],'Pos':k[1],'Ref':k[2],'Alt':k[3],**{name:value for name,value in a.items() if not name.startswith('_')}}
                cw.writerow(rec)
                fmt=p[8].split(':')
                status={}; reasons=collections.Counter()
                for i,s in enumerate(seqnames):
                    st,why=parse_gt(fmt,p[i+9],args.min_dp,args.min_gq,args.ref_max_ab,
                                    args.het_min_ab,args.het_max_ab,args.alt_min_reads)
                    if s in contam:
                        st='NO_CALL';why='contamination_proxy_review'
                    status[s]=st
                    if why: reasons[why]+=1
                familys=collections.defaultdict(lambda:collections.Counter())
                for sample,st in status.items():
                    t=records[sample]['td']
                    if t=='unknown':continue
                    familiy=sample.rsplit('-',1)[0]
                    familys[familiy][(t,st)]+=1
                n_carrier_fam=n_consistent=n_inconsistent=0
                for fam,fc in familys.items():
                    na=fc[('affected','ALT')]; nu=fc[('unaffected','ALT')]
                    ar=fc[('affected','REF')]; ur=fc[('unaffected','REF')]
                    am=fc[('affected','NO_CALL')]; um=fc[('unaffected','NO_CALL')]
                    if na==0: continue
                    n_carrier_fam+=1
                    if nu or ar:
                        lab='INCONSISTENT_WITH_STRICT_RULE';n_inconsistent+=1
                    elif am or um:
                        lab='INCOMPLETE_GENOTYPES_REVIEW'
                    elif ur:
                        lab='COMPATIBLE_STRICT_SCREEN';n_consistent+=1
                    else:
                        lab='AFFECTED_CARRIER_ONLY'
                    fw.writerow({'Variant_ID':display(k),'Gene_refGene':a['Gene_refGene'],
                                 'Family_ID':fam,'Affected_ALT':na,'Unaffected_ALT':nu,
                                 'Affected_REF':ar,'Unaffected_REF':ur,
                                 'Affected_missing_or_low_quality':am,
                                 'Unaffected_missing_or_low_quality':um,'Strict_screen_status':lab})
                    counts['family_events_'+lab]+=1
                de_novo_n=0
                if k[0].isdigit() and 1<=int(k[0])<=22:
                    for child,fa,mo in trios:
                        if status[child]=='ALT' and status[fa]=='REF' and status[mo]=='REF':
                            de_novo_n+=1
                            dw.writerow({'Variant_ID':display(k),'Gene_refGene':a['Gene_refGene'],
                                         'Child_ID':child,'Father_ID':fa,'Mother_ID':mo,
                                         'Evidence_status':'POTENTIAL_DE_NOVO_REVIEW_ONLY',
                                         'QC_NOTE':f'GT/DP/GQ/AD passed current thresholds; site quality conservative; require manual read and parental confirmation'})
                            counts['possible_de_novo_child_variant_events']+=1
                rw.writerow({**rec,'Original_FILTER':p[6],
                             'Affected_QC_ALT_samples':sum(1 for s,st in status.items() if st=='ALT' and records[s]['td']=='affected'),
                             'Unaffected_QC_ALT_samples':sum(1 for s,st in status.items() if st=='ALT' and records[s]['td']=='unaffected'),
                             'Affected_QC_REF_samples':sum(1 for s,st in status.items() if st=='REF' and records[s]['td']=='affected'),
                             'Family_count_affected_carrier':n_carrier_fam,
                             'Family_count_strict_compatible':n_consistent,
                             'Family_count_strict_inconsistent':n_inconsistent,
                             'Possible_de_novo_children_review':de_novo_n,
                             'QC_no_call_sample_count':sum(1 for x in status.values() if x=='NO_CALL'),
                             'VCF_ALT_frequency_NOTE':'NOT_GNOMAD_FILTERED'})
                counts['variants_with_affected_confident_carrier']+=bool(n_carrier_fam)
                counts['variants_with_possible_de_novo_event']+=bool(de_novo_n)
                if counts['matched_alleles']%1000==0:
                    print('Processed LoF ALT variants:',counts['matched_alleles'],flush=True)
            q.stdout.close()
            bgzip.stdin.close()
            if bgzip.wait()!=0:
                fail('bcftools BGZF compression failed')
            if q.wait()!=0:
                fail('bcftools view failed; inspect run log for details')
        except Exception:
            try:q.kill()
            except Exception:pass
            try:bgzip.kill()
            except Exception:pass
            raise
        finally:
            outbinary.close()
    if counts['matched_alleles'] != len(annotations):
        fail(f'Only {counts["matched_alleles"]} / {len(annotations)} LoF annotated ALT sites matched to normalized VCF. STOP: inspect variant normalization.')
    subprocess.run(['bcftools','index','-t',str(out_path)],check=True)
    return counts


def main():
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--vcf',type=Path,default=VCF)
    pa.add_argument('--annovar-vcf',type=Path,default=ANNOT)
    pa.add_argument('--phenotype',type=Path,default=PHENO)
    pa.add_argument('--data-dir',type=Path,default=DATA)
    pa.add_argument('--outdir',type=Path,default=OUT)
    pa.add_argument('--min-dp',type=int,default=10)
    pa.add_argument('--min-gq',type=int,default=20)
    pa.add_argument('--alt-min-reads',type=int,default=4)
    pa.add_argument('--het-min-ab',type=float,default=0.25)
    pa.add_argument('--het-max-ab',type=float,default=0.75)
    pa.add_argument('--ref-max-ab',type=float,default=0.05)
    args=pa.parse_args()
    for p in (args.vcf,args.annovar_vcf,args.phenotype,
              args.data_dir/'sample_qc/relatedness_pairs.tsv',args.data_dir/'sample_qc/sample_qc.tsv'):
        check_exists(p)
    for tool in ('bcftools',):
        if shutil.which(tool) is None:
            fail('Required executable missing: '+tool)
    args.outdir.mkdir(parents=True,exist_ok=True)
    output={'bed':args.outdir/'lof_sites_for_vcf.bed',
            'catalog':args.outdir/'batch1114_lof_annovar_variant_catalog.tsv.gz',
            'rollup':args.outdir/'batch1114_lof_variant_rollup.tsv.gz',
            'family':args.outdir/'batch1114_lof_family_variant_screen.tsv.gz',
            'dnv':args.outdir/'batch1114_lof_potential_de_novo_review.tsv.gz',
            'vcf':args.outdir/'batch1114_lof_QCsite_candidates.with_genotypes.vcf.gz',
            'summary':args.outdir/'batch1114_lof_family_screen_summary.txt'}
    for p in output.values():
        if p.exists():fail('Refusing to overwrite existing analysis output: '+str(p))
    if Path(str(output['vcf'])+'.tbi').exists():
        fail('Existing candidate VCF index; refusing to overwrite')
    annotated,n_seen=load_annovar(args.annovar_vcf)
    q=subprocess.run(['bcftools','query','-l',str(args.vcf)],capture_output=True,text=True,check=True)
    samples=q.stdout.splitlines()
    if len(samples)!=1245 or len(set(samples))!=len(samples):
        fail(f'Expected 1245 unique VCF samples, found {len(samples)}')
    records=load_phenotypes(args.phenotype,set(samples))
    counts_td=collections.Counter(r['td'] for r in records.values())
    relation=load_relationships(args.data_dir/'sample_qc/relatedness_pairs.tsv')
    contam=load_sample_warnings(args.data_dir/'sample_qc/sample_qc.tsv')
    trios,trio_stats=trios_from_pedigree(records,set(samples),relation)
    print('Phenotype counts:',dict(counts_td),flush=True)
    print('Pedigree trio reconciliation:',dict(trio_stats),flush=True)
    print('Affected-child, unaffected-parent complete supported trios:',len(trios),flush=True)
    print('Contamination_proxy samples marked unusable for interpretation:',len(contam),flush=True)
    positions={(k[0],k[1]) for k in annotated}
    counts=write_vcf_bgzip(args.vcf,positions,annotated,output['vcf'],records,samples,trios,contam,args,output)
    with open(output['summary'],'w') as f:
        f.write('BATCHES 11-14 JOINT CALL: PRELIMINARY LoF/SPLICING QUALITY & FAMILY REVIEW\n')
        f.write(f'UTC date: {datetime.now(timezone.utc).isoformat()}\n')
        f.write('Build: GRCh38; Twist +/-20 bp; site filter: VQSR SNP99.7/indel95 intersection with Yale original PASS\n')
        f.write(f'Inputs: {args.vcf} ; {args.annovar_vcf} ; {args.phenotype}\n')
        f.write(f'ANNOVAR site records: {n_seen}; LoF/splicing unique ALT: {len(annotated)}\n')
        f.write(f'Phenotypes: {dict(counts_td)}\n')
        f.write(f'Trio reconciliation: {dict(trio_stats)}\n')
        f.write(f'Supported affected-child/two-unaffected-parent trios: {len(trios)}\n')
        f.write(f'Contamination-proxy flagged samples: {len(contam)}\n')
        f.write(f'GT QC: min DP={args.min_dp}; min GQ={args.min_gq}; ALT reads>={args.alt_min_reads}; heterozygous allele balance [{args.het_min_ab},{args.het_max_ab}]; reference max alternate fraction {args.ref_max_ab}; hom-alt min alt fraction 0.80\n')
        for k,v in counts.most_common(): f.write(f'{k}: {v}\n')
        f.write('IMPORTANT: These candidate variants have NOT been run through the gnomAD v4 API filter. Affected carrier is not proven segregation.\n')
        f.write('Potential de novo is NOT validated: check parental depth/read evidence, allelic balance, complex regions and sample contamination before interpretation.\n')
        f.write('All prior original files remain unchanged.\n')
    print('SUCCESS: '+str(output['summary']),flush=True)
    for n,p in output.items():print(n,p)


if __name__=='__main__':
    try:main()
    except Exception as e:
        print('FAILED:',str(e),file=sys.stderr)
        raise
