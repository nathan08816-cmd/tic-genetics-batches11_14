#!/usr/bin/env python3
"""Audit all families in the phenotype workbook.

No sequencing membership is consulted. All subjects in the source are included,
including Evaluated=N; an additional Evaluated=Y-only summary is generated.
Structural anomalies are REVIEW FLAGS, not automatically corrected relationships.
Standalone Python 3, standard library only; source Excel is never modified.

Usage:
 python3 full_phenotype_family_audit_20261009.py \
   --phenotype '/path/to/phenotypes_current.xlsx' \
   --old-phenotype '/path/to/2025-06-18_Clinical_Phenotype_Data.xlsx' \
   --output-dir /path/to/results
"""
import argparse, csv, os, re, zipfile, xml.etree.ElementTree as ET
from collections import defaultdict, Counter

DATE = '20261009'
NS = {'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}

class XLSXReader:
    def __init__(self,path):
        self.z=zipfile.ZipFile(path)
        self.strings=[]
        if 'xl/sharedStrings.xml' in self.z.namelist():
            root=ET.fromstring(self.z.read('xl/sharedStrings.xml'))
            for si in root.findall('m:si',NS):
                self.strings.append(''.join(el.text or '' for el in si.iter('{%s}t'%NS['m'])))
        book=ET.fromstring(self.z.read('xl/workbook.xml'))
        rel=ET.fromstring(self.z.read('xl/_rels/workbook.xml.rels'))
        paths={x.attrib['Id']:x.attrib['Target'] for x in rel}
        self.sheets={}
        for s in book.find('m:sheets',NS):
            target=paths[s.attrib['{%s}id'%NS['r']]]
            if target.startswith('/'):target=target.lstrip('/')
            elif not target.startswith('xl/'):target='xl/'+target
            self.sheets[s.attrib['name']]=target
    def rows(self,sheet):
        root=ET.fromstring(self.z.read(self.sheets[sheet]))
        hdr=None
        for r in root.findall('.//m:sheetData/m:row',NS):
            cells={}
            for c in r.findall('m:c',NS):
                s=re.match('([A-Z]+)',c.attrib.get('r',''))
                if not s:continue
                col=0
                for char in s.group(1):col=col*26+ord(char)-64
                v=c.find('m:v',NS)
                inline=c.find('m:is',NS)
                if c.attrib.get('t')=='s' and v is not None: val=self.strings[int(v.text)]
                elif c.attrib.get('t')=='inlineStr' and inline is not None:
                    val=''.join(t.text or '' for t in inline.iter('{%s}t'%NS['m']))
                elif v is not None:val=v.text or ''
                else:val=''
                cells[col]=val
            if hdr is None:
                hdr=[str(cells.get(i,'')).strip() for i in range(1,max(cells)+1)]
                continue
            yield dict(zip(hdr,[str(cells.get(i,'')).strip() for i in range(1,len(hdr)+1)]))

def family_of(sid):
    parts=sid.split('-')
    return '-'.join(parts[:3]) if len(parts)>=4 else ''

def local(sid):return sid.split('-')[-1]

def status(r):return {'1':'Affected','0':'Unaffected'}.get(r.get('Tourette DisorderID',''),'Unknown')

def write_tsv(path,fields,rows):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        for row in rows:w.writerow({k:row.get(k,'') for k in fields})

def pedigree_analysis(family,rows,ids):
    """Pedigree relations from recorded FatherID/MotherID, never from naming alone."""
    sample={r['SubjID_Display']:r for r in rows}
    loc={r['SubjID']:r['SubjID_Display'] for r in rows if r.get('SubjID')}
    triples=[]
    graph={x:set() for x in sample}
    edge_info=[]
    for r in rows:
        child=r['SubjID_Display']
        pf=r.get('FatherID','');pm=r.get('MotherID','')
        father=loc.get(pf,'');mother=loc.get(pm,'')
        if father and mother and father!=mother and child not in (father,mother):
            triples.append((child,father,mother))
        for role,plink,person in [('Father',pf,father),('Mother',pm,mother)]:
            if plink:
                edge_info.append((child,role,plink,person))
            if person:
                graph[child].add(person)
                graph[person].add(child)
    comps=[];seen=set()
    for node in sorted(ids):
        if node in seen:continue
        stack=[node];seen.add(node);comp=[]
        while stack:
            a=stack.pop();comp.append(a)
            for b in graph[a]:
                if b not in seen:seen.add(b);stack.append(b)
        comps.append(sorted(comp))
    comps.sort(key=lambda x:x[0])
    confirmed=[t for t in triples if set(t).issubset(ids)]
    def is_simplex(t):
        c,f,m=t
        return status(sample[c])=='Affected' and status(sample[f])==status(sample[m])=='Unaffected'
    return confirmed,comps,edge_info,is_simplex

def derive_type(ids,rows,confirmed,is_simplex):
    n=len(ids)
    if n==0:return 'No Data'
    if n==1:
        s=status(rows[0])
        return {'Affected':'Singleton','Unaffected':'Singleton-Unaffected','Unknown':'Singleton-Unknown'}[s]
    if n==3 and confirmed:
        return 'Simplex' if any(is_simplex(t) for t in confirmed) else 'Other Trio'
    if sum(status(r)=='Affected' for r in rows)>=2:return 'Multiplex'
    return 'Other Family'

def potential_trios(fam,rows,confirmed):
    """Conservative heuristic for exactly three available members only.
    Require a younger apparent child, adjacent-generation male/female adults,
    and missing parent fields. Does not prove any biological relationship.
    """
    if confirmed or len(rows)!=3:return []
    candidates=[]
    for r in rows:
        child=r.get('SubjID','')
        if not re.fullmatch(r'[3-8]\d{3}',child):continue
        gen=int(child[0]);parents=[p for p in rows if re.fullmatch(r'[3-8]\d{3}',p.get('SubjID',''))
                 and int(p['SubjID'][0])==gen-1]
        males=[p for p in parents if p.get('GenderDescr')=='Male']
        females=[p for p in parents if p.get('GenderDescr')=='Female']
        if not males or not females:continue
        fa=r.get('FatherID','');mo=r.get('MotherID','')
        plausible=(not fa or not mo)
        if not plausible:continue
        if (fa and not any(p['SubjID']==fa for p in males)) or (mo and not any(p['SubjID']==mo for p in females)):
            continue
        if len(males)>2 or len(females)>2:continue # unusually broad combinations remain in other audits
        candidates.append({'Family_ID':fam,'Child_ID':r['SubjID_Display'],
          'Possible_Fathers':';'.join(p['SubjID_Display'] for p in males),
          'Possible_Mothers':';'.join(p['SubjID_Display'] for p in females),
          'Recorded_FatherID':fa,'Recorded_MotherID':mo,
          'Reason':'Heuristic adjacent generation + sexes; parental field(s) missing; not a confirmed trio'})
    return candidates

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--phenotype',required=True)
    ap.add_argument('--old-phenotype',default='')
    ap.add_argument('--output-dir',required=True)
    args=ap.parse_args();os.makedirs(args.output_dir,exist_ok=True)
    new=XLSXReader(args.phenotype)
    sheet='Sheet1' if 'Sheet1' in new.sheets else next(iter(new.sheets))
    rows=list(new.rows(sheet))
    required=['SubjID_Display','FamilyID','SubjID','FatherID','MotherID','ProbandID','Evaluated','Tourette DisorderID','GenderDescr']
    for col in required:
        if not rows or col not in rows[0]:raise ValueError('Missing column '+col)
    byfamily=defaultdict(list)
    seen={};issues=[];possible=[];disconnected=[];families=[]
    issuecounts=Counter()
    def flag(fid,sid,code,evidence,severity='Review'):
        issues.append({'Family_ID':fid,'Sample_ID':sid,'Issue_Type':code,'Severity':severity,'Evidence':evidence})
        issuecounts[code]+=1
    for i,r in enumerate(rows,2):
        sid=r['SubjID_Display'];fam=family_of(sid)
        if not sid or len(sid.split('-'))!=4:
            flag(fam,sid,'MALFORMED_SAMPLE_ID',f'Excel row {i}: {sid!r}','Data error')
            continue
        if sid in seen:flag(fam,sid,'DUPLICATE_SUBJECT_ID',f'rows {seen[sid]} and {i}','Data error')
        seen[sid]=i
        if r['SubjID']!=local(sid) and not (r['SubjID'].isdigit() and local(sid).isdigit() and int(r['SubjID'])==int(local(sid))):
            flag(fam,sid,'LOCAL_ID_MISMATCH',f"SubjID={r['SubjID']} vs ID suffix={local(sid)}",'Data error')
        if r['FamilyID']!=sid.split('-')[2]:flag(fam,sid,'FAMILY_ID_MISMATCH',f"FamilyID={r['FamilyID']}",'Data error')
        try:
            expected=f"{int(r['SiteID']):03d}-{int(r['SubsiteID']):02d}-{r['FamilyID']}-{r['SubjID']}"
            if sid!=expected and not (sid.rsplit('-',1)[0]==expected.rsplit('-',1)[0] and local(sid).isdigit() and r['SubjID'].isdigit() and int(local(sid))==int(r['SubjID'])):
                flag(fam,sid,'SOURCE_ID_FIELDS_MISMATCH',f'expected {expected}','Data error')
        except (ValueError,KeyError):flag(fam,sid,'BAD_SITE_IDENTIFIERS','Non-numeric or missing SiteID/SubsiteID','Review')
        if r['Tourette DisorderID'] not in ('0','1','9'):flag(fam,sid,'UNEXPECTED_TD_CODE',str(r['Tourette DisorderID']),'Review')
        if r['Evaluated'] not in ('Y','N'):flag(fam,sid,'UNEXPECTED_EVALUATED_CODE',str(r['Evaluated']),'Review')
        if r['FatherID'] and r['FatherID']==r['MotherID']:
            flag(fam,sid,'FATHER_EQUALS_MOTHER',f"ID={r['FatherID']}",'Data error')
        if r['SubjID'] and (r['SubjID']==r['FatherID'] or r['SubjID']==r['MotherID']):
            flag(fam,sid,'SELF_AS_PARENT',f"father={r['FatherID']} mother={r['MotherID']}",'Data error')
        byfamily[fam].append(r)

    allfamilytypes=Counter();evalfamilytypes=Counter();full_count=Counter();n_missing_parent_ref=Counter()
    for fam,frows in sorted(byfamily.items()):
        fullids={r['SubjID_Display'] for r in frows}
        locals_by_id=defaultdict(list)
        for r in frows:locals_by_id[r['SubjID']].append(r)
        for id_,a in locals_by_id.items():
            if len(a)>1:flag(fam,'','DUPLICATE_LOCAL_SUBJECT_ID',f'{id_} repeats {len(a)} times','Data error')
        probands=[r for r in frows if r['ProbandID']=='1']
        if len(probands)!=1:
            flag(fam,'','PROBAND_COUNT_NOT_ONE',f'ProbandID=1 count {len(probands)} in {len(frows)}-member family')
        confirmed,comps,links,is_simplex=pedigree_analysis(fam,frows,fullids)
        for child,role,localid,person in links:
            if not person:
                n_missing_parent_ref[role]+=1
                flag(fam,child,'REFERENCED_PARENT_NOT_IN_DATASET',f'{role}ID={localid} not found in this family; does not prove an error')
            else:
                par=next(r for r in frows if r['SubjID_Display']==person)
                expected='Male' if role=='Father' else 'Female'
                if par.get('GenderDescr')!=expected:
                    flag(fam,child,'PARENT_SEX_CONFLICT',f'{role} {person} is {par.get("GenderDescr")}','Data error')
                clocal=local(child)
                if re.fullmatch(r'\d{4}',clocal) and re.fullmatch(r'\d{4}',localid):
                    if int(localid[0])!=int(clocal[0])-1:
                        flag(fam,child,'UNUSUAL_GENERATION_NUMBERING',f'{role}ID={localid} child={clocal}; numbering is only heuristic')
        # Directed cycle check with recorded in-family parent links.
        graph={r['SubjID_Display']:[] for r in frows}
        for c,_,_,par in links:
            if par and par!=c:graph[c].append(par) # self-parent already reported independently
        colors={};cycle_seen=set()
        def visit(n,path):
            colors[n]=1
            for p in graph[n]:
                if colors.get(p)==1:cycle_seen.add(n)
                elif not colors.get(p):visit(p,path+[n])
            colors[n]=2
        for s in graph:
            if not colors.get(s):visit(s,[])
        if cycle_seen:flag(fam,'','PEDIGREE_CYCLE_NONSELF','Directed child-to-parent cycle excluding self-parent loops','Data error')
        if len(comps)>1:
            parts=' | '.join(';'.join(comp) for comp in comps)
            disconnected.append({'Family_ID':fam,'Member_Count':len(frows),'Component_Count':len(comps),
                                 'Components':parts,'Has_Recorded_Trio':'Yes' if confirmed else 'No',
                                 'Review_Note':'Disconnected in recorded pedigree graph; may be missing ancestors or independent branches'})
            flag(fam,'','DISCONNECTED_COMPONENTS',f'{len(comps)} components among {len(frows)} members')
        poss=potential_trios(fam,frows,confirmed)
        possible.extend(poss)
        for p in poss:flag(fam,p['Child_ID'],'POSSIBLE_MISSING_PARENT_LINK',f"possible father(s) {p['Possible_Fathers']} and mother(s) {p['Possible_Mothers']}; fields {p['Recorded_FatherID']}/{p['Recorded_MotherID']}")
        fulltype=derive_type(fullids,frows,confirmed,is_simplex)
        allfamilytypes[fulltype]+=1
        ys=[r for r in frows if r['Evaluated']=='Y']
        yids={r['SubjID_Display'] for r in ys}
        confy=[tr for tr in confirmed if set(tr).issubset(yids)]
        ytype=derive_type(yids,ys,confy,is_simplex)
        evalfamilytypes[ytype]+=1
        ct=Counter(status(r) for r in frows)
        cy=Counter(status(r) for r in ys)
        if ct['Affected']==0:full_count['families_no_TD_affected']+=1
        if len(confirmed)>0 and len(frows)>3:full_count['families_with_trio_inside_larger_family']+=1
        if len(confirmed)>1:full_count['families_with_multiple_trios']+=1
        families.append({'Family_ID':fam,'All_Individual_Count':len(frows),
           'Evaluated_Y_Count':len(ys),'All_Affected':ct['Affected'],'All_Unaffected':ct['Unaffected'],'All_Unknown':ct['Unknown'],
           'Full_Family_Type':fulltype,'Evaluated_Y_Family_Type':ytype,
           'Evaluated_Affected':cy['Affected'],'Evaluated_Unaffected':cy['Unaffected'],'Evaluated_Unknown':cy['Unknown'],
           'Full_Confirmed_Trio_Count':len(confirmed),'Confirmed_Trio_Children':';'.join(t[0] for t in confirmed),
           'Connected_Component_Count':len(comps),'Possible_Missing_Parent_Child_Count':len(poss),
           'ProbandID_1_Count':len(probands),
           'Sample_IDs':';'.join(sorted(fullids))})

    prefix=lambda name:os.path.join(args.output_dir,name+'_20261009.tsv')
    fam_fields=list(families[0])
    write_tsv(prefix('full_phenotype_family_structure'),fam_fields,families)
    type_rows=[]
    order=['No Data','Singleton','Singleton-Unaffected','Singleton-Unknown','Simplex','Other Trio','Multiplex','Other Family']
    for scope,counts in [('All individuals (including Evaluated=N)',allfamilytypes),('Evaluated=Y only',evalfamilytypes)]:
        for typ in order:type_rows.append({'Scope':scope,'Family_Type':typ,'Families':counts[typ]})
        type_rows.append({'Scope':scope,'Family_Type':'TOTAL','Families':sum(counts.values())})
    write_tsv(prefix('full_phenotype_family_type_summary'),['Scope','Family_Type','Families'],type_rows)
    write_tsv(prefix('full_phenotype_audit_issues'),['Family_ID','Sample_ID','Issue_Type','Severity','Evidence'],issues)
    write_tsv(prefix('full_phenotype_disconnected_families'),list(disconnected[0]) if disconnected else ['Family_ID','Member_Count','Component_Count','Components','Has_Recorded_Trio','Review_Note'],disconnected)
    write_tsv(prefix('full_phenotype_possible_missing_parent_trios'),['Family_ID','Child_ID','Possible_Fathers','Possible_Mothers','Recorded_FatherID','Recorded_MotherID','Reason'],possible)

    historical=[];nhist={}
    if args.old_phenotype:
        old=XLSXReader(args.old_phenotype)
        osheet='DistributionData' if 'DistributionData' in old.sheets else next(iter(old.sheets))
        orows=list(old.rows(osheet))
        oldmap={r['SubjID_Display']:r for r in orows}
        newmap={r['SubjID_Display']:r for r in rows}
        hist_cols=['FatherID','MotherID','ProbandID','GenderDescr','Evaluated','Tourette DisorderID']
        for sid in sorted(oldmap.keys() & newmap.keys()):
            for field in hist_cols:
                oldval=oldmap[sid].get(field,'');newval=newmap[sid].get(field,'')
                if oldval!=newval:
                    historical.append({'Family_ID':family_of(sid),'Sample_ID':sid,'Field':field,'Old_Value':oldval,'New_Value':newval})
        nhist={'old_unique_subjects':len(oldmap),'new_unique_subjects':len(newmap),
               'new_subjects':len(newmap.keys()-oldmap.keys()),
               'old_subjects_not_in_new':len(oldmap.keys()-newmap.keys()),
               'shared_subjects':len(oldmap.keys()&newmap.keys())}
        write_tsv(prefix('old_vs_new_phenotype_changes'),['Family_ID','Sample_ID','Field','Old_Value','New_Value'],historical)

    priority_map={
      'SELF_AS_PARENT':1,'PEDIGREE_CYCLE_NONSELF':1,'FATHER_EQUALS_MOTHER':1,
      'PARENT_SEX_CONFLICT':1,'DUPLICATE_SUBJECT_ID':1,'DUPLICATE_LOCAL_SUBJECT_ID':1,
      'LOCAL_ID_MISMATCH':1,'FAMILY_ID_MISMATCH':1,'SOURCE_ID_FIELDS_MISMATCH':1,
      'POSSIBLE_MISSING_PARENT_LINK':2,'DISCONNECTED_COMPONENTS':2,
      'PROBAND_COUNT_NOT_ONE':3,'REFERENCED_PARENT_NOT_IN_DATASET':3,
      'UNUSUAL_GENERATION_NUMBERING':3}
    oldidset=set(oldmap) if args.old_phenotype else set()
    before_by_family=defaultdict(int)
    for sid in oldidset:before_by_family[family_of(sid)]+=1
    now_by_family=Counter(family_of(r['SubjID_Display']) for r in rows)
    priority_rows=[]
    for rec in issues:
        q=dict(rec)
        q['Priority_Level']=priority_map.get(rec['Issue_Type'],3)
        q['Source_Cohort']=rec['Family_ID'].split('-')[2][:2] if len(rec['Family_ID'].split('-'))>=3 else ''
        q['Row_Present_In_Older_Workbook']=('Yes' if rec['Sample_ID'] in oldidset else ('No' if rec['Sample_ID'] else 'Not a sample-level flag'))
        q['Family_Has_New_Records']=('Yes' if now_by_family[rec['Family_ID']]>before_by_family[rec['Family_ID']] else 'No') if args.old_phenotype else ''
        priority_rows.append(q)
    priority_rows.sort(key=lambda x:(x['Priority_Level'],x['Issue_Type'],x['Family_ID'],x['Sample_ID']))
    write_tsv(prefix('full_phenotype_priority_review'),
        ['Priority_Level','Family_ID','Sample_ID','Issue_Type','Severity','Evidence','Source_Cohort','Row_Present_In_Older_Workbook','Family_Has_New_Records'],priority_rows)

    text=[]
    def add(s=''):text.append(s)
    add('FULL PHENOTYPE PEDIGREE AUDIT — October 9, 2026')
    add('Source: '+args.phenotype)
    add(f'Individuals: {len(rows)}; unique IDs: {len(set(r["SubjID_Display"] for r in rows))}; families: {len(families)}')
    add(f'Evaluated=Y: {sum(r["Evaluated"]=="Y" for r in rows)}; Evaluated=N: {sum(r["Evaluated"]=="N" for r in rows)}')
    add('Tourette DisorderID raw counts: '+repr(dict(sorted(Counter(r['Tourette DisorderID'] for r in rows).items()))))
    add('Status = Tourette DisorderID only (1 affected, 0 unaffected, otherwise unknown)')
    for scope,counts in [('All individuals',allfamilytypes),('Evaluated=Y only',evalfamilytypes)]:
        add('\nFAMILY TYPES — '+scope)
        for k in order:add(f'  {k}: {counts[k]}')
        add(f'  TOTAL: {sum(counts.values())}')
    add('\nREVIEW FLAGS (not necessarily confirmed source errors)')
    for typ,n in issuecounts.most_common():add(f'  {typ}: {n}')
    add(f'Families with disconnected components: {len(disconnected)}')
    add(f'Heuristic possible missing-parent child rows: {len(possible)}')
    for k,v in full_count.items():add(f'{k}: {v}')
    add('Notes: parent references missing from the workbook may reflect relatives not enrolled; disconnected branches are not automatically errors.')
    add('Possible missing-parent trios use subject numbering and parent sexes as heuristic only; they are not confirmed relationships.')
    add('Exact 3-member trio types are based on FatherID/MotherID; families with >3 people remain Multiplex if >=2 affected.')
    if nhist:
        add('\nOLD VS NEW FILE')
        for k,v in nhist.items():add(f'{k}: {v}')
        add('Changes by field: '+repr(dict(Counter(x['Field'] for x in historical))))
        add('Previous family-subset comparison is a separate result; do not infer old family types from this full file.')
    summary_path=os.path.join(args.output_dir,'full_phenotype_audit_summary_20261009.txt')
    open(summary_path,'w',encoding='utf-8').write('\n'.join(text)+'\n')
    print('\n'.join(text))
    print('\nOutputs:',args.output_dir)

if __name__=='__main__':main()
