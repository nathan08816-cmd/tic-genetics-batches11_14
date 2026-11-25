#!/usr/bin/env python3
"""Synthetic unit and integration tests; no participant-level research data."""
import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def module(path):
    spec=importlib.util.spec_from_file_location('tic_test_'+Path(path).stem, ROOT/path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

stage02=module('scripts/02_lof_and_family/02_lof_family_screen.py')
stage03=module('scripts/03_local_annotation/gnomad_flag_evidence_from_vcf.py')
stage04=module('scripts/04_local_v41_filter/refinalize_local_gnomad_v41.py')
stage05=module('scripts/05_clinvar_priority/filter_clinvar_plp_pass.py')
trio=module('scripts/06_trio_validation/validate_five_trios_from_vcf.py')

class PipelineTests(unittest.TestCase):
    def test_split_genes(self):
        self.assertEqual(stage05.split_genes('GENE1,GENE2;GENE3|GENE1'), ['GENE1','GENE2','GENE3'])
        self.assertEqual(stage04.genes_of('GENE1,GENE2;GENE3|GENE1'), ['GENE1','GENE2','GENE3'])

    def test_clinical_bucket(self):
        self.assertEqual(stage04.clinical_bucket('Pathogenic [111]'), 'P/LP')
        self.assertEqual(stage04.clinical_bucket('Conflicting_classifications_of_pathogenicity [222]'), 'Conflicting')
        self.assertEqual(stage04.clinical_bucket('Likely_benign'), 'Benign/LB')
        self.assertEqual(stage05.label({'ClinVar':'Pathogenic/Likely_pathogenic [4]'}), 'Pathogenic/Likely_pathogenic')

    def test_genotype_qc_existing_thresholds(self):
        fmt=['GT','DP','GQ','AD']
        gt=lambda s:stage02.parse_gt(fmt,s,10,20,0.05,0.25,0.75,4)[0]
        self.assertEqual(gt('0/1:30:99:15,15'),'ALT')
        self.assertEqual(gt('0/0:20:99:19,1'),'REF') # QC REF but not conclusive de novo
        self.assertEqual(gt('0/0:20:99:18,2'),'NO_CALL')
        self.assertEqual(gt('0/1:30:99:26,4'),'NO_CALL')
        self.assertEqual(gt('./.:30:99:30,0'),'NO_CALL')
        self.assertEqual(gt('1:30:99:0,30'),'ALT') # haploid, ploidy REVIEW in stage 04

    def test_parent_alt_support(self):
        fmt=['GT','DP','GQ','AD']
        r=trio.parse_sample(fmt,'0/0:20:99:19,1')
        self.assertEqual(r['QC_call'],'REF')
        self.assertEqual(len(trio.parent_reference_warnings('father_unaffected',r)),1)
        self.assertEqual(trio.parent_reference_warnings('affected_child',r),[])
        r2=trio.parse_sample(fmt,'0/0:20:99:20,0')
        self.assertEqual(trio.parent_reference_warnings('mother_unaffected',r2),[])
        hap=trio.parse_sample(fmt,'1:20:99:0,20')
        self.assertIn('Haploid',hap['QC_note'])

    def test_local_assessment_pass_common_ac0_missing(self):
        def make(filt,af,an,flag):
            return {'Flag_Source':'echtvar_gnomad_v4.1+vcfanno_exact_allele',
                    'Local_gnomAD_AF_joint':str(af) if af is not None else '',
                    'Local_gnomAD_AF_grpmax':'', 'Local_gnomAD_AF_exomes':'',
                    'Local_gnomAD_AF_genomes':'',
                    'Local_gnomAD_AF':str(af) if af is not None else '',
                    'Local_gnomAD_AN':str(an),'Local_gnomAD_Filter':filt,
                    'Flag_Assessment':flag}
        self.assertEqual(stage04.local_assessment(make('PASS',0.0002,1000,'PASS'),'1-10-A-G')[0], 'PASS')
        self.assertEqual(stage04.local_assessment(make('PASS',0.001,1000,'EXCLUDE_AF'),'1-11-A-G')[0],'EXCLUDE_AF')
        self.assertEqual(stage04.local_assessment(make('AC0',0,1000,'UNRESOLVED'),'1-12-A-G')[0],'AC0_REVIEW')
        self.assertEqual(stage04.local_assessment(make('MISSING',None,-1,'UNRESOLVED'),'1-13-A-G')[0],'UNOBSERVED')
        self.assertEqual(stage04.local_assessment(make('EXOMES_FILTERED',0.0001,1000,'FLAGGED'),'1-14-A-G')[0],'FLAGGED')
        with self.assertRaises(ValueError):
            stage04.local_assessment(make('PASS',0.001,1000,'PASS'),'1-15-A-G')

    def test_local_ac0_annotation_extraction(self):
        required=['gnomad_filter','gnomad_AN','gnomad_AF_grpmax','clinvar_sig','clinvar_id','gnomad_lcr','gnomad_segdup']
        with tempfile.TemporaryDirectory() as td:
            x=Path(td); inp=x/'test.vcf';out=x/'out.tsv'
            hdr='##fileformat=VCFv4.2\n'+''.join('##INFO=<ID='+k+',Number=1,Type=String,Description="test">\n' for k in required)
            hdr+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            inp.write_text(hdr+'1\t10\t.\tA\tG\t.\tPASS\tgnomad_filter=AC0;gnomad_AN=200;gnomad_AF_grpmax=0;clinvar_sig=Pathogenic;clinvar_id=9;gnomad_lcr=0;gnomad_segdup=0\n')
            run=subprocess.run([sys.executable,str(ROOT/'scripts/03_local_annotation/gnomad_flag_evidence_from_vcf.py'),
                                '--vcf',str(inp),'--output',str(out),'--expected-variants','1'],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            with out.open() as f:row=next(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(row['Flag_Assessment'],'UNRESOLVED')
            self.assertIn('AC0',row['Flag_Detail'])

    def test_clinvar_split_and_review_separation(self):
        fields=['Variant_ID','Gene','Function','ClinVar','ClinVar_Source','ClinVar_Disease','ClinVar_Review_Status',
                'ClinVar_Stars','gnomAD_AF','Internal_AF','Affected_Carriers','Unaffected_Carriers',
                'Affected_Families','Strict_Families','Inconsistent_Families','Segregation','Filter_Status',
                'Filter_Reason','Local_gnomAD_Filter','Local_LCR','Local_SegDup']
        def entry(i,gene,status,cl='Pathogenic [42]'):
            return dict(Variant_ID=f'1-{i}-A-G',Gene=gene,Function='stopgain',ClinVar=cl,
                        ClinVar_Source='exact_allele_vcfanno',ClinVar_Disease='example disease',
                        ClinVar_Review_Status='single submitter',ClinVar_Stars='1',gnomAD_AF='0.0001',
                        Internal_AF='0.0005',Affected_Carriers='1',Unaffected_Carriers='0',
                        Affected_Families='1',Strict_Families='1',Inconsistent_Families='0',
                        Segregation='Strict-compatible',Filter_Status=status,Filter_Reason='unobserved',
                        Local_gnomAD_Filter='MISSING',Local_LCR='0',Local_SegDup='0')
        rows=[entry(10,'GENE1,GENE2','PASS'),entry(11,'GENE2;GENE3','PASS', 'Likely_pathogenic [43]'),
              entry(12,'GENE4','REVIEW'),entry(13,'GENE5','EXCLUDE')]
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);a=root/'audit.tsv';e=root/'events.tsv';dest=root/'output'
            with a.open('w',newline='') as fh:
                w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
            with e.open('w',newline='') as fh:
                w=csv.DictWriter(fh,fieldnames=['Variant_ID','Family_ID','Affected_ALT'],delimiter='\t');w.writeheader()
                for x in rows:w.writerow({'Variant_ID':x['Variant_ID'],'Family_ID':'F'+x['Variant_ID'],'Affected_ALT':'1'})
            run=subprocess.run([sys.executable,str(ROOT/'scripts/05_clinvar_priority/filter_clinvar_plp_pass.py'),
                                '--audit',str(a),'--events',str(e),'--outdir',str(dest),'--include-review'],
                                capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            with (dest/'ClinVar_P_or_LP_PASS_Genes.tsv').open() as fh:g=list(csv.DictReader(fh,delimiter='\t'))
            self.assertEqual(set(x['Gene'] for x in g),{'GENE1','GENE2','GENE3'})
            with (dest/'ClinVar_P_or_LP_REVIEW_Variants.tsv').open() as fh:v=list(csv.DictReader(fh,delimiter='\t'))
            self.assertEqual(len(v),1)
            self.assertEqual(v[0]['Filter_Status'],'REVIEW')
            self.assertEqual(v[0]['Gene'],'GENE4')

if __name__=='__main__': unittest.main()
