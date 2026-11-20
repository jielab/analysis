#!/usr/bin/env python3
"""Focused validation regressions. Uses synthetic reads/SAM; no reference downloads.

Run: python f/validation_selftest.py
These tests exercise panel accounting and alignment interpretation, not biological
accuracy or performance of a real BWA-MEM2/minimap2 installation.
"""
from collections import Counter
from contextlib import ExitStack, redirect_stdout
import csv
import gzip
import io
import json
from pathlib import Path
import sqlite3
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import human_validate as human
import plant_validate as plant
import qc_human_screen
import species_refine
from fastq_stream import fragments


class Tax:
    def __init__(self):
        self.rank={1:'no rank',10:'family',11:'genus',20:'family'}
        self.name={}
    def resolve(self,t):return int(t)
    def lineage(self,t):
        t=int(t)
        if t==0:return iter(())
        if t==1:return iter((1,))
        if t in (10,20):return iter((t,1))
        if t==11:return iter((11,10,1))
        return iter((t,11,10,1))


def fastqs(directory,count,base='plant/species_other',low_last=False):
    files=[directory/(base+f'.R{mate}.fastq.gz') for mate in (1,2)]
    files[0].parent.mkdir(parents=True,exist_ok=True)
    for mate,path in enumerate(files,1):
        with gzip.open(path,'wt') as out:
            for n in range(1,count+1):
                seq=('A'*100 if low_last and n==count else 'ACGTTGCAAGTCGATCAGGA'*5)
                out.write(f'@r{n}/{mate}\n{seq}\n+\n'+('I'*100)+'\n')
    return files


def sample_fixture(root,tids=(101,102,103),group='plant'):
    folder=root/'split/d/s';meta=folder/'qc';meta.mkdir(parents=True)
    base=group+'/species_other';reads=fastqs(folder,len(tids),base)
    rows=[];assignments=[]
    for n,tid in enumerate(tids,1):
        row=dict(assigned_taxid=str(tid),species_taxid=str(tid),group=group,
                 status='assigned_species',bin=base,fragments='1')
        rows.append(row)
        assignments.append(dict(ordinal=str(n),r1_read_id=f'r{n}/1',r2_read_id=f'r{n}/2',
                                **{k:v for k,v in row.items() if k!='fragments'}))
    path=meta/'read-assignments.tsv.gz'
    with gzip.open(path,'wt',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(assignments[0]),delimiter='\t')
        writer.writeheader();writer.writerows(assignments)
    outputs={str(p.relative_to(folder)):dict(state=dict(size=p.stat().st_size)) for p in [path,*reads]}
    (meta/'complete.json').write_text(json.dumps(dict(outputs=outputs)))
    summary=dict(sample_id='d/s',paired=True,fragments=len(tids))
    return folder,summary,rows


def rep_fixture(root,tids,group='plant',bases=100):
    reps={}
    for tid in tids:
        path=root/f'ref{tid}.fa'
        path.write_text('>chr chromosome 1\n'+('ACGTTGCAAGTCGATCAGGA'*5)+'\n')
        reps[tid]=dict(species_taxid=tid,name=str(tid),group=group,bases=bases,
                       path=str(path),resolved_path=str(path))
    return reps


def sam_record(q,mate,contig='T101|R1',mapq=60,secondary=False,score=200):
    flag=(99 if mate==1 else 147)+(256 if secondary else 0)
    position=1 if mate==1 else 151
    partner=151 if mate==1 else 1
    return '\t'.join(map(str,[q,flag,contig,position,mapq,'100M','=',partner,
                              250 if mate==1 else -250,'A'*100,'I'*100,
                              'NM:i:0','AS:i:'+str(score)]))+'\n'


class PanelTests(unittest.TestCase):
    def test_all_77_targets_are_accounted_for_without_one_giant_panel(self):
        with tempfile.TemporaryDirectory() as d:
            reps=rep_fixture(Path(d),range(100,177))
            claims=[dict(ordinal=i+1,species_taxid=t,assigned_taxid=t) for i,t in enumerate(reps)]
            panels,assigned,untested=plant.plan_panels(claims,Counter({t:1 for t in reps}),reps,Tax(),[],60,100000)
            self.assertGreater(len(panels),1)
            self.assertEqual(len(assigned),77)
            self.assertTrue(all(assigned.values()))
            self.assertFalse(untested)
            self.assertTrue(all(len(p['selected'])<=60 and p['declared_bases']<=100000 for p in panels))

    def test_oversize_missing_and_broad_are_visible(self):
        with tempfile.TemporaryDirectory() as d:
            reps=rep_fixture(Path(d),[101,102]);reps[102]['bases']=2000
            claims=[dict(ordinal=1,species_taxid=101,assigned_taxid=101),
                    dict(ordinal=2,species_taxid=102,assigned_taxid=102),
                    dict(ordinal=3,species_taxid=999,assigned_taxid=999),
                    dict(ordinal=4,species_taxid=0,assigned_taxid=10),
                    dict(ordinal=5,species_taxid=0,assigned_taxid=20)]
            _,assigned,untested=plant.plan_panels(claims,Counter(),reps,Tax(),[],10,1000)
            self.assertTrue(assigned[1]);self.assertTrue(assigned[4])
            self.assertIsNone(assigned[2]);self.assertIsNone(assigned[3]);self.assertIsNone(assigned[5])
            self.assertEqual(untested[102],'single_reference_exceeds_panel_budget')
            self.assertEqual(untested[999],'representative_missing')
            self.assertIn('descendant',untested[20])

    def test_shared_species_bin_selection_preserves_exact_pairs(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,summary,rows=sample_fixture(root)
            out=root/'out';out.mkdir()
            claims,reads,_,_=plant.extract_candidates(folder,rows,folder/'qc/read-assignments.tsv.gz',
                out,{102},False,10,3)
            self.assertEqual([c['species_taxid'] for c in claims],[102])
            self.assertEqual(claims[0]['original_ordinal'],'2')
            self.assertEqual(len(list(fragments(*reads))),1)

    def test_sam_near_tie_and_unavailable_mapq_do_not_support_species(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'test.sam'
            path.write_text(sam_record('1',1,mapq=255)+sam_record('1',2)+
                            sam_record('1',2,'T102|R1',secondary=True,score=198))
            contigs={f'T{t}|R1':dict(species_taxid=t,contig_class='nuclear',repeat_annotation=False) for t in (101,102)}
            align=plant.parse_sam(path,1,contigs,20,.95,.8)
            self.assertFalse(align[1,1]['passed']);self.assertTrue(align[1,2]['near_tied'])
            claim=dict(assigned_taxid=101,species_taxid=101,low_complexity=False)
            self.assertEqual(plant.category(claim,[align[1,1],align[1,2]],Tax()),'ambiguous')

    def test_duplicate_mapper_primaries_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'test.sam';path.write_text(sam_record('1',1)*2+sam_record('1',2))
            contigs={'T101|R1':dict(species_taxid=101,contig_class='nuclear',repeat_annotation=False)}
            with self.assertRaisesRegex(RuntimeError,'Duplicate primary'):
                plant.parse_sam(path,1,contigs,20,.95,.8)

    def test_higher_rank_support_stays_higher_rank(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'test.sam';path.write_text(sam_record('1',1)+sam_record('1',2))
            contigs={'T101|R1':dict(species_taxid=101,contig_class='nuclear',repeat_annotation=False)}
            align=plant.parse_sam(path,1,contigs,20,.95,.8)
            claim=dict(ordinal=1,assigned_taxid=10,species_taxid=0,low_complexity=False)
            row=plant.evidence_row(claim,[align[1,1],align[1,2]],Tax(),'001')
            summary=plant.aggregate_evidence([row],Tax())[0]
            self.assertEqual(row['category'],'higher_rank_compatible')
            self.assertEqual(summary['higher_rank_compatible_mates'],2)
            self.assertEqual(summary['exact_species_support_mates'],0)


class HumanTests(unittest.TestCase):
    def test_mapq_255_is_unavailable_in_both_alignment_engines(self):
        fields=sam_record('1',1,mapq=255).strip().split('\t')
        context=dict(contig_regex=None,min_mapq=20,min_aligned_fraction=.8,min_identity=.95)
        self.assertFalse(qc_human_screen._sam_positive(fields,context)[0])
        self.assertFalse(species_refine._sam_positive(fields,min_mapq=20,min_fraction=.8,min_identity=.95)[0])
        fields[2]='chrEBV';fields[4]='60'
        self.assertFalse(qc_human_screen._sam_positive(fields,context)[0])

    def test_candidates_report_duplicates_complexity_and_loci(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);reads=fastqs(root,3,base='human',low_last=True)
            calls=root/'calls.bin';calls.write_bytes(b'\3'*3)
            db=root/'evidence.sqlite'
            with sqlite3.connect(db) as conn:
                conn.execute('CREATE TABLE alignments (ordinal INTEGER,mate INTEGER,flag INTEGER,contig TEXT,'
                             'position INTEGER,cigar TEXT,mapq INTEGER,aligned_fraction REAL,identity REAL,passed INTEGER)')
                for ordinal in (3,1,2):
                    for mate in (2,1):
                        pos=(1 if mate==1 else 151)+(500 if ordinal==3 else 0)
                        conn.execute('INSERT INTO alignments VALUES (?,?,?,?,?,?,?,?,?,?)',
                                     (ordinal,mate,99 if mate==1 else 147,'chr1',pos,'100M',60,1,1,1))
            result=human.write_candidate_evidence([dict(r1=str(reads[0]),r2=str(reads[1]),bin='human')],
                calls,3,'any_mate',root/'evidence.tsv',db)[0]
            self.assertEqual(result['duplicate_sequence_fragments'],1)
            self.assertEqual(result['low_complexity_fragments'],1)
            self.assertEqual(result['passing_distinct_contigs'],1)
            self.assertEqual(result['passing_covered_reference_bases'],400)
            self.assertEqual(result['passing_nonoverlapping_regions'],4)
            self.assertEqual(result['both_mates_mapper_proper_pair'],3)


class EndToEndTests(unittest.TestCase):
    def test_zero_candidates_do_not_require_mapper_or_reference_panel(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);cfg=dict(work_dir=str(root/'work'),taxonomy_dir=str(root/'taxonomy'))
            with ExitStack() as stack:
                for name,replacement in [('require_compute',lambda x:None),('load_config',lambda *x:cfg),
                    ('code_identity',lambda:'test'),('configured_split_root',lambda c:root/'split'),
                    ('load_sample',lambda *x:(dict(paired=True,fragments=10),[],{}))]:
                    stack.enter_context(patch.object(plant,name,replacement))
                plant.main(['--sample','d/s'])
            report=json.loads((root/'work/qc/plants/d/s/plant-validation.json').read_text())
            self.assertEqual(report['status'],'no_candidates')
            self.assertEqual(report['candidate_fragments'],0)

    def run_fixture(self,root,fail_first=False,group='plant'):
        tids=(101,102,103) if group=='plant' else (9606,)
        folder,summary,rows=sample_fixture(root,tids,group)
        reps=rep_fixture(root,tids,group)
        if group=='human':reps.update(rep_fixture(root,[9598,9595],'mammal'))
        repfile=root/'representatives.tsv'
        with repfile.open('w',newline='') as out:
            writer=csv.DictWriter(out,fieldnames=list(next(iter(reps.values()))),delimiter='\t')
            writer.writeheader();writer.writerows(reps.values())
        cfg=dict(work_dir=str(root/'work'),taxonomy_dir=str(root/'taxonomy'),reference_raw=str(root/'refs'),tools_dir=str(root/'tools'))
        calls=[]
        def mapper(command,stdout,stderr,check):
            calls.append(command)
            if fail_first and len(calls)==1:raise subprocess.CalledProcessError(1,command)
            target=Path(command[-3]).read_text().splitlines()[0][1:]
            for n,_ in enumerate(fragments(command[-2],command[-1]),1):
                stdout.write(sam_record(str(n),1,target)+sam_record(str(n),2,target))
            return subprocess.CompletedProcess(command,0)
        with ExitStack() as stack:
            for name,replacement in [('require_compute',lambda x:None),('load_config',lambda *x:cfg),
                    ('code_identity',lambda:'test'),('configured_split_root',lambda c:root/'split'),
                    ('load_sample',lambda *x:(summary,rows,{})),('Taxonomy',lambda x:Tax()),
                    ('_executable',lambda x:'/bin/true')]:
                stack.enter_context(patch.object(plant,name,replacement))
            stack.enter_context(patch.object(plant.subprocess,'run',mapper))
            stack.enter_context(redirect_stdout(io.StringIO()))
            args=['--sample','d/s','--representatives',str(repfile),'--group',group,
                  '--max-species','2' if group=='plant' else '10','--max-reference-gib','.000001']
            if fail_first:
                with self.assertRaises(SystemExit):plant.main(args)
            else:plant.main(args)
        output=root/'work/qc'/('plants' if group=='plant' else 'human-competitive')/'d/s'
        report=json.loads((output/('plant-validation.json' if group=='plant' else 'human-competitive-validation.json')).read_text())
        return report,output,calls

    def test_serial_panels_complete_and_remove_large_intermediates(self):
        with tempfile.TemporaryDirectory() as d:
            report,out,calls=self.run_fixture(Path(d))
            self.assertEqual(report['tested_fragments'],3)
            self.assertEqual(report['categories'],{'supported_same_species':3})
            self.assertEqual(len(calls),3)
            self.assertFalse(list(out.glob('.panel.*')))
            self.assertTrue((out/'fragment-evidence.tsv').is_file())

    def test_panel_failure_retains_other_results_and_is_not_reusable_complete(self):
        with tempfile.TemporaryDirectory() as d:
            report,out,_=self.run_fixture(Path(d),fail_first=True)
            self.assertEqual(report['tested_fragments'],2)
            self.assertEqual(report['untested_fragments'],1)
            self.assertEqual(report['status'],'partially_tested')
            self.assertEqual(json.loads((out/'complete.json').read_text())['status'],'incomplete')

    def test_competitive_human_includes_primate_competitors(self):
        with tempfile.TemporaryDirectory() as d:
            report,out,_=self.run_fixture(Path(d),group='human')
            self.assertEqual(report['tested_fragments'],1)
            self.assertEqual(report['target_group'],'human')
            selected=json.loads((out/'panel.json').read_text())['panels'][0]['selected']
            self.assertTrue({9606,9598,9595}<={x['species_taxid'] for x in selected})


if __name__=='__main__':unittest.main(verbosity=2)
