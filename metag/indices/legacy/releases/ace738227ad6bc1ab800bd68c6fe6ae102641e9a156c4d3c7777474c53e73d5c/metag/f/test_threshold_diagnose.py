#!/usr/bin/env python3
"""Focused regression checks for whole-input threshold sensitivity."""
import csv
import gzip
import json
from pathlib import Path
import random
import struct
import tempfile
import types
import unittest
from unittest import mock

import threshold_diagnose as diagnostic
from fastq_stream import fragments, new_digest, update_digest


class TinyTaxonomy:
    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

    def validate_current(self):
        pass

    def describe(self, taxid):
        status, group, species = {0:("unclassified", "unresolved", 0),
             2:("higher_rank", "bacteria", 0), 9606:("assigned_species", "human", 9606),
             33090:("higher_rank", "plant", 0), 3517:("assigned_species", "plant", 3517)}[taxid]
        return dict(status=status, group=group, species_taxid=species)


def fixture(folder, count=37, paired=True):
    folder = Path(folder)
    paths = [folder / "raw.R1.fastq.gz"] + ([folder / "raw.R2.fastq.gz"] if paired else [])
    for mate,path in enumerate(paths,1):
        with gzip.open(path,"wb") as handle:
            for number in range(count):
                sequence = (b"ACGT" * 8) + bytes([b"ACGT"[number % 4]])
                handle.write(b"@read" + str(number).encode() + b"/" + str(mate).encode() + b" comment\r\n")
                handle.write(sequence + b"\r\n+annotated\r\n" + b"I"*len(sequence) + b"\r\n")
    sample = dict(sample_id="example/sample",dataset="example",sample="sample",r1=str(paths[0]),
                  r2=str(paths[1]) if paired else None)
    digest = new_digest()
    for a,b in fragments(sample['r1'],sample.get('r2')):
        update_digest(digest,a,b)
    return sample,digest.hexdigest()


def fake_query_job(command, **kwargs):
    payload = json.loads(Path(command[-1]).read_text())
    folder = Path(payload['out'])
    count = payload['expected_fragments']
    ids = [0,2,9606,3517,33090] * (count//5+1)
    ids = ids[:count]
    name = payload['profile'][0]
    if name == 'c0g2':
        ids[0] = 33090
    elif name == 'c0g1':
        ids[0] = 3517
    (folder/'taxids.bin').write_bytes(b''.join(struct.pack('<Q',x) for x in ids))
    (folder/'query.json').write_text(json.dumps(dict(fragments=count, input_digest=payload['expected_digest'],
          taxids=str(folder/'taxids.bin'), seconds=1.0, memory_mapping=True, batch_index_load_seconds=0.5)))
    return types.SimpleNamespace(returncode=0)


class ThresholdTests(unittest.TestCase):
    def test_seeded_selection_covers_full_input_and_preserves_bytes(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); sample,digest=fixture(root,101)
            result=diagnostic.extract_sample(sample,101,digest,13,1234,root/'sampled')
            expected=sorted(x+1 for x in random.Random(1234).sample(range(101),13))
            self.assertEqual(result['original_ordinals'],expected)
            self.assertEqual(len(set(expected)),13)
            self.assertGreater(max(expected),13)
            all_pairs=list(fragments(sample['r1'],sample['r2']))
            selected=list(fragments(result['r1'],result['r2']))
            self.assertEqual(selected,[all_pairs[i-1] for i in expected])
            self.assertTrue(result['full_input_digest_verified'])
            self.assertIn(b'\r\n',selected[0][0][0])

    def test_requested_size_clamps_to_population_and_single_end(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); sample,digest=fixture(root,7,False)
            result=diagnostic.extract_sample(sample,7,digest,100,23,root/'sampled')
            self.assertEqual(result['original_ordinals'],list(range(1,8)))
            self.assertIsNone(result['r2'])
            self.assertEqual(result['input_digest'],digest)

    def test_full_scan_rejects_changed_bytes_not_just_ids(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); sample,digest=fixture(root)
            with gzip.open(sample['r1'],'rb') as h:
                raw=h.read()
            with gzip.open(sample['r1'],'wb') as h:
                h.write(raw.replace(b'IIII',b'HHHH',1))
            with self.assertRaisesRegex(RuntimeError,'bytes/order/count'):
                diagnostic.extract_sample(sample,37,digest,4,8,root/'sampled')
            self.assertTrue((root/'sampled/sample.R1.fastq').exists())

    def test_full_scan_rejects_wrong_counts_and_mate_order(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); sample,digest=fixture(root,11)
            with self.assertRaisesRegex(RuntimeError,'more fragments'):
                diagnostic.extract_sample(sample,10,digest,3,8,root/'wrong_count')
            sample['r1'],sample['r2']=sample['r2'],sample['r1']
            with self.assertRaisesRegex(RuntimeError,'Reversed/duplicated'):
                diagnostic.extract_sample(sample,11,digest,3,8,root/'wrong_mates')

    def test_transition_conservation_and_target_counts(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            sampled=dict(fragments=5,original_ordinals=[1,5,12,23,34],input_digest='sample')
            records=[]
            arrays=[[0,2,9606,3517,33090],[2,33090,9606,3517,33090],[3517,3517,9606,3517,33090]]
            for (name,conf,g),ids in zip(diagnostic.PROFILES,arrays):
                path=root/(name+'.bin'); path.write_bytes(b''.join(struct.pack('<Q',x) for x in ids))
                records.append(dict(profile=name,confidence=conf,minimum_hit_groups=g,status='complete',error='',
                                    query=dict(taxids=str(path),seconds=1,memory_mapping=True)))
            summary=diagnostic.summarize(records,sampled,TinyTaxonomy(),root)
            self.assertEqual([r['plant_fragments'] for r in summary],[2,3,4])
            self.assertEqual([r['human_fragments'] for r in summary],[1,1,1])
            with open(root/'threshold-transitions.tsv') as h:
                rows=list(csv.DictReader(h,delimiter='\t'))
            pairs=set((r['left_profile'],r['right_profile']) for r in rows)
            self.assertEqual(len(pairs),3)
            for left,right in pairs:
                self.assertEqual(sum(int(r['fragments']) for r in rows if r['left_profile']==left and r['right_profile']==right),5)
            (root/'c0g1.bin').write_bytes(b'bad')
            diagnostic.summarize(records,sampled,TinyTaxonomy(),root)
            self.assertEqual(records[-1]['status'],'failed')

    def _main_fixture(self,root):
        sample,digest=fixture(root,17)
        source=root/'species/example/sample/qc';source.mkdir(parents=True)
        (source/'complete.json').write_text(json.dumps(dict(status='complete')))
        (source/'summary.json').write_text(json.dumps(dict(sample_id='example/sample',fragments=17,paired=True,
                  input_digest=digest,validated_scope='entire_input',requested_limit=0)))
        (source/'provenance.json').write_text(json.dumps(dict(context=dict(sample=sample,database_identity='union'))))
        cfg=dict(database_layout='partitioned',output_dir=str(root/'species'),work_dir=str(root/'work'),storage_layout='separated')
        return cfg

    def test_serial_workers_continue_failure_and_keep_sample(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);cfg=self._main_fixture(root);seen=[]
            def runner(command,**kwargs):
                job=json.loads(Path(command[-1]).read_text());seen.append(job)
                self.assertTrue(job['config']['union_memory_mapping'])
                if job['profile'][0]=='c0g2':return types.SimpleNamespace(returncode=7)
                return fake_query_job(command,**kwargs)
            with mock.patch.object(diagnostic,'load_config',return_value=cfg), \
                 mock.patch.object(diagnostic,'engine_snapshot',return_value=dict(identity='union')), \
                 mock.patch.object(diagnostic,'validate_snapshot'), \
                 mock.patch.object(diagnostic,'load_taxonomy',return_value=TinyTaxonomy()), \
                 mock.patch.object(diagnostic.subprocess,'run',side_effect=runner):
                with self.assertRaisesRegex(RuntimeError,'incomplete'):
                    diagnostic.main(['--sample','example/sample','--out',str(root/'out'),'--allow-local','--fragments','5'])
            self.assertEqual([r['profile'][0] for r in seen],['c10g2','c0g2','c0g1'])
            self.assertEqual(len({r['sample']['r1'] for r in seen}),1)
            self.assertEqual(len({r['expected_digest'] for r in seen}),1)
            self.assertTrue((root/'out/sample.R1.fastq').exists())
            self.assertFalse((root/'out/complete.json').exists())
            report=json.loads((root/'out/threshold-diagnostic.json').read_text())
            self.assertEqual(report['status'],'incomplete')

    def test_success_commits_and_releases_only_sample_fastqs(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);cfg=self._main_fixture(root)
            with mock.patch.object(diagnostic,'load_config',return_value=cfg), \
                 mock.patch.object(diagnostic,'engine_snapshot',return_value=dict(identity='union')), \
                 mock.patch.object(diagnostic,'validate_snapshot'), \
                 mock.patch.object(diagnostic,'load_taxonomy',return_value=TinyTaxonomy()), \
                 mock.patch.object(diagnostic.subprocess,'run',side_effect=fake_query_job):
                diagnostic.main(['--sample','example/sample','--out',str(root/'out'),'--allow-local','--fragments','5'])
            self.assertTrue((root/'out/complete.json').exists())
            self.assertFalse((root/'out/sample.R1.fastq').exists())
            self.assertTrue((root/'raw.R1.fastq.gz').exists())
            self.assertTrue((root/'out/sample-ordinals.tsv').exists())
            self.assertTrue((root/'out/c10g2/taxids.bin').exists())


if __name__=='__main__':
    unittest.main()
