"""End-to-end native relocation, estimate reuse and full-routing regression."""
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from common import sha256
import group_library as group
import hierarchical_split as h
import partitioned_build as union
from library_paths import create_bundle,select_library
from reference import GIB
from test_group_library import fixture
from test_hierarchy_native import build_detailed_categories,read_rows
from test_hierarchy_evaluate_native import normalize_fixture_category_padding,install_test_union,integration_available


@unittest.skipUnless(integration_available(),'Native Kraken/source/minimap2 required')
class NativeBundleTests(unittest.TestCase):
    def test_moved_baselines_new_group_estimate_and_full_routing(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);cfg,catalog,sequences,_=fixture(d,os.environ['METAG_TEST_KRAKEN_BIN'])
            cfg.update(output_dir=str(root/'species'),species_groups=['human','plant'],
                       group_confidence=0.,group_minimum_hit_groups=2,evaluation_profile='basic',
                       evaluation_anchor_references=1,evaluation_anchor_bases=8000,evaluation_anchor_pairs=4,
                       taxonomy_dir=catalog['taxonomy_dir'],split_gzip_level=3)
            compact=group.build_group(cfg,max_gib=48*1024/GIB)
            build_detailed_categories(cfg,catalog,compact);normalize_fixture_category_padding(cfg,catalog)
            with patch('partitioned_build._preparation',return_value=catalog):full=union.publish(cfg)
            install_test_union(cfg,os.environ['METAG_TEST_KRAKEN_SOURCE'],root/'union-source')
            original=Path(cfg['library_dir']);before={str(p.relative_to(original)):sha256(p) for p in original.rglob('*') if p.is_file()}
            moved=root/'lib0';original.rename(moved)
            old=select_library(cfg,'lib0')
            snapshot=union.database_snapshot(old,full=True)
            self.assertEqual(snapshot['identity'],full['identity'])
            self.assertEqual(group.compact_snapshot(old)['identity'],compact['identity'])
            self.assertTrue(all(Path(x).is_file() for x in Path(snapshot['_hashes_path']).read_text().splitlines()))
            forged=copy.deepcopy(snapshot);forged['taxonomy_dir']=str(root/'missing')
            with self.assertRaises(RuntimeError):union.validate_snapshot(forged)
            create_bundle(original,moved);new=select_library(cfg,'lib')
            real_stream=group._run_stream;calls=[]
            def stream(cmd,*args):calls.append(str(cmd[0]));return real_stream(cmd,*args)
            with patch.object(group,'_run_stream',side_effect=stream):
                built=group.build_group(new,max_gib=96*1024/GIB,estimate='reuse')
            self.assertEqual(len(calls),1);self.assertTrue(calls[0].endswith('build_db'))
            self.assertGreater(built['downsampling']['native_options']['selected_hash_fraction'],compact['downsampling']['native_options']['selected_hash_fraction'])
            sample=dict(sample_id='d/native',sample='native',dataset='d',r1='',r2='')
            for mate in (1,2):
                p=root/f'R{mate}.fastq';seq=sequences[9606][500:650]
                p.write_text(f'@r/{mate}\n{seq}\n+\n'+('I'*len(seq))+'\n');sample[f'r{mate}']=str(p)
            # Full union routing must use relocated hash paths, not stale /lib paths.
            new['group_routing']='full';result=h.split_group(new,sample)
            summary=json.loads((result/'group/qc/summary.json').read_text())
            self.assertEqual(summary['group_counts'],{'human':1})
            with patch('hierarchical_split.classifier_agreement',side_effect=AssertionError('basic must not scan union')),patch('stage_alignment.evaluate_alignment',side_effect=AssertionError('basic must not align')):
                evaluated=h.evaluate(result,'group',new)
            self.assertNotEqual(evaluated['status'],'failed',evaluated)
            self.assertEqual(evaluated['alignment']['status'],'skipped')
            self.assertEqual(before,{str(p.relative_to(moved)):sha256(p) for p in moved.rglob('*') if p.is_file()})
            target=moved/'index/plants/hash.k2d';data=target.read_bytes();target.write_bytes(data[:-1])
            with self.assertRaises(RuntimeError):union.database_snapshot(old)

if __name__=='__main__':unittest.main()
