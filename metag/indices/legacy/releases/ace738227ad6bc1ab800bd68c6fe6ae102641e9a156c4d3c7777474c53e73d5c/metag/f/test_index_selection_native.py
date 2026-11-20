"""Real native group relocation, automatic estimate reuse and atomic replacement."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from common import sha256
import group_library as group
import partitioned_build as union
import hierarchical_split as h
from index_selection import select_index
from reference import GIB
from test_group_library import fixture
from test_hierarchy_native import build_detailed_categories
from test_hierarchy_evaluate_native import normalize_fixture_category_padding,integration_available


@unittest.skipUnless(integration_available(),'Pinned native tools required')
class NamedIndexNativeTests(unittest.TestCase):
    def test_single_library_legacy_and_optimized_can_both_classify(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);cfg,catalog,sequences,_=fixture(d,os.environ['METAG_TEST_KRAKEN_BIN'])
            cfg.update(output_dir=str(root/'species'),taxonomy_dir=catalog['taxonomy_dir'],group_confidence=0.,group_minimum_hit_groups=2,split_gzip_level=3)
            # Original standalone library, then simulate the user's lib0 rename.
            old=group.build_group(cfg,max_gib=48*1024/GIB)
            build_detailed_categories(cfg,catalog,old);normalize_fixture_category_padding(cfg,catalog)
            with patch('partitioned_build._preparation',return_value=catalog):full=union.publish(cfg)
            base=Path(cfg['library_dir']);saved={str(p.relative_to(base)):sha256(p) for p in base.rglob('*') if p.is_file()}
            renamed=base.with_name('lib0');base.rename(renamed);renamed.rename(base)
            (base/'group').rename(base/'.old-group');(base/'group').mkdir();(base/'.old-group').rename(base/'group/legacy')
            legacy=select_index(cfg,'legacy');optimized=select_index(cfg)
            self.assertEqual(group.compact_snapshot(legacy)['identity'],old['identity'])
            self.assertEqual(union.database_snapshot(legacy)['identity'],full['identity'])
            calls=[];real=group._run_stream
            def run(command,*args):calls.append(str(command[0]));return real(command,*args)
            with patch.object(group,'_run_stream',side_effect=run):
                built=group.build_group(optimized,max_gib=96*1024/GIB,estimate='auto',replace=True)
            self.assertEqual([Path(x).name for x in calls],['build_db'])
            self.assertGreater(built['downsampling']['capacity_fraction'],old['downsampling']['capacity_fraction'])
            with patch.object(group,'_run_stream',side_effect=AssertionError('Unchanged index rebuilt')):
                again=group.build_group(optimized,max_gib=96*1024/GIB,estimate='auto',replace=True)
                self.assertEqual(again['identity'],built['identity'])
            # A failed replacement keeps the last published optimized index.
            with patch.object(group,'_run_stream',side_effect=RuntimeError('injected interruption')):
                with self.assertRaises(RuntimeError):group.build_group(optimized,max_gib=100*1024/GIB,estimate='auto',replace=True)
            self.assertEqual(group.compact_snapshot(optimized)['identity'],built['identity'])
            sample=dict(sample_id='d/s',sample='s',dataset='d')
            for mate in (1,2):
                p=root/f'R{mate}.fq';seq=sequences[9606][500:650]
                p.write_text(f'@r/{mate}\n{seq}\n+\n'+('I'*150)+'\n');sample[f'r{mate}']=str(p)
            for selected in (legacy,optimized):
                out=h.split_group(selected,sample)
                summary=json.loads((out/'group/qc/summary.json').read_text())
                self.assertEqual(summary['group_counts'],{'human':1})
            self.assertNotEqual(legacy['output_dir'],optimized['output_dir'])
            for relative,digest in saved.items():
                actual=base/(relative.replace('group/','group/legacy/',1) if relative.startswith('group/') else relative)
                self.assertEqual(sha256(actual),digest,relative)

    def test_auto_estimates_when_no_legacy_index_exists(self):
        with tempfile.TemporaryDirectory() as d:
            cfg,_,_,_=fixture(d,os.environ['METAG_TEST_KRAKEN_BIN']);cfg['output_dir']=str(Path(d)/'species')
            cfg=select_index(cfg);calls=[];real=group._run_stream
            def run(command,*args):calls.append(Path(command[0]).name);return real(command,*args)
            with patch.object(group,'_run_stream',side_effect=run):group.build_group(cfg,max_gib=48*1024/GIB,estimate='auto')
            self.assertEqual(calls,['estimate_capacity','build_db'])

if __name__=='__main__':unittest.main()
