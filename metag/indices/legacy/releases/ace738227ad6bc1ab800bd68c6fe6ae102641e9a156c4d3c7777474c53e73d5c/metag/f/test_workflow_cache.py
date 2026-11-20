import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import project_cache as cache
from query import discover_samples
from test_group_library import fixture
from taxonomy import Taxonomy


class DatasetCacheTests(unittest.TestCase):
    def make_reads(self,root,dataset,name,length=20):
        folder=Path(root)/'raw'/dataset
        folder.mkdir(parents=True,exist_ok=True)
        for mate in (1,2):
            with gzip.open(folder/f'{name}_{mate}.fastq.gz','wb') as f:
                f.write(f'@{name}/{mate}\n'.encode()+b'A'*length+b'\n+\n'+b'I'*length+b'\n')

    def test_reuse_changes_new_files_and_batch_size(self):
        with tempfile.TemporaryDirectory() as temp:
            cfg=dict(reads_raw=str(Path(temp)/'raw'),work_dir=str(Path(temp)/'work'))
            self.make_reads(temp,'a','s1')
            samples=discover_samples(cfg,datasets='a')
            first=cache.batch_plan(cfg,samples,10,'all')
            with patch.object(cache,'_probe',side_effect=AssertionError('unexpected FASTQ reread')):
                again=discover_samples(cfg,datasets='a')
                self.assertTrue(cfg['project_cache_discovery']['a']['inventory_hit'])
                self.assertEqual(cfg['project_cache_discovery']['a']['probe_hits'],1)
                self.assertEqual(cache.batch_plan(cfg,again,10,'all'),first)
            folder=Path(temp)/'work/cache/projects/a/batches'
            self.assertEqual(len(list(folder.glob('*.json'))),1)
            cache.batch_plan(cfg,again,1,'all')
            self.assertEqual(len(list(folder.glob('*.json'))),2)
            self.make_reads(temp,'a','s1',30)
            again=discover_samples(cfg,datasets='a')
            self.assertEqual(cfg['project_cache_discovery']['a']['probed_samples'],1)
            cache.batch_plan(cfg,again,10,'all')
            self.assertEqual(len(list(folder.glob('*.json'))),3)
            self.make_reads(temp,'a','s2')
            again=discover_samples(cfg,datasets='a')
            self.assertEqual(len(again),2)
            self.assertEqual(cfg['project_cache_discovery']['a']['probe_hits'],1)
            self.assertEqual(len(cache.batch_plan(cfg,again,1,'all')),2)

    def test_dataset_isolation_parameter_invalidation_and_corrupt_cache(self):
        with tempfile.TemporaryDirectory() as temp:
            cfg=dict(reads_raw=str(Path(temp)/'raw'),work_dir=str(Path(temp)/'work'))
            for d in ('a','b'):self.make_reads(temp,d,'s1')
            samples=discover_samples(cfg,datasets='a,b')
            cache.batch_plan(cfg,samples,10,'all')
            root=Path(temp)/'work/cache/projects'
            self.assertTrue((root/'a/batches').is_dir())
            self.assertTrue((root/'b/batches').is_dir())
            p=next((root/'a/batches').glob('*.json'))
            saved=json.loads(p.read_text());saved['batches']=[];p.write_text(json.dumps(saved))
            self.assertEqual(len(cache.batch_plan(cfg,samples,10,'all')),2)
            cfg['group_confidence']=0.2
            cache.batch_plan(cfg,samples,10,'all')
            self.assertEqual(len(list((root/'a/batches').glob('*.json'))),2)

    def test_taxonomy_reuse_matches_uncached_and_invalidates(self):
        with tempfile.TemporaryDirectory() as temp:
            cfg,catalog,_,_=fixture(temp)
            folder=catalog['taxonomy_dir']
            expected=Taxonomy(folder)
            first=cache.routing_taxonomy(cfg,folder)
            for tid in (0,9606,3702,562):
                original=expected.describe(tid)
                self.assertEqual({k:first.describe(tid)[k] for k in original},original)
            self.assertIs(cache.routing_taxonomy(cfg,folder),first)
            cache.close_routing_taxonomies()
            with patch.object(cache,'_compile',side_effect=AssertionError('unexpected rebuild')):
                reused=cache.routing_taxonomy(cfg,folder)
                self.assertTrue(reused.preparation['taxonomy_cache_hit'])
            names=Path(folder)/'names.dmp'
            names.write_text(names.read_text().replace('Homo sapiens','Human changed'))
            changed=cache.routing_taxonomy(cfg,folder)
            self.assertIsNot(changed,reused)
            self.assertFalse(changed.preparation['taxonomy_cache_hit'])
            self.assertEqual(changed.describe(9606)['name'],'Human changed')
            cache.close_routing_taxonomies()


if __name__=='__main__':unittest.main()
