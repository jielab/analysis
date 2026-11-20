"""Cache hits must preserve competitive results and reject stale inputs."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import alignment_cache as cache
from plant_validate import build_panel
from stage_alignment import evaluate_alignment
import test_stage_alignment as stage_tests
from test_stage_alignment import minimap_binary


class ReferenceCacheTests(unittest.TestCase):
    def test_reference_cache_preserves_bytes_and_invalidates(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);cfg={'work_dir':str(root/'work')}
            ref=root/'a.fa';ref.write_text('>chr1 chromosome nuclear\nACGTNacgt\n>cp chloroplast\nTTTAAA\n')
            rows=[dict(species_taxid=101,resolved_path=str(ref),bases=15)]
            expected,actual=build_panel(rows,root/'old.fa',1000)
            got,n,metrics=cache.build_cached_panel(cfg,rows,root/'new.fa',1000)
            self.assertEqual((got,n),(expected,actual));self.assertEqual((root/'old.fa').read_bytes(),(root/'new.fa').read_bytes())
            with patch.object(cache,'build_panel',side_effect=AssertionError('unexpected decode')):
                _,_,metrics=cache.build_cached_panel(cfg,rows,root/'hit.fa',1000)
            self.assertEqual(metrics['reference_hits'],1)
            ref.write_text(ref.read_text().replace('TTTAAA','CCCGGG'))
            _,_,metrics=cache.build_cached_panel(cfg,rows,root/'changed.fa',1000)
            self.assertEqual(metrics['reference_builds'],1)
            self.assertIn(b'CCCGGG',(root/'changed.fa').read_bytes())
            with self.assertRaises(RuntimeError):cache.build_cached_panel(cfg,rows,root/'small.fa',2)

    def test_zero_budget_falls_back_and_corruption_rebuilds(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ref=root/'a.fa';ref.write_text('>x\nACGT\n')
            rows=[dict(species_taxid=1,resolved_path=str(ref),bases=4)]
            cfg=dict(work_dir=str(root/'work'),stage_alignment={'reference_cache_gib':0})
            _,_,m=cache.build_cached_panel(cfg,rows,root/'zero.fa',100)
            self.assertEqual(m['reference_fallbacks'],1)
            cfg['stage_alignment']['reference_cache_gib']=1
            cache.build_cached_panel(cfg,rows,root/'one.fa',100)
            entry=next((root/'work/cache/alignment/references').glob('*/reference.fa'));entry.write_text('broken')
            cache.build_cached_panel(cfg,rows,root/'fixed.fa',100)
            self.assertEqual((root/'zero.fa').read_bytes(),(root/'fixed.fa').read_bytes())


class NativeCacheTests(unittest.TestCase):
    def setUp(self):
        if not minimap_binary():self.skipTest('minimap2 unavailable')
        self.fixture=stage_tests.StageAlignmentTests();self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)

    def test_native_hits_checkpoints_and_invalidation(self):
        f=self.fixture
        cold=copy.deepcopy(f.cfg);cold['stage_alignment']['cache_enabled']=False
        a=evaluate_alignment(f.sample,'group',cold)
        self.assertEqual(a['status'],'complete',a)
        out=f.qc/'alignment';before=(out/'sampled-evidence.tsv').read_text()
        (out/'summary.json').unlink()
        b=evaluate_alignment(f.sample,'group',f.cfg)
        self.assertEqual(b['status'],'complete',b)
        self.assertEqual(before,(out/'sampled-evidence.tsv').read_text())
        # Avoid whole-result hits and force a fresh mapping to an existing MMI.
        for p in (out/'checkpoints').glob('*.json'):p.unlink()
        summary=json.loads((out/'summary.json').read_text());summary['status']='running';(out/'summary.json').write_text(json.dumps(summary))
        c=evaluate_alignment(f.sample,'group',f.cfg)
        self.assertEqual(c['status'],'complete',c)
        self.assertTrue(all(x['cache']['index_hit'] for x in c['commands']))
        self.assertEqual(before,(out/'sampled-evidence.tsv').read_text())
        summary=json.loads((out/'summary.json').read_text());summary['status']='running';(out/'summary.json').write_text(json.dumps(summary))
        with patch('stage_alignment.panel_index',side_effect=AssertionError('unexpected mapping')):
            d=evaluate_alignment(f.sample,'group',f.cfg)
        self.assertEqual(d['status'],'complete',d);self.assertEqual(d['cached_panels'],d['planned_panels'])
        with patch('stage_alignment.panel_index',side_effect=AssertionError('unexpected mapping')):
            e=evaluate_alignment(f.sample,'group',f.cfg)
        self.assertTrue(e.get('cache_reused'),e)
        f.cfg['stage_alignment']['min_mapq']=30
        changed=evaluate_alignment(f.sample,'group',f.cfg)
        self.assertEqual(changed['status'],'complete',changed)
        self.assertFalse(changed.get('cache_reused',False));self.assertEqual(changed['cached_panels'],0)

    def test_different_reference_order_and_mapper_state_change_key(self):
        f=self.fixture;root=f.base
        rows=[dict(species_taxid=tid,resolved_path=str(root/'references'/f'{tid}.fa'),bases=9000) for tid in (101,102)]
        exe=minimap_binary()
        with cache.panel_index(f.cfg,rows,exe,'test',2,100000,root/'panel.fa') as (_,_,_,m):
            first=m['index_key']
        with cache.panel_index(f.cfg,list(reversed(rows)),exe,'test',2,100000,root/'panel.fa') as (_,_,_,m):
            self.assertNotEqual(m['index_key'],first)
        with cache.panel_index(f.cfg,rows,exe,'different-version',2,100000,root/'panel.fa') as (_,_,_,m):
            self.assertNotEqual(m['index_key'],first)


def _concurrent_reference(args):
    cfg,row,out=args
    return cache.build_cached_panel(cfg,[row],out,1000)[2]


class ConcurrencyTests(unittest.TestCase):
    def test_shared_reference_is_built_once(self):
        import multiprocessing
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ref=root/'a.fa';ref.write_text('>x\nACGT\n')
            row=dict(species_taxid=1,resolved_path=str(ref),bases=4)
            cfg=dict(work_dir=str(root/'work'))
            with multiprocessing.get_context('spawn').Pool(2) as pool:
                results=pool.map(_concurrent_reference,[(cfg,row,str(root/f'out{i}')) for i in range(2)])
            self.assertEqual(sum(r['reference_builds'] for r in results),1)
            self.assertEqual(sum(r['reference_hits'] for r in results),1)


class CompatibilityTests(unittest.TestCase):
    def test_only_audited_code_change_can_reuse_classification(self):
        import cache_compatibility as c
        old=dict(code='old',sample='x',confidence=.1)
        new=dict(old,code='new')
        receipt=json.dumps(dict(current_code='new',compatible_codes=['old']))
        with patch.object(Path,'read_text',return_value=receipt):
            self.assertTrue(c.compatible_request(old,new))
            self.assertFalse(c.compatible_request(old,dict(new,confidence=.2)))
            self.assertFalse(c.compatible_request(old,dict(new,code='unverified')))
            self.assertFalse(c.compatible_request(dict(old,code='other'),new))

    def test_agreement_cache_reuses_success_and_rejects_changed_sample(self):
        import hierarchical_split as h
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);qc=root/'group/qc';qc.mkdir(parents=True)
            for name in ('sample.tsv','sample.R1.fastq.gz','sample.R2.fastq.gz'):(qc/name).write_text('fixture')
            request=dict(confidence=.1,minimum_hit_groups=2,minimum_base_quality=0)
            (qc/'summary.json').write_text(json.dumps(dict(request=request)))
            def compute(*args):
                out=qc/'agreement';out.mkdir(exist_ok=True)
                result=dict(status='complete')
                (out/'summary.json').write_text(json.dumps(result))
                for name in ('pairs.tsv','transitions.tsv'):(out/name).write_text('fixture')
                return result
            with patch('partitioned_build.database_snapshot',return_value={'identity':'db'}),patch.object(h,'_classifier_agreement_uncached',side_effect=compute) as run:
                h.classifier_agreement(root,'group',{})
                h.classifier_agreement(root,'group',{})
                self.assertEqual(run.call_count,1)
                (qc/'sample.tsv').write_text('new fixture')
                h.classifier_agreement(root,'group',{})
                self.assertEqual(run.call_count,2)
