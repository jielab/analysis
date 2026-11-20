"""Shared cache leases, moved library baselines and stage ordering."""
import copy
import json
import multiprocessing
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import alignment_cache as cache
from common import atomic_json, fingerprint
from library_paths import create_bundle, select_library, prepared_catalog
from test_group_library import fixture


def wait_for_reference(cfg,row,started,acquired):
    started.set()
    with cache.normalized_reference(cfg,row,1000):acquired.set()


class SharedCacheTests(unittest.TestCase):
    def test_old_device_number_does_not_invalidate_same_shared_files(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'file').write_text('data');cache.publish(root,'key',['file'])
            rec=json.loads((root/'complete.json').read_text());rec['files']['file']['device']=999999
            atomic_json(root/'complete.json',rec);self.assertIsNotNone(cache.valid_entry(root,'key'))
            (root/'file').write_text('changed');self.assertIsNone(cache.valid_entry(root,'key'))

    def test_reference_reader_keeps_lock_until_bytes_consumed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ref=root/'a.fa';ref.write_text('>x\nACGT\n')
            row=dict(species_taxid=1,resolved_path=str(ref),bases=4)
            cfg=dict(work_dir=str(root/'work'));ctx=multiprocessing.get_context('spawn')
            started,acquired=ctx.Event(),ctx.Event();process=None
            try:
                with cache.normalized_reference(cfg,row,1000) as (folder,_,_):
                    process=ctx.Process(target=wait_for_reference,args=(cfg,row,started,acquired));process.start()
                    self.assertTrue(started.wait(15))
                    self.assertFalse(acquired.wait(.3),'Second reader acquired key during active lease')
                    self.assertIn(b'ACGT',(folder/'reference.fa').read_bytes())
                self.assertTrue(acquired.wait(15));process.join(15);self.assertEqual(process.exitcode,0)
            finally:
                if process is not None and process.is_alive():process.terminate();process.join()


class PortableProjectCacheTests(unittest.TestCase):
    def test_node_local_device_does_not_change_project_or_taxonomy_key(self):
        from types import SimpleNamespace
        import project_cache
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'file';path.write_text('fixed');st=path.stat();values=[]
            for device in (49,46):
                observed=SimpleNamespace(st_size=st.st_size,st_mtime_ns=st.st_mtime_ns,
                    st_ctime_ns=st.st_ctime_ns,st_ino=st.st_ino,st_dev=device)
                with patch.object(Path,'stat',return_value=observed):values.append(project_cache._stamp(path))
            self.assertEqual(values[0],values[1]);self.assertNotIn('device',values[0])

    def test_unselected_fastq_is_not_opened_during_batch_discovery(self):
        import project_cache
        from query import discover_samples
        from test_workflow_cache import DatasetCacheTests
        with tempfile.TemporaryDirectory() as d:
            f=DatasetCacheTests();f.make_reads(d,'a','one');f.make_reads(d,'a','two')
            cfg=dict(reads_raw=str(Path(d)/'raw'),work_dir=str(Path(d)/'work'))
            with patch.object(project_cache,'_probe',side_effect=AssertionError('unselected FASTQ opened')):
                samples=discover_samples(cfg,datasets='a',probe=False)
            calls=[]
            with patch.object(project_cache,'_probe',side_effect=lambda sample:calls.append(sample['sample_id']) or dict(ok=True)):
                project_cache.probe_samples(cfg,samples[:1])
            self.assertEqual(calls,['a/one'])


class LibrarySelectionTests(unittest.TestCase):
    def test_relocated_prepared_catalog_and_bundle_keep_source_immutable(self):
        with tempfile.TemporaryDirectory() as d:
            cfg,catalog,_,_=fixture(d)
            original=Path(cfg['library_dir']);saved=(original/'catalog.json').read_bytes()
            full={'status':'complete','context':{'fixture':True}};full['identity']=fingerprint(full['context'])
            atomic_json(original/'index/database.json',full)
            moved=original.with_name('lib0');original.rename(moved)
            actual=prepared_catalog(moved);self.assertEqual(actual['identity'],catalog['identity'])
            for group in actual['groups'].values():self.assertTrue(all(Path(row['path']).is_file() for row in group['inputs']))
            cfg.update(library_dir=str(original),output_dir=str(Path(d)/'species'))
            create_bundle(original,moved);selected=select_library(cfg,str(original))
            self.assertEqual(selected['library_dir'],str(moved));self.assertEqual(selected['group_library_dir'],str(original/'group'))
            self.assertEqual((moved/'catalog.json').read_bytes(),saved)
            self.assertEqual(select_library(selected)['work_dir'],selected['work_dir'])
            changed=copy.deepcopy(full);changed['context']['fixture']=False;changed['identity']=fingerprint(changed['context'])
            atomic_json(moved/'index/database.json',changed)
            with self.assertRaisesRegex(RuntimeError,'full index changed'):select_library(cfg,str(original))

    def test_lib0_outputs_and_work_do_not_overlap_new_library(self):
        cfg=dict(library_dir='/data/annot/metag/lib',output_dir='/data/metag/species',work_dir='/analysis/metag')
        a,b=select_library(cfg,'lib'),select_library(cfg,'lib0')
        self.assertNotEqual(a['output_dir'],b['output_dir']);self.assertNotEqual(a['work_dir'],b['work_dir'])
        with self.assertRaises(ValueError):select_library(cfg,'../raw')


class BatchOrderTests(unittest.TestCase):
    def test_all_samples_finish_split_before_any_evaluation(self):
        import hierarchical_split as h
        with tempfile.TemporaryDirectory() as d:
            cfg=dict(output_dir=d,work_dir=d);samples=[dict(sample_id='d/a'),dict(sample_id='d/b')];events=[]
            def split(kind):
                def run(cfg,sample,**kw):
                    events.append((kind,sample['sample_id']))
                    qc=h.qc_folder(Path(d)/sample['sample_id'],kind);qc.mkdir(parents=True,exist_ok=True);atomic_json(qc/'summary.json',{})
                return run
            def evaluate(root,stage,cfg,**kw):events.append(('evaluate',root.name));return {'status':'complete'}
            with patch.object(h,'load_config',return_value=cfg),patch.object(h,'require_compute'),patch.object(h,'discover_samples',return_value=samples),patch.object(h,'split_group',side_effect=split('group')),patch.object(h,'split_species',side_effect=split('species')),patch.object(h,'evaluate',side_effect=evaluate),patch.object(h,'report_samples'):
                h.main(['all'])
            self.assertEqual(events[:4],[('group','d/a'),('group','d/b'),('species','d/a'),('species','d/b')]);self.assertEqual(len(events),8)

if __name__=='__main__':unittest.main()
