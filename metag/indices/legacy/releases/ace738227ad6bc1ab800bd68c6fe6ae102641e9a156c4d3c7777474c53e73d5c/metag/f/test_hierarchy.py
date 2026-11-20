"""Focused regressions for hierarchy publication, sampling and job dependencies."""
import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from common import atomic_json, fingerprint
import hierarchical_split as h
import qc_submit as scheduler


class HierarchyTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()

    def publication_fixture(self, foreign=False):
        parent=self.root/'group';parent.mkdir()
        for name in ('human','plant'):(parent/name).mkdir()
        (parent/'human/species').mkdir();(parent/'human/species/old').write_text('old-human')
        if foreign:
            (parent/'plant/species').mkdir();(parent/'plant/species/foreign').write_text('foreign')
        qc=parent/'species.qc';qc.mkdir();atomic_json(qc/'summary.json',{'request':{'version':1},'processed_groups':['human']})
        stage=parent/'.new';stage.mkdir()
        for name in ('human','plant'):
            (stage/name/'species').mkdir(parents=True);(stage/name/'species/new').write_text('new-'+name)
        (stage/'species.qc').mkdir();atomic_json(stage/'species.qc/summary.json',{'request':{'version':2}})
        return parent,stage

    def test_foreign_child_preflight_preserves_every_old_output(self):
        parent,stage=self.publication_fixture(foreign=True)
        with self.assertRaisesRegex(RuntimeError,'unowned'):
            h.publish_species(stage,parent,['human','plant'],{'human'},{'version':2},True)
        self.assertEqual((parent/'human/species/old').read_text(),'old-human')
        self.assertEqual((parent/'plant/species/foreign').read_text(),'foreign')
        self.assertEqual(h.read_json(parent/'species.qc/summary.json')['request'],{'version':1})
        self.assertFalse(list(parent.glob('.species-backup-*')))

    def test_partial_rename_failure_rolls_back_children_and_receipt(self):
        parent,stage=self.publication_fixture();rename=Path.rename
        def injected(path,target):
            if path==stage/'plant/species':raise OSError('injected second-child failure')
            return rename(path,target)
        with patch.object(Path,'rename',injected):
            with self.assertRaisesRegex(OSError,'injected'):
                h.publish_species(stage,parent,['human','plant'],{'human'},{'version':2},True)
        self.assertEqual((parent/'human/species/old').read_text(),'old-human')
        self.assertFalse((parent/'human/species/new').exists())
        self.assertFalse((parent/'plant/species').exists())
        self.assertEqual(h.read_json(parent/'species.qc/summary.json')['request'],{'version':1})

    def test_interrupted_publication_recovery(self):
        parent,stage=self.publication_fixture();backup=parent/'.species-backup-interrupted';backup.mkdir()
        atomic_json(backup/'transaction.json',dict(schema=h.SCHEMA,prior=['human'],chosen=['human','plant'],
                    request_identity=fingerprint({'version':2})))
        (parent/'species.qc').rename(backup/'species.qc')
        (parent/'human/species').rename(backup/'human')
        (stage/'human/species').rename(parent/'human/species')
        h.recover_species_publication(parent)
        self.assertEqual((parent/'human/species/old').read_text(),'old-human')
        self.assertTrue((parent/'species.qc/summary.json').is_file())
        self.assertFalse(backup.exists())

    def test_successful_publication_keeps_coarse_inputs(self):
        parent,stage=self.publication_fixture();(parent/'human/R1.fastq.gz').write_bytes(b'coarse')
        h.publish_species(stage,parent,['human','plant'],{'human'},{'version':2},True)
        self.assertEqual((parent/'human/R1.fastq.gz').read_bytes(),b'coarse')
        self.assertEqual((parent/'plant/species/new').read_text(),'new-plant')
        self.assertFalse(list(parent.glob('.species-backup-*')))

    def test_reservoir_is_reproducible_not_a_prefix_and_records_population(self):
        cfg={'evaluation_pairs_per_stratum':25,'evaluation_target_pairs_per_stratum':25,'evaluation_seed':23}
        chosen=[]
        for version in (1,2):
            r=h.Reservoir(cfg)
            for n in range(1,1001):
                row=dict(zip(h.FIELDS,[n,'id','id/1','id/2','bacteria',562,562,'assigned_species','group/bacteria','',562]))
                a=(b'@id/1\n',b'ACGT\n',b'+\n',b'IIII\n');b=(b'@id/2\n',b'ACGT\n',b'+\n',b'IIII\n')
                r.add(row,a,b)
            out=self.root/str(version);r.publish(out)
            with (out/'sample.tsv').open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
            chosen.append([int(x['ordinal']) for x in rows])
            self.assertEqual({int(x['stratum_pairs']) for x in rows},{1000})
            self.assertEqual({int(x['sampled_pairs']) for x in rows},{25})
        self.assertEqual(chosen[0],chosen[1]);self.assertGreater(max(chosen[0]),800)

    def test_quality_parameter_changes_request_identity(self):
        sample={'sample_id':'d/s'}
        for mate in (1,2):
            p=self.root/f'{mate}.fq';p.write_text('x');sample[f'r{mate}']=str(p)
        with patch('hierarchical_split.code_identity',return_value='same-code'):
            a=h.stage_request({'minimum_base_quality':0},sample,'group',{})
            b=h.stage_request({'minimum_base_quality':20},sample,'group',{})
        self.assertNotEqual(a,b)

    def test_uncertain_refinement_can_promote_target_but_known_conflict_cannot(self):
        info={'taxid':9606,'species_taxid':9606,'group':'human'}
        rescued=h.species_route(info,'unclassified',0,{},20)
        self.assertEqual((rescued['group'],rescued['species_taxid'],rescued['leaf']),('human',9606,'9606'))
        conflict=h.species_route(info,'plant',33090,{},20)
        self.assertEqual((conflict['group'],conflict['species_taxid'],conflict['status']),('ambiguous',0,'group_conflict'))

    def test_dependency_scope_and_latest_parent(self):
        cfg={'work_dir':str(self.root)};parent=self.root/'jobs/qc/hierarchical_split'
        scope={'module':'hierarchical_split','action':'species','samples':['england/A'],'datasets':[],
               'from_existing':False,'output_dir':'/output','group_library_dir':'/lib/group'}
        for run,sample,job in [('run','england/B','101'),('run2','england/A','102')]:
            atomic_json(parent/run/'submission.json',dict(status='submitted',job_id=job,
                        scope={**scope,'action':'group','samples':[sample]}))
        with patch('qc_submit.active_job',return_value=True):
            self.assertEqual(scheduler.pending_hierarchy_jobs(cfg,scope),['102'])
        self.assertFalse(scheduler.scope_overlap(scope,{'samples':['england/B']}))

    def test_replacement_library_build_is_a_dependency_even_if_old_index_exists(self):
        lib=self.root/'lib';lib.mkdir();(lib/'database.json').write_text('{}')
        cfg={'work_dir':str(self.root)}
        scope={'module':'hierarchical_split','action':'group','samples':['d/s'],'datasets':[],
               'from_existing':False,'output_dir':'/output','group_library_dir':str(lib)}
        atomic_json(self.root/'jobs/qc/group_library/run/submission.json',dict(status='submitted',job_id='123',scope=scope))
        with patch('qc_submit.active_job',return_value=True):
            self.assertEqual(scheduler.pending_hierarchy_jobs(cfg,scope),['123'])


if __name__=='__main__':unittest.main()
