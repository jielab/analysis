"""One library, independent named indexes and explicit public defaults."""
import argparse
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from common import atomic_json, fingerprint
from index_selection import select_index
import library_cli
from qc_submit import frozen_args, pending_hierarchy_jobs
from test_group_library import fixture


class IndexSelectionTests(unittest.TestCase):
    def test_names_share_inputs_and_keep_outputs_and_dependencies_separate(self):
        cfg=dict(library_dir='/references/lib',work_dir='/analysis/metag',output_dir='/reads/species')
        old=select_index(cfg,'legacy');new=select_index(cfg)
        self.assertEqual(old['library_dir'],new['library_dir'])
        self.assertEqual(old['db_dir'],new['db_dir'])
        self.assertEqual(new['group_library_dir'],'/references/lib/group/optimized')
        self.assertEqual(old['group_library_dir'],'/references/lib/group/legacy')
        self.assertEqual(new['output_dir'],'/reads/species/optimized')
        self.assertEqual(old['output_dir'],'/reads/species/legacy')
        self.assertNotEqual(old['work_dir'],new['work_dir'])
        self.assertEqual(select_index(new),new)
        self.assertEqual(select_index(old,'optimized'),new)
        self.assertEqual(select_index(cfg,'full')['group_routing'],'full')
        for name in ('128','384','../raw','lib0'):
            with self.assertRaises(ValueError):select_index(cfg,name)

    def test_default_build_is_one_optimized_job_in_same_library(self):
        with tempfile.TemporaryDirectory() as d:
            cfg,catalog,_,_=fixture(d);cfg['output_dir']=str(Path(d)/'species')
            context=dict(preparation_identity=catalog['identity'])
            atomic_json(Path(cfg['library_dir'])/'index/database.json',dict(status='complete',context=context,identity=fingerprint(context)))
            with patch.object(library_cli,'load_config',return_value=cfg),patch('partitioned_build.database_snapshot'),patch('group_library.source_catalog'),patch.object(library_cli.subprocess,'run') as run:
                library_cli.main(['build','--plan-only'])
            run.assert_called_once();cmd=run.call_args.args[0]
            self.assertEqual(cmd[cmd.index('--index-use')+1],'optimized')
            self.assertEqual(cmd[cmd.index('--estimate')+1],'auto')
            self.assertEqual(cmd[cmd.index('--replace')+1],'TRUE')
            self.assertNotIn('--library',cmd);self.assertFalse((Path(cfg['library_dir'])/'library.json').exists())
            with patch.object(library_cli,'load_config',return_value=cfg),patch('partitioned_build.database_snapshot'),patch('group_library.source_catalog'),patch('group_library.compact_snapshot') as check,patch.object(library_cli.subprocess,'run') as run:
                library_cli.main(['build','--index-use','legacy'])
                check.assert_called_once();run.assert_not_called()

    def test_frozen_worker_uses_resolved_config_not_a_later_default(self):
        args=frozen_args(['all','--index-use','legacy','--config','old.json','--routing','compact'],Path('/frozen/config.json'),Path('/work'))
        self.assertNotIn('--index-use',args);self.assertIn('/frozen/config.json',args)

    def test_build_dependency_tracks_selected_index(self):
        with tempfile.TemporaryDirectory() as d:
            cfg=dict(library_dir=str(Path(d)/'lib'),work_dir=str(Path(d)/'work'),output_dir=str(Path(d)/'species'))
            chosen=select_index(cfg)
            atomic_json(Path(chosen['work_dir'])/'jobs/qc/group_library/run/submission.json',dict(status='submitted',job_id='123',scope=dict(group_library_dir=chosen['group_library_dir'])))
            scope=dict(module='hierarchical_split',action='group,species,evaluate,report',group_library_dir=chosen['group_library_dir'])
            with patch('qc_submit.active_job',return_value=True):
                self.assertEqual(pending_hierarchy_jobs(chosen,scope),['123'])
                legacy=select_index(cfg,'legacy');scope['group_library_dir']=legacy['group_library_dir']
                self.assertEqual(pending_hierarchy_jobs(legacy,scope),[])

if __name__=='__main__':unittest.main()
