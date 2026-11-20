import argparse
import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import hierarchical_split as h
import metag_submit
import qc_submit
from workflow_selection import actions,batches,select_batches


class WorkflowTests(unittest.TestCase):
    def test_alias_order_and_rejections(self):
        self.assertEqual(actions('all'),actions('split-group,split-species,evaluate,report'))
        self.assertEqual(actions('split'),['group','species','report'])
        for text in ('species,group','group,group','all,report'):
            with self.assertRaises(argparse.ArgumentTypeError):actions(text)

    def test_dataset_batches_and_tail(self):
        samples=[dict(dataset=d,sample_id=f'{d}/{n}') for d in ('a','b') for n in range(5)]
        jobs=batches(samples,2,'2-')
        self.assertEqual([(d,n,len(s)) for d,n,s in jobs],[('a',2,2),('a',3,1),('b',2,2),('b',3,1)])
        for spec in ('0','4','2-1','1,1','1-2,2'):
            with self.assertRaises(ValueError):select_batches(spec,3)

    def run_workflow(self, action, fail_group=False):
        events=[]
        with tempfile.TemporaryDirectory() as temp:
            cfg=dict(output_dir=temp,work_dir=temp)
            sample=dict(sample_id='d/s',dataset='d')
            root=Path(temp)/'d/s'
            def split_group(*args,**kwargs):
                events.append('group')
                if fail_group:raise RuntimeError('fixture classification failure')
                p=root/'group/qc';p.mkdir(parents=True);(p/'summary.json').write_text('{}')
            def split_species(*args,**kwargs):
                events.append('species')
                p=root/'group/species.qc';p.mkdir();(p/'summary.json').write_text('{}')
            def evaluate(root,stage,cfg,**kwargs):
                events.append('evaluate:'+stage);return dict(status='complete')
            with patch.object(h,'load_config',return_value=cfg),patch.object(h,'require_compute'), \
                 patch.object(h,'discover_samples',return_value=[sample]), \
                 patch.object(h,'split_group',side_effect=split_group), \
                 patch.object(h,'split_species',side_effect=split_species), \
                 patch.object(h,'evaluate',side_effect=evaluate), \
                 patch.object(h,'report_samples',side_effect=lambda *a:events.append('report')):
                if fail_group:
                    with self.assertRaises(SystemExit):h.main([action])
                else:h.main([action])
        return events

    def test_full_sequence_evaluates_each_stage_once(self):
        self.assertEqual(self.run_workflow('all'),['group','report','species','report','evaluate:group','evaluate:species','report'])

    def test_single_split_has_no_hidden_evaluation(self):
        self.assertEqual(self.run_workflow('group'),['group','report'])

    def test_failed_group_prevents_species_and_evaluation(self):
        self.assertEqual(self.run_workflow('all',True),['group','report'])

    def test_batch_submission_is_one_workflow_per_batch(self):
        samples=[dict(dataset='d',sample_id=f'd/{n}') for n in range(3)]
        with patch.object(metag_submit,'load_config',return_value={}), \
             patch.object(metag_submit,'discover_samples',return_value=samples), \
             patch.object(metag_submit,'batch_plan',side_effect=lambda cfg,s,n,b:batches(s,n,b)), \
             patch.object(metag_submit.subprocess,'run') as run:
            metag_submit.main(['all','--batch','all','--batch-size','2','--plan-only','--submit','hpc'])
        self.assertEqual(run.call_count,2)
        for number,call in enumerate(run.call_args_list,1):
            cmd=call.args[0]
            self.assertIn('group,species,evaluate,report',cmd)
            self.assertEqual(cmd[cmd.index('--job-name')+1],f'metag.d.batch{number}')
            self.assertIn('--plan-only',cmd)
        self.assertEqual(run.call_args_list[0].args[0].count('--sample'),2)
        self.assertEqual(run.call_args_list[1].args[0].count('--sample'),1)

    def test_sample_submission_and_invalid_scope(self):
        with patch.object(metag_submit,'load_config',return_value={}), \
             patch.object(metag_submit,'discover_samples',return_value=[dict(dataset='d',sample_id='d/s')]), \
             patch.object(metag_submit.subprocess,'run') as run:
            metag_submit.main(['all','--sample','d/s','--plan-only'])
            self.assertEqual(run.call_args.args[0][run.call_args.args[0].index('--job-name')+1],'metag.sample')
            run.reset_mock()
            with self.assertRaises(SystemExit):metag_submit.main(['all','--sample','d/s','--batch','1'])
            run.assert_not_called()

    def test_full_workflow_waits_for_group_index(self):
        with tempfile.TemporaryDirectory() as temp:
            import json
            directory=Path(temp)/'jobs/qc/group_library/run'
            directory.mkdir(parents=True)
            (directory/'submission.json').write_text(json.dumps(dict(status='submitted',job_id='123',scope=dict(group_library_dir='/lib/group'))))
            scope=dict(module='hierarchical_split',action='group,species,evaluate,report',group_library_dir='/lib/group')
            with patch.object(qc_submit,'active_job',return_value=True):
                self.assertEqual(qc_submit.pending_hierarchy_jobs(dict(work_dir=temp),scope),['123'])


if __name__=='__main__':unittest.main()
