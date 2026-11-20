#!/usr/bin/env python3
"""Synthetic regressions for automatic stage evaluation; no large database needed."""
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import random
import tempfile
import unittest
import zipfile

import stage_evaluate as evaluation


def tsv(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(rows[0])
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'wt',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def pairs(left,right,records):
    Path(left).parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(left,'wt') as a,gzip.open(right,'wt') as b:
        for name,sequence,quality in records:
            for mate,out in ((1,a),(2,b)):
                out.write(f'@{name}/{mate}\n{sequence}\n+\n'+quality*len(sequence)+'\n')


def assignment(ordinal,group,species=0):
    return dict(ordinal=ordinal,read_id=f'r{ordinal}',group=group,taxid=species,species_taxid=species,
                status='assigned' if group!='unclassified' else 'unclassified',bin=group,source_group=group)


def fixture(root,with_references=True):
    sample=root/'dataset/sample';qc=evaluation.stage_qc(sample,'group');qc.mkdir(parents=True)
    rng=random.Random(1)
    seqs=[''.join(rng.choices('ACGT',weights=weights,k=150)) for weights in ((1,1,1,1),(3,1,1,3))]
    source=[('r1',seqs[0],'I'),('r2',seqs[1],'!')]
    samples=[];counts=[];assignments=[]
    for n,group in enumerate(('human','plant'),1):
        path=sample/'group'/group;path.mkdir(parents=True)
        left,right=path/'R1.fastq.gz',path/'R2.fastq.gz'
        pairs(left,right,[source[n-1]])
        row=assignment(n,group,9606 if group=='human' else 100)
        row['source_group']='';assignments.append(row)
        samples.append(row|dict(stratum=group,stratum_pairs=1,sampled_pairs=1))
        counts.append(row|dict(pairs=1,r1=str(left.relative_to(sample)),r2=str(right.relative_to(sample)),
                              bytes_r1=left.stat().st_size,bytes_r2=right.stat().st_size))
    pairs(qc/'sample.R1.fastq.gz',qc/'sample.R2.fastq.gz',source)
    tsv(qc/'sample.tsv',samples);tsv(qc/'counts.tsv',counts);tsv(qc/'assignments.tsv.gz',assignments)
    (qc/'summary.json').write_text(json.dumps(dict(sample_id='dataset/sample',input_pairs=2,written_pairs=2,
                                                assignment_rows=2,expected_pairs=2,paired=True)))
    cfg=dict(group_library_dir=str(root/'library'),work_dir=str(root/'analysis'),reference_raw=str(root/'refs'),
             evaluation_anchor_pairs=6,evaluation_anchor_references=1,evaluation_anchor_bases=2000)
    if with_references:
        rows=[]
        for tid,group,weights in ((9606,'human',(1,1,1,1)),(100,'plant',(3,1,1,3))):
            path=root/f'{group}.fa'
            path.write_text('>chr1\n'+''.join(rng.choices('ACGT',weights=weights,k=2500))+'\n')
            rows.append(dict(path=str(path),taxid=tid,source_taxid=tid,species_taxid=tid,group=group))
        tsv(root/'library/references.tsv',rows)
    return sample,qc,cfg


class FeatureTests(unittest.TestCase):
    def test_reverse_complement_and_mate_boundaries(self):
        self.assertEqual(len(set(evaluation._CANONICAL)),136)
        seq=b'AACGTGCTAGCAACGT'
        reverse=seq.translate(bytes.maketrans(b'ACGT',b'TGCA'))[::-1]
        self.assertEqual(evaluation.tetranucleotide_features([seq]),evaluation.tetranucleotide_features([reverse]))
        self.assertEqual(evaluation.tetranucleotide_features([b'AAA',b'CCC']),(None,0))
        self.assertEqual(evaluation.tetranucleotide_features([b'AAAANAAAA'])[1],2)

    def test_reference_prefix_sampling_preserves_contig_and_orientation(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'ref.fa'
            path.write_text('>short\nAAAA\n>chr1\nACGTACGTACGTNNNNACGTACGTACGT\n>later\n'+'G'*1000)
            rows,bases=evaluation.bounded_reference_windows(path,bases_limit=32,read_lengths=(4,4),
                insert_size=10,maximum=100,seed=7)
            self.assertEqual(bases,32)
            self.assertTrue(rows)
            self.assertTrue(all(row['contig']=='chr1' and len(row['sequence1'])==4 and len(row['sequence2'])==4 for row in rows))
            self.assertTrue(all('N' not in row['sequence1']+row['sequence2'] for row in rows))
            self.assertTrue(all(row['sampling_scope']=='bounded_reference_prefix_windows' for row in rows))

    def test_stratification_weights_are_used(self):
        with tempfile.TemporaryDirectory() as d:
            qc=Path(d);sequence='ACGTTGCAAGTCGATCAGGA'*5
            pairs(qc/'sample.R1.fastq.gz',qc/'sample.R2.fastq.gz',[('r2',sequence,'I'),('r10',sequence,'!')])
            rows=[assignment(2,'human')|dict(stratum='human',stratum_pairs=10,sampled_pairs=1),
                  assignment(10,'unclassified')|dict(stratum='unclassified',stratum_pairs=90,sampled_pairs=1)]
            tsv(qc/'sample.tsv',rows)
            report,observations,_=evaluation.sampled_data(qc)
            self.assertAlmostEqual(report['weighted_population_estimates']['q30_fraction'],.1)
            self.assertEqual(report['represented_pairs'],100)
            self.assertEqual([x['weight'] for x in observations],[10,90])

    def test_sample_identity_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            qc=Path(d);pairs(qc/'sample.R1.fastq.gz',qc/'sample.R2.fastq.gz',[('wrong','ACGT'*30,'I')])
            tsv(qc/'sample.tsv',[assignment(1,'human')|dict(stratum='human',stratum_pairs=1,sampled_pairs=1)])
            with self.assertRaisesRegex(RuntimeError,'ID mismatch'):evaluation.sampled_data(qc)


class TruthTests(unittest.TestCase):
    def test_strain_truth_is_compared_at_species_rank(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);taxonomy=root/'taxonomy';taxonomy.mkdir()
            records=[(1,1,'no rank'),(2,1,'superkingdom'),(100,2,'species'),(101,100,'strain')]
            (taxonomy/'nodes.dmp').write_text(''.join(f'{tid}\t|\t{parent}\t|\t{rank}\t|\n' for tid,parent,rank in records))
            (taxonomy/'names.dmp').write_text(''.join(f'{tid}\t|\tTaxon {tid}\t|\t\t|\tscientific name\t|\n' for tid,_,_ in records))
            tsv(root/'truth.tsv',[dict(ordinal=1,true_taxid=101,true_group='bacteria')])
            tsv(root/'assignments.tsv.gz',[assignment(1,'bacteria',100)])
            result=evaluation.evaluate_truth(root/'assignments.tsv.gz',root/'truth.tsv',root,'species',
                                             cfg=dict(taxonomy_dir=str(taxonomy)))
            self.assertEqual(result['recall'],1)
            self.assertEqual(result['taxid_normalization'],{'frozen_taxonomy_species_ancestor':1})

    def test_unsorted_assignments_missing_calls_and_unlabeled_reads(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            tsv(root/'truth.tsv',[dict(ordinal=1,true_group='human'),dict(ordinal=2,true_group='plant'),
                                 dict(ordinal=3,true_group='plant'),dict(ordinal=4,true_group='human')])
            tsv(root/'assignments.tsv.gz',[assignment(3,'unclassified'),assignment(1,'human'),
                                           assignment(2,'human'),assignment(5,'bacteria')])
            result=evaluation.evaluate_truth(root/'assignments.tsv.gz',root/'truth.tsv',root,'group')
            self.assertEqual(result['correct_pairs'],1)
            self.assertEqual(result['missing_predictions'],1)
            self.assertEqual(result['assignment_rows_without_truth'],1)
            self.assertAlmostEqual(result['precision_assigned'],.5)
            self.assertAlmostEqual(result['recall'],.25)
            self.assertFalse((root/'.truth.sqlite').exists())
            with evaluation.text_table(root/'truth-confusion.tsv') as rows:
                self.assertTrue(any(r['predicted_label']==evaluation.MISSING for r in rows))

    def test_duplicate_truth_or_matched_assignment_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);tsv(root/'truth.tsv',[dict(ordinal=1,true_group='human')]*2)
            tsv(root/'assignments.tsv.gz',[assignment(1,'human')])
            with self.assertRaisesRegex(ValueError,'Duplicate truth'):
                evaluation.evaluate_truth(root/'assignments.tsv.gz',root/'truth.tsv',root,'group')
            tsv(root/'truth.tsv',[dict(ordinal=1,true_group='human')])
            tsv(root/'assignments.tsv.gz',[assignment(1,'human')]*2)
            with self.assertRaisesRegex(ValueError,'Duplicate assignment'):
                evaluation.evaluate_truth(root/'assignments.tsv.gz',root/'truth.tsv',root,'group')

    def test_target_scope_counts_routing_misses_and_off_target_false_positives(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            tsv(root/'truth.tsv',[dict(ordinal=1,true_taxid=9606,true_group='human'),
                                 dict(ordinal=2,true_taxid=100,true_group='plant'),
                                 dict(ordinal=3,true_taxid=200,true_group='bacteria')])
            tsv(root/'assignments.tsv.gz',[assignment(1,'human',9606),assignment(3,'human',9606)])
            result=evaluation.evaluate_truth(root/'assignments.tsv.gz',root/'truth.tsv',root,'species',
                                             selected_groups=['human','plant'])
            target=result['selected_group_scope']
            self.assertEqual(target['true_pairs'],2)
            self.assertEqual(target['predicted_pairs'],2)
            self.assertEqual(target['precision'],.5)
            self.assertEqual(target['recall'],.5)


class WorkflowTests(unittest.TestCase):
    def test_existing_result_reuses_alignment_reference_audit_without_group_index(self):
        with tempfile.TemporaryDirectory() as d:
            sample,qc,cfg=fixture(Path(d))
            (qc/'alignment').mkdir()
            table=qc/'alignment/representatives.tsv'
            (Path(cfg['group_library_dir'])/'references.tsv').rename(table)
            (qc/'alignment/summary.json').write_text(json.dumps(dict(status='complete',reference_table=str(table))))
            result=evaluation.evaluate_stage(sample,'group',cfg)
            self.assertEqual(result['anchors']['status'],'complete')
            self.assertEqual(result['anchors']['source'],str(table))
            self.assertEqual(result['composition']['fit_scope'],'reference_anchors_only')
            self.assertTrue((qc/'controls.truth.tsv').is_file())

    def test_complete_evaluation_reference_fit_and_bounded_bundle(self):
        with tempfile.TemporaryDirectory() as d:
            sample,qc,cfg=fixture(Path(d))
            result=evaluation.evaluate_stage(sample,'group',cfg)
            self.assertEqual(result['status'],'complete_with_limitations')
            self.assertEqual(result['composition']['fit_scope'],'reference_anchors_only')
            self.assertEqual(result['controls']['status'],'generated_not_run')
            self.assertEqual(result['truth']['status'],'not_supplied')
            self.assertNotIn('accuracy',result['composition'])
            with zipfile.ZipFile(qc/'diagnostic.zip') as archive:
                self.assertIn('sample.R1.fastq.gz',archive.namelist())
                self.assertIn('composition.svg',archive.namelist())
                self.assertNotIn('assignments.tsv.gz',archive.namelist())
            self.assertTrue((qc/'evaluation.html').is_file())

    def test_missing_references_still_report_sample_only_pca(self):
        with tempfile.TemporaryDirectory() as d:
            sample,qc,cfg=fixture(Path(d),with_references=False)
            result=evaluation.evaluate_stage(sample,'group',cfg)
            self.assertEqual(result['status'],'complete_with_limitations')
            self.assertEqual(result['anchors']['status'],'not_available')
            self.assertEqual(result['composition']['fit_scope'],'sample_only')

    def test_file_count_failure_is_not_green(self):
        with tempfile.TemporaryDirectory() as d:
            sample,qc,cfg=fixture(Path(d));summary=json.loads((qc/'summary.json').read_text())
            summary['written_pairs']=3;(qc/'summary.json').write_text(json.dumps(summary))
            result=evaluation.evaluate_stage(sample,'group',cfg)
            self.assertEqual(result['status'],'incomplete')
            self.assertEqual(result['accounting']['status'],'failed')

    def test_species_accounting_uses_processed_inputs_including_uncertain(self):
        with tempfile.TemporaryDirectory() as d:
            sample,parent,cfg=fixture(Path(d),with_references=False)
            with evaluation.text_table(parent/'counts.tsv') as reader:parent_counts=list(reader)
            # Only human is a requested biological target, but the uncertain
            # branch is also refined; the plant parent pair remains preserved.
            parent_counts.append(dict(bin='group/unclassified',pairs=1,r1='',r2='',
                                      bytes_r1='',bytes_r2='',ordinal='',read_id='',
                                      group='unclassified',taxid=0,species_taxid=0,
                                      status='unclassified',source_group=''))
            tsv(parent/'counts.tsv',parent_counts)
            qc=evaluation.stage_qc(sample,'species');qc.mkdir()
            counts=[]
            for source,ordinal in (('human',1),('unclassified',3)):
                folder=sample/'group'/source/'species/9606'
                left,right=folder/'R1.fastq.gz',folder/'R2.fastq.gz'
                pairs(left,right,[(f'r{ordinal}','ACGT'*30,'I')])
                counts.append(dict(bin=str(folder.relative_to(sample)),pairs=1,
                    source_group=source,group='human',r1=str(left.relative_to(sample)),
                    r2=str(right.relative_to(sample)),bytes_r1=left.stat().st_size,bytes_r2=right.stat().st_size))
            tsv(qc/'counts.tsv',counts)
            tsv(qc/'assignments.tsv.gz',[assignment(1,'human',9606),
                assignment(3,'human',9606)|dict(source_group='unclassified')])
            summary=dict(input_pairs=2,written_pairs=2,assignment_rows=2,expected_pairs=3,
                         selected_groups=['human'],processed_groups=['human','unclassified'],unprocessed_pairs=1)
            (qc/'summary.json').write_text(json.dumps(summary))
            result,_,_=evaluation.accounting(sample,qc)
            self.assertEqual(result['status'],'complete')
            checked={row['check']:row for row in result['checks']}
            self.assertEqual(checked['species_input_vs_processed_parent_bins']['expected'],2)
            summary['processed_groups']=['human'];(qc/'summary.json').write_text(json.dumps(summary))
            result,_,_=evaluation.accounting(sample,qc)
            self.assertEqual(result['status'],'failed')

    def test_positive_control_recovery_requires_actual_predictions(self):
        with tempfile.TemporaryDirectory() as d:
            sample,qc,cfg=fixture(Path(d));evaluation.evaluate_stage(sample,'group',cfg)
            with evaluation.text_table(qc/'controls.truth.tsv') as reader:truth=list(reader)
            predictions=[assignment(int(row['ordinal']),row['true_group'],int(row['true_taxid']))|dict(read_id=row['read_id']) for row in truth]
            tsv(qc/'controls.assignments.tsv.gz',predictions)
            result=evaluation.evaluate_controls(qc,'group',cfg)
            self.assertEqual(result['precision_assigned'],1)
            self.assertEqual(result['recall'],1)
            self.assertTrue(result['source_in_classification_reference'])
            self.assertIn('not held-out',result['interpretation'])
            evaluation.prepare_reference_anchors(qc,cfg)
            self.assertEqual(evaluation.evaluate_controls(qc,'group',cfg)['status'],'generated_not_run')

    def test_species_controls_are_explicitly_conditional(self):
        with tempfile.TemporaryDirectory() as d:
            sample,qc,cfg=fixture(Path(d));evaluation.prepare_reference_anchors(qc,cfg)
            summary=json.loads((qc/'summary.json').read_text());summary['selected_groups']=['human']
            (qc/'summary.json').write_text(json.dumps(summary))
            with evaluation.text_table(qc/'controls.truth.tsv') as reader:truth=list(reader)
            predictions=[assignment(int(row['ordinal']),'human',9606)|dict(read_id=row['read_id'])
                         for row in truth if row['true_group']=='human']
            tsv(qc/'controls.assignments.tsv.gz',predictions)
            result=evaluation.evaluate_controls(qc,'species',cfg)
            self.assertEqual(result['control_scope'],'conditional_on_known_source_group')
            self.assertEqual(result['recall'],1)

    def test_alignment_and_agreement_are_included_without_accuracy_label(self):
        with tempfile.TemporaryDirectory() as d:
            sample,qc,cfg=fixture(Path(d))
            (qc/'alignment').mkdir();(qc/'agreement').mkdir()
            (qc/'alignment/summary.json').write_text(json.dumps(dict(status='complete',accuracy_estimated=False)))
            (qc/'alignment/panels.json').write_text(json.dumps(dict(panels=[dict(panel_id='p1')],
                interpretation='Competition is limited to the selected reference panels.')))
            tsv(qc/'alignment/omitted-references.tsv',[dict(species_taxid=9598,reason='reference_unavailable')])
            (qc/'agreement/summary.json').write_text(json.dumps(dict(status='complete',agreement_fraction=.95)))
            tsv(qc/'agreement/pairs.tsv',[dict(ordinal=1,agreement=True)])
            result=evaluation.evaluate_stage(sample,'group',cfg)
            self.assertFalse(result['alignment']['accuracy_estimated'])
            self.assertEqual(result['agreement']['agreement_fraction'],.95)
            with zipfile.ZipFile(qc/'diagnostic.zip') as archive:
                self.assertIn('agreement/pairs.tsv',archive.namelist())
                self.assertIn('alignment/summary.json',archive.namelist())
                self.assertIn('alignment/panels.json',archive.namelist())
                self.assertIn('alignment/omitted-references.tsv',archive.namelist())
                self.assertIn(b'reference_unavailable',archive.read('alignment/omitted-references.tsv'))

    def test_corrupt_optional_alignment_summary_does_not_discard_other_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            sample,qc,cfg=fixture(Path(d))
            (qc/'alignment').mkdir();(qc/'alignment/summary.json').write_text('{broken')
            result=evaluation.evaluate_stage(sample,'group',cfg)
            self.assertEqual(result['status'],'incomplete')
            self.assertEqual(result['alignment']['status'],'failed')
            self.assertEqual(result['accounting']['status'],'complete')
            self.assertTrue((qc/'diagnostic.zip').is_file())


if __name__=='__main__':unittest.main(verbosity=2)
