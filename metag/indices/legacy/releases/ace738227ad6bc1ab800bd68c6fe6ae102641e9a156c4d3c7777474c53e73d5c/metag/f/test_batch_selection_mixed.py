"""Unequal dataset sizes must not block valid batches in other datasets."""
import tempfile
import unittest
from pathlib import Path
from workflow_selection import batches
from project_cache import batch_plan


class MixedBatchTests(unittest.TestCase):
    def test_mixed_sizes_cached_and_uncached(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reads = root / 'reads'; reads.mkdir()
            r1 = reads / 'r1'; r1.write_text('fixture')
            r2 = reads / 'r2'; r2.write_text('fixture')
            cfg = dict(reads_raw=str(reads), work_dir=str(root / 'work'))
            samples = [dict(dataset=d, sample_id=f'{d}/{n:03}', r1=str(r1), r2=str(r2))
                       for d,count in [('small',10),('large',25)] for n in range(count)]
            for selection,expected in [
                ('1-2',[('large',1,10),('large',2,10),('small',1,10)]),
                ('2',[('large',2,10)]),
                ('2-',[('large',2,10),('large',3,5)]),
                ('all',[('large',1,10),('large',2,10),('large',3,5),('small',1,10)]),
            ]:
                for make in (lambda: batches(samples,10,selection),
                             lambda: batch_plan(cfg,samples,10,selection)):
                    with self.subTest(selection=selection):
                        self.assertEqual([(d,n,len(s)) for d,n,s in make()],expected)
                        self.assertEqual([(d,n,len(s)) for d,n,s in make()],expected)
            for selection in ('2-1','1,1','1-2,2','0'):
                with self.subTest(invalid=selection):
                    with self.assertRaises(ValueError): batches(samples,10,selection)
                    with self.assertRaises(ValueError): batch_plan(cfg,samples,10,selection)


class AvailableOnlyTests(unittest.TestCase):
    def test_requested_two_batches_uses_whatever_is_available(self):
        for count in (1, 10, 15, 20, 25):
            with self.subTest(count=count), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                reads = root / 'reads'; reads.mkdir()
                r1 = reads / 'r1'; r1.write_text('fixture')
                r2 = reads / 'r2'; r2.write_text('fixture')
                cfg = dict(reads_raw=str(reads), work_dir=str(root / 'work'))
                samples = [dict(dataset='d', sample_id=f'd/{n:03}', r1=str(r1), r2=str(r2))
                           for n in range(count)]
                for plan in (batches(samples,10,'1-2'), batch_plan(cfg,samples,10,'1-2')):
                    self.assertEqual(sum(len(members) for _,_,members in plan),min(count,20))
                    self.assertTrue(all(len(members) <= 10 for _,_,members in plan))
                self.assertEqual(batches(samples,10,'100-'),[])
                self.assertEqual(batch_plan(cfg,samples,10,'100-'),[])
