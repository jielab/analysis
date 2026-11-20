"""Regression coverage for active hierarchy dependencies; no scheduler access."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import qc_submit


class PendingHierarchyDependencyTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cfg = {"work_dir": str(self.root)}
        self.scope = dict(module="hierarchical_split", action="species",
                          samples=["cohort_a/sample_1"], datasets=["cohort_a"],
                          output_dir="/analysis/current", library_dir="/library/main",
                          group_library_dir="/library/group")
        self.sequence = 0

    def receipt(self, module, job, *, action="group", status="submitted", **scope):
        self.sequence += 1
        path = self.root / "jobs/qc" / module / f"run{self.sequence:04d}" / "submission.json"
        path.parent.mkdir(parents=True)
        prior = dict(self.scope, module=module, action=action)
        prior.update(scope)
        path.write_text(json.dumps(dict(scope=prior, status=status, job_id=job)))
        # Deterministic newest-first ordering, independent of filesystem timing.
        stamp = 1_700_000_000_000_000_000 + self.sequence * 1_000_000
        os.utime(path, ns=(stamp, stamp))

    def pending(self, states, **scope):
        queried = []

        def bjobs(command, **kwargs):
            self.assertEqual(command[:4], ["bjobs", "-noheader", "-o", "stat"])
            job = command[4]
            self.assertIn(job, states, "Irrelevant submission queried: " + job)
            queried.append(job)
            return subprocess.CompletedProcess(command, 0, states[job] + "\n", "")

        with patch.object(qc_submit.subprocess, "run", side_effect=bjobs):
            jobs = qc_submit.pending_hierarchy_jobs(self.cfg, dict(self.scope, **scope))
        return jobs, queried

    def test_group_library_keeps_all_active_ids_and_deduplicates_receipts(self):
        self.receipt("group_library", "101")
        self.receipt("group_library", "102")
        self.receipt("group_library", "101")
        self.receipt("group_library", "901", group_library_dir="/library/other")
        self.receipt("group_library", "902", status="planned")
        jobs, queried = self.pending({"101": "RUN", "102": "PEND"}, action="group")
        self.assertEqual(jobs, ["101", "102"])
        self.assertEqual(set(queried), {"101", "102"})

    def test_completed_latest_group_library_job_does_not_hide_older_active_job(self):
        self.receipt("group_library", "201")
        self.receipt("group_library", "202")
        jobs, queried = self.pending({"201": "RUN", "202": "DONE"}, action="group")
        self.assertEqual(jobs, ["201"])
        self.assertEqual(queried, ["202", "201"])

    def test_hierarchy_keeps_all_active_jobs_for_identical_scope(self):
        self.receipt("hierarchical_split", "301")
        self.receipt("hierarchical_split", "302")
        self.receipt("hierarchical_split", "301")
        jobs, queried = self.pending({"301": "RUN", "302": "PEND"})
        self.assertEqual(jobs, ["301", "302"])
        self.assertEqual(set(queried), {"301", "302"})

    def test_completed_latest_hierarchy_job_does_not_hide_older_active_job(self):
        self.receipt("hierarchical_split", "401")
        self.receipt("hierarchical_split", "402")
        jobs, queried = self.pending({"401": "RUN", "402": "DONE"})
        self.assertEqual(jobs, ["401"])
        self.assertEqual(queried, ["402", "401"])

    def test_complete_receipt_does_not_shadow_older_submitted_job(self):
        self.receipt("hierarchical_split", "501")
        self.receipt("hierarchical_split", "502", status="complete")
        jobs, queried = self.pending({"501": "RUN"})
        self.assertEqual(jobs, ["501"])
        self.assertEqual(queried, ["501"])

    def test_hierarchy_preserves_output_dataset_sample_and_stage_scope(self):
        self.receipt("hierarchical_split", "601", action="group")
        self.receipt("hierarchical_split", "602", action="species", samples=[])
        self.receipt("hierarchical_split", "601", action="group,species")
        self.receipt("hierarchical_split", "901", output_dir="/analysis/other")
        self.receipt("hierarchical_split", "902", datasets=["cohort_b"], samples=[])
        self.receipt("hierarchical_split", "903", samples=["cohort_a/sample_2"])
        self.receipt("hierarchical_split", "904", action="report")
        self.receipt("hierarchical_split", "905", status="planned")
        jobs, queried = self.pending({"601": "RUN", "602": "PEND"}, action="evaluate")
        self.assertEqual(jobs, ["601", "602"])
        self.assertEqual(set(queried), {"601", "602"})

    def test_report_and_unrelated_module_do_not_wait_for_hierarchy_jobs(self):
        self.receipt("hierarchical_split", "701")
        self.receipt("group_library", "702")
        for scope in (dict(action="report"), dict(module="group_library", action="group")):
            with self.subTest(scope=scope):
                jobs, queried = self.pending({}, **scope)
                self.assertEqual(jobs, [])
                self.assertEqual(queried, [])


if __name__ == "__main__":
    unittest.main()
