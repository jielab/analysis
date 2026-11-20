#!/usr/bin/env python3
"""Focused regression tests for low-memory admission and suite completion."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from common import effective_workers
import validation_suite as suite
from parallel_cli import parser


class SuiteTests(unittest.TestCase):
    def test_large_mapped_index_has_no_full_ram_floor(self):
        # Sparse files exercise admission using realistic sizes, not real RAM.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, size in [("a.k2d", 1200 * 2**30), ("b.k2d", 350 * 2**30),
                               ("taxo.k2d", 2**20), ("opts.k2d", 64)]:
                with (root / name).open("wb") as handle:
                    handle.truncate(size)
            snapshot = dict(layout="partitioned", hashes=["a.k2d", "b.k2d"], db_dir=str(root))
            cfg = dict(classifier_engine="kraken", union_memory_mapping=True,
                       classify_memory_gib=128, system_reserve_gib=16, worker_memory_gib=8,
                       classify_cores=24, classify_threads=22, classify_workers=3,
                       memory_accounting="conservative")
            with contextlib.redirect_stdout(io.StringIO()):
                admitted = effective_workers(cfg, snapshot)
            self.assertEqual(admitted["classify_workers"], 1)
            self.assertFalse(admitted["classification_memory_policy"]["rss_limit_enforced_by_mmap"])
            with self.assertRaisesRegex(RuntimeError, "Complete index"):
                effective_workers(dict(cfg, union_memory_mapping=False), snapshot)
            with self.assertRaisesRegex(RuntimeError, "still needs"):
                effective_workers(dict(cfg, classify_memory_gib=20), snapshot)

    def test_new_split_flags_parse(self):
        args = parser().parse_args(["split", "--sample", "england/ERR14374050", "--memory-mode", "mapped", "--mem-gib", "128"])
        self.assertEqual((args.memory_mode, args.mem_gib), ("mapped", 128))

    def test_failures_are_retained_and_later_steps_execute(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            def runner(command, stdout, stderr):
                stdout.write("test step\n")
                return subprocess.CompletedProcess(command, 3 if command[0] == "fail" else 0)
            first = suite.execute_step(dict(name="first", command=["fail"], out=root / "first"), logs=root / "logs", runner=runner)
            second = suite.execute_step(dict(name="second", command=["ok"], out=root / "second"), logs=root / "logs", runner=runner)
            self.assertEqual(first["status"], "failed")
            self.assertEqual(second["status"], "complete")
            self.assertIn("test step", first["note"])

    def test_partial_coverage_cannot_be_reported_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "plant-validation.json").write_text(json.dumps(dict(status="not_tested", untested_fragments=77)))
            result = dict(output=str(root), status="complete")
            suite.annotate_coverage(result)
            self.assertEqual(result["status"], "partially_tested")
            self.assertIn("77", result["note"])

    def test_reuse_requires_all_declared_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "x.tsv").write_text("abc\n")
            (root / "complete.json").write_text(json.dumps(dict(status="complete", files={"x.tsv":4})))
            self.assertTrue(suite.intact_result(root))
            (root / "x.tsv").write_text("bad")
            self.assertFalse(suite.intact_result(root))
            (root / "complete.json").write_text(json.dumps(dict(status="complete", files={"../outside":0})))
            self.assertFalse(suite.intact_result(root))

    def test_blocked_panel_does_not_run_with_stale_references(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with patch.object(suite.subprocess, "run", side_effect=AssertionError("must not run")):
                result = suite.execute_step(dict(name="plant", command=["unsafe"], out=root / "plant",
                                                blocked="metadata failed"), logs=root / "logs")
            self.assertEqual(result["status"], "blocked")

    def test_reuse_detects_changed_reference_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "ref.fa"
            reference.write_text(">r\nACGT\n")
            stats = reference.stat()
            (root / "evidence.tsv").write_text("ok\n")
            (root / "complete.json").write_text(json.dumps(dict(status="complete", files={"evidence.tsv":3},
                input_states=[dict(path=str(reference), size=stats.st_size, mtime_ns=stats.st_mtime_ns, ctime_ns=stats.st_ctime_ns)])))
            self.assertTrue(suite.intact_result(root))
            reference.write_text(">r\nTGCA\n")
            # Filesystem timestamp ticks can coalesce immediate equal-size writes.
            import os
            os.utime(reference, ns=(stats.st_atime_ns, stats.st_mtime_ns + 1_000_000_000))
            self.assertFalse(suite.intact_result(root))

    def test_identical_reference_audit_preserves_table_state(self):
        from reference_audit import _write_tsv
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "representatives.tsv"
            _write_tsv(path, [{"taxid":9606}], ["taxid"])
            original = path.stat().st_mtime_ns
            _write_tsv(path, [{"taxid":9606}], ["taxid"])
            self.assertEqual(original, path.stat().st_mtime_ns)
            _write_tsv(path, [{"taxid":9598}], ["taxid"])
            self.assertIn("9598", path.read_text())

    def test_active_human_index_is_a_success_dependency(self):
        from qc_submit import pending_human_index
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "index.jobid").write_text("12345\n")
            cfg = {"human_reference_dir":str(root)}
            with patch("qc_submit.subprocess.run", return_value=subprocess.CompletedProcess([], 0, "RUN\n", "")):
                self.assertEqual(pending_human_index(cfg), "12345")
            with patch("qc_submit.subprocess.run", return_value=subprocess.CompletedProcess([], 0, "EXIT\n", "")):
                self.assertEqual(pending_human_index(cfg), "")
            fasta = root / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
            for file in [fasta, root / "BUILD.COMPLETE", *[Path(str(fasta)+s) for s in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")]]:
                file.write_text("ready")
            with patch("qc_submit.subprocess.run", side_effect=AssertionError("completed index needs no scheduler query")):
                self.assertEqual(pending_human_index(cfg), "")

    def test_report_escapes_identifiers_and_error_messages(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report = dict(sample_id="dataset/<sample>", status="incomplete", split_summary={"fragments":4},
                          status_fragments={"unclassified":4}, steps=[dict(step="plant", status="failed",
                          output=str(root / "plant"), log=str(root / "missing"), note="<script>bad</script>")])
            suite.render_report(root, report)
            page = (root / "report.html").read_text()
            self.assertIn("&lt;sample&gt;", page)
            self.assertNotIn("<script>", page)
            self.assertIn("100.000000%", page)


if __name__ == "__main__":
    unittest.main(verbosity=2)
