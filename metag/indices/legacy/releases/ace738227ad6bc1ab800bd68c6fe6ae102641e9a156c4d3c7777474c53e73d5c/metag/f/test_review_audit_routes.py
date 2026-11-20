#!/usr/bin/env python3
"""Reference-audit shell routes on tiny frozen metadata; never build a database."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


CODE = Path(__file__).resolve().parents[1]


class ReferenceAuditRouteTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="metag audit routes ")
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.root = self.base / "data root"
        self.config = self.base / "audit config.json"
        self.config.write_text(json.dumps({
            "schema": 1,
            # The explicit --root must override this deliberately nonexistent root.
            "root": str(self.base / "wrong root"),
            "storage_layout": "separated",
            "database_layout": "partitioned",
            "classifier_engine": "kraken",
            "build_kunpeng": False,
            "library_dir": "library",
            "metadata_dir": "metadata",
            "reference_raw": "raw",
            "reads_raw": "reads",
            "output_dir": "output",
            "work_dir": "state",
            "tools_dir": "tools",
        }))
        self.env = dict(os.environ)
        for key in ("METAG_ENV_FILE", "METAG_HPC_CONFIG", "METAG_CONFIG", "AUDIT_DIR",
                    "SYLPH_DIR", "SYLPH_TOOL_DIR", "parallel_reference_dir"):
            self.env.pop(key, None)
        self.env.update(METAG_PYTHON=sys.executable, submit="local", SUBMIT="local",
                        SUBMIT_WAIT="FALSE", jobs_N="1", PYTHONDONTWRITEBYTECODE="1")

    def run_shell(self, script, *args):
        return subprocess.run(["bash", str(CODE / script), *map(str, args)],
                              cwd=self.base, env=self.env, text=True,
                              capture_output=True, timeout=30)

    def fixture(self, human=True):
        library = self.root / "library"
        taxonomy = self.root / "metadata/taxonomy"
        taxonomy.mkdir(parents=True)
        nodes = [(1, 1, "no rank"), (2, 1, "superkingdom"), (100, 2, "species")]
        if human:
            nodes.append((9606, 1, "species"))
        (taxonomy / "nodes.dmp").write_text("".join(
            f"{tid}\t|\t{parent}\t|\t{rank}\t|\n" for tid, parent, rank in nodes))
        (taxonomy / "names.dmp").write_text("".join(
            f"{tid}\t|\tTaxon {tid}\t|\t\t|\tscientific name\t|\n" for tid, _, _ in nodes))
        rows = [{"path": "bacteria/reference.fna.gz", "taxid": 100}]
        if human:
            rows.append({"path": "human/reference.fna.gz", "taxid": 9606})
        with (self.root / "metadata/reference_taxids.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["path", "taxid"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        receipt_dir = library / "prepared/task000001"
        receipt_dir.mkdir(parents=True)
        (library / "plan.json").write_text(json.dumps({
            "tasks": [{"task_id": 1, "manifest": {"path": "task.json"}}]}))
        (library / "task.json").write_text(json.dumps({"files": rows}))
        (receipt_dir / "receipt.json").write_text(json.dumps({"context": {"stats": {
            "sources": [{"path": row["path"], "records": 1, "bases": 100}
                        for row in rows]}}}))
        # No reference FASTA exists: successful audit must use receipts only.
        return len(rows)

    def test_reference_help_uses_the_reference_worker(self):
        for flag in ("--help", "-h", "--help-full"):
            with self.subTest(flag=flag):
                result = self.run_shell("f/analysis_entry.sh", "check", "--references", flag)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("--require-human", result.stdout)
                self.assertIn("--skip-assembly-metadata", result.stdout)

    def test_public_reference_check_forwards_config_root_out_and_options(self):
        expected = self.fixture()
        out = self.base / "audit result"
        result = self.run_shell("1.lib.sh", "check", "--references",
                                "--config", self.config, "--root", self.root,
                                "--out", out, "--skip-assembly-metadata", "--require-human",
                                "--submit", "local")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads((out / "reference-audit.json").read_text())
        self.assertEqual(report["references"], expected)
        self.assertEqual(report["bases"], expected * 100)
        self.assertTrue(report["human"]["present"])
        self.assertNotIn("assembly_metadata_matches", report)
        self.assertTrue((out / "species-representatives.tsv").is_file())
        self.assertFalse((self.root / "raw").exists())

    def test_default_audit_out_is_derived_from_the_overridden_root(self):
        self.fixture()
        result = self.run_shell("f/analysis_entry.sh", "check", "--references",
                                "--config", self.config, "--root", self.root)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((self.root / "state/qc/reference-audit/reference-audit.json").is_file())
        self.assertFalse((self.base / "wrong root").exists())

    def test_reference_worker_failure_is_propagated(self):
        self.fixture(human=False)
        out = self.base / "missing human report"
        result = self.run_shell("f/analysis_entry.sh", "check", "--references",
                                "--config", self.config, "--root", self.root,
                                "--out", out, "--require-human")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no frozen human reference lineage was found", result.stderr)
        self.assertFalse(json.loads((out / "reference-audit.json").read_text())["human"]["present"])

    def test_sylph_runs_the_reference_audit_before_the_database_planner(self):
        expected = self.fixture()
        sylph = self.root / "tools/sylph/1.0.0/sylph"
        sylph.parent.mkdir(parents=True)
        sylph.write_text("#!/usr/bin/env bash\n"
                         "case ${1:-} in\n"
                         "  --version) echo 'sylph 1.0.0' ;;\n"
                         "  -h) exit 0 ;;\n"
                         "  *) echo 'Unexpected Sylph build invocation' >&2; exit 99 ;;\n"
                         "esac\n")
        sylph.chmod(0o755)
        calls = self.base / "python calls.jsonl"
        planner_args = self.base / "planner args.json"
        selected_python = self.base / "selected-python"
        selected_python.write_text(
            f"#!{sys.executable}\n"
            "import json, os, pathlib, sys\n"
            "with open(os.environ['AUDIT_TEST_CALLS'], 'a') as handle:\n"
            "    handle.write(json.dumps(sys.argv[1:]) + '\\n')\n"
            "if pathlib.Path(sys.argv[1]).name == 'sylph_profile.py':\n"
            "    pathlib.Path(os.environ['AUDIT_TEST_PLANNER']).write_text(json.dumps(sys.argv[1:]))\n"
            "    raise SystemExit(73)  # Stop before any database planning/build/submission.\n"
            "real = os.environ['AUDIT_TEST_REAL_PYTHON']\n"
            "os.execv(real, [real, *sys.argv[1:]])\n")
        selected_python.chmod(0o755)
        self.env.update(METAG_PYTHON=str(selected_python), AUDIT_TEST_REAL_PYTHON=sys.executable,
                        AUDIT_TEST_CALLS=str(calls), AUDIT_TEST_PLANNER=str(planner_args))
        result = self.run_shell("1.lib.sh", "sylph", "--config", self.config,
                                "--root", self.root, "--submit", "local")
        self.assertEqual(result.returncode, 73, result.stdout + result.stderr)
        audit_dir = self.root / "state/qc/reference-audit"
        self.assertEqual(json.loads((audit_dir / "reference-audit.json").read_text())["references"], expected)
        argv = json.loads(planner_args.read_text())
        self.assertEqual(argv[1], "plan-db")
        self.assertEqual(argv[argv.index("--config") + 1], str(self.config))
        self.assertEqual(argv[argv.index("--root") + 1], str(self.root))
        self.assertEqual(argv[argv.index("--representatives") + 1], str(audit_dir / "species-representatives.tsv"))
        dispatched = [json.loads(line) for line in calls.read_text().splitlines()]
        self.assertTrue(any(Path(args[0]).name == "reference_audit.py" for args in dispatched))
        self.assertFalse((self.root / "raw").exists())
        self.assertFalse(list(self.root.rglob("*.BUILD.COMPLETE")))


if __name__ == "__main__":
    unittest.main()
