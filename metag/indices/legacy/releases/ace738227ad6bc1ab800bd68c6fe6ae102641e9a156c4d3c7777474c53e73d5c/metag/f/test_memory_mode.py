#!/usr/bin/env python3
"""Focused regression tests for mapped execution and its memory measurements.

Run: python3 f/test_memory_mode.py
The streaming integration uses a deterministic native-process stand-in to test
ordinal ordering and telemetry. Actual union lookup equivalence is separately
tested by test_union_memory.cc with the locally installed Kraken source.
"""
import gzip
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ordinal_query import (_read_process_metrics, NativeTelemetry,
                           classifier_command, query_many)


class MemoryModeTests(unittest.TestCase):
    def test_only_loading_flag_changes(self):
        snapshot = {"db_dir": "/example", "layout": "partitioned", "hashes": ["a/hash.k2d"]}
        cfg = {"database_layout": "partitioned", "classify_threads": 3,
               "confidence": 0.1, "minimum_hit_groups": 2, "minimum_base_quality": 7}
        ram = classifier_command("classifier", snapshot, dict(cfg, union_memory_mapping=False), True)
        mapped = classifier_command("classifier", snapshot, dict(cfg, union_memory_mapping=True), True)
        self.assertNotIn("-M", ram)
        self.assertEqual(ram, [x for x in mapped if x != "-M"])
        self.assertIn("--metag-hashes", mapped)
        self.assertIn("-S", mapped)
        with self.assertRaises(ValueError):
            classifier_command("classifier", snapshot, dict(cfg, union_memory_mapping="false"), True)

    def test_procfs_parsing_and_missing_permissions(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc = Path(tmp) / "123"
            proc.mkdir()
            (proc / "status").write_text("VmRSS:\t50 kB\nVmHWM:\t100 kB\nVmPeak:\t200 kB\n"
                                         "VmSize:\t150 kB\nRssAnon:\t40 kB\nRssFile:\t10 kB\n")
            # PID and comm, then fields 3..12. comm deliberately contains ')'.
            (proc / "stat").write_text("123 (name ) with spaces) S 1 2 3 4 5 6 77 8 99\n")
            (proc / "io").write_text("rchar: 5000\nread_bytes: 4096\n")
            data = _read_process_metrics(123, tmp)
            self.assertEqual(data["rss_peak_observed_bytes"], 100 * 1024)
            self.assertEqual(data["virtual_peak_observed_bytes"], 200 * 1024)
            self.assertEqual(data["minor_faults_observed"], 77)
            self.assertEqual(data["major_faults_observed"], 99)
            self.assertEqual(data["read_bytes_observed"], 4096)
            (proc / "io").unlink()
            self.assertNotIn("read_bytes_observed", _read_process_metrics(123, tmp))
            self.assertEqual(_read_process_metrics(456, tmp), {})

    @unittest.skipUnless(Path("/proc/self/status").is_file(), "Linux procfs is optional")
    def test_live_native_memory(self):
        child = subprocess.Popen([sys.executable, "-c",
            "import sys; data=bytearray(24*1024*1024); print('ready', flush=True); sys.stdin.readline()"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
        try:
            self.assertEqual(child.stdout.readline().strip(), "ready")
            monitor = NativeTelemetry(child.pid, interval=0.02)
            time.sleep(0.06)
            metrics = monitor.finish()
            self.assertTrue(metrics["available"])
            self.assertGreater(metrics["rss_peak_observed_bytes"], 24 * 1024 * 1024)
            self.assertGreaterEqual(metrics["virtual_peak_observed_bytes"], metrics["rss_peak_observed_bytes"])
            self.assertGreater(metrics["samples"], 0)
        finally:
            child.communicate("stop\n", timeout=5)

    def test_mapped_batch_stream_preserves_exact_pair_ordinals(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            native = root / "native.py"
            native.write_text("#!" + sys.executable + "\n" +
                "import sys\n"
                "rows=sys.stdin.buffer.readlines()\n"
                "for i in range(len(rows)-8, -1, -8):\n"
                " name=rows[i][1:].strip().decode(); a=rows[i+1].strip(); b=rows[i+5].strip()\n"
                " tax=9606 if a[:1]==b'A' else 2 if a[:1]==b'C' else 0\n"
                " status='C' if tax else 'U'\n"
                " print(f'{status}\\t{name}\\t{tax}\\t{len(a)}|{len(b)}\\t0:1')\n")
            native.chmod(0o755)
            samples = []
            for name, bases in (("one", ["A", "N"]), ("two", ["C"])):
                sample = {"sample_id": "fixture/" + name, "sample": name, "dataset": "fixture"}
                for mate in (1, 2):
                    path = root / (name + ".R" + str(mate) + ".fastq.gz")
                    with gzip.open(path, "wt") as handle:
                        for base in bases:
                            # Reused IDs across samples/records are safe because
                            # classifier/split communication is ordinal based.
                            handle.write(f"@duplicate/{mate}\n{base * 8}\n+\nIIIIIIII\n")
                    sample["r" + str(mate)] = str(path)
                samples.append(sample)
            snapshot = {"db_dir": str(root), "layout": "partitioned", "hashes": ["hash.k2d"]}
            cfg = {"database_layout": "partitioned", "classify_threads": 1, "min_free_gib": 0}
            tax = SimpleNamespace(describe=lambda tid: {"taxid": tid})
            result = []
            with patch("ordinal_query.validate_snapshot"), patch("ordinal_query.resolve_native", return_value=native):
                for mapped in (False, True):
                    work = root / ("mapped" if mapped else "ram")
                    work.mkdir()
                    queries = query_many(dict(cfg, union_memory_mapping=mapped), samples, snapshot, work, taxonomy=tax)
                    result.append([Path(q["taxids"]).read_bytes() for q in queries])
                    self.assertEqual([q["fragments"] for q in queries], [2, 1])
                    self.assertTrue(all(q["memory_mode"] == ("mmap" if mapped else "ram") for q in queries))
                    self.assertTrue(all("native_metrics" in q for q in queries))
                    log = (work / "query.log").read_text()
                    self.assertIn("NATIVE_METRICS ", log)
                    self.assertEqual("MEMORY_NOTE " in log, mapped)
            self.assertEqual(result[0], result[1])
            self.assertEqual(struct.unpack("<QQ", result[0][0]), (9606, 0))
            self.assertEqual(struct.unpack("<Q", result[0][1]), (2,))


if __name__ == "__main__":
    unittest.main(verbosity=2)
