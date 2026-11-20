#!/usr/bin/env python3
"""Group-index safety and real pinned-native integration tests.

python3 f/test_group_library.py
METAG_TEST_KRAKEN_BIN=/path/to/pinned/kraken/bin python3 f/test_group_library.py

The optional integration actually estimates and builds capped and uncapped
indexes, then classifies paired human, plant, bacterial, shared and absent DNA.
It never downloads, masks, or modifies an existing library.
"""
import csv
import json
import os
from pathlib import Path
import random
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from common import atomic_json, fingerprint, sha256
import group_library as group
from ordinal_query import query_many
from reference import GIB, SOURCE_COMMIT


def fixture(root, native=None):
    root = Path(root)
    lib, meta, raw, work, tools = [root / name for name in ("lib", "metadata", "raw", "work", "tools")]
    tax = lib / "global" / "taxonomy"
    tax.mkdir(parents=True)
    for folder in (meta, raw, work, tools, lib / "manifests", lib / "prepared"):
        folder.mkdir(parents=True, exist_ok=True)
    nodes = [(1, 1, "no rank", "root"), (2, 1, "superkingdom", "Bacteria"),
             (2759, 1, "superkingdom", "Eukaryota"), (33090, 2759, "kingdom", "Viridiplantae"),
             (33208, 2759, "kingdom", "Metazoa"), (7742, 33208, "subphylum", "Vertebrata"),
             (40674, 7742, "class", "Mammalia"), (9606, 40674, "species", "Homo sapiens"),
             (3702, 33090, "species", "Arabidopsis thaliana"), (562, 2, "species", "Escherichia coli")]
    (tax / "nodes.dmp").write_text("".join(f"{tid}\t|\t{parent}\t|\t{rank}\t|\n" for tid, parent, rank, _ in nodes))
    (tax / "names.dmp").write_text("".join(f"{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n" for tid, _, _, name in nodes))
    (tax / "merged.dmp").write_text("")
    shutil.copytree(tax, meta / "taxonomy")
    generator = random.Random(418)
    random_dna = lambda n: "".join(generator.choices("ACGT", k=n))
    sequences = {9606: random_dna(24000), 3702: random_dna(24000), 562: random_dna(24000)}
    shared = random_dna(2000)
    groups, tasks, mapping, receipt_ids = {}, [], [], []
    for task_id, (category, tid) in enumerate((("vertebrates", 9606), ("plants", 3702), ("microbes", 562)), 1):
        parts = [sequences[tid]] + ([shared] if tid in (9606, 3702) else [])
        records = []
        for seq_id, sequence in enumerate(parts, 1):
            name = f"metag_f{task_id:010d}_r{seq_id:010d}|kraken:taxid|{tid}"
            records.append(f">{name}\n{sequence}")
            mapping.append(f"{name}\t{tid}\n")
        contents = "\n".join(records)  # Deliberately no final newline.
        task = lib / "prepared" / f"task{task_id:06d}"
        task.mkdir()
        fasta = task / "sequences.fna"
        fasta.write_text(contents)
        reference = raw / f"reference{task_id}.fna"
        reference.write_text(contents)
        bases = sum(map(len, parts))
        row = {"path": reference.name, "taxid": tid, "accession": f"GCF_00000000{task_id}.1",
               "category": category, "state": {"size": reference.stat().st_size}}
        manifest_path = f"manifests/task{task_id:06d}.json"
        atomic_json(lib / manifest_path, {"task_id": task_id, "files": [row]})
        receipt_context = {"stats": {"sources": [
            {"path": reference.name, "records": len(parts), "bases": bases}]}}
        receipt_id = fingerprint(receipt_context)
        atomic_json(task / "receipt.json", {"status": "complete", "context": receipt_context, "identity": receipt_id})
        receipt_ids.append(receipt_id)
        tasks.append({"task_id": task_id, "manifest": {"path": manifest_path,
                      "size": (lib / manifest_path).stat().st_size, "sha256": sha256(lib / manifest_path)}})
        item = {"path": str(fasta), "size": fasta.stat().st_size, "sha256": sha256(fasta),
                "records": len(parts), "bases": bases}
        groups[category] = {"inputs": [item], "inputs_tsv": str(lib / (category + ".inputs.tsv")),
                            "records": len(parts), "bases": bases, "longest_record": max(map(len, parts))}
    (lib / "global" / "seqid2taxid.map").write_text("".join(mapping))
    plan_context = {"tasks": tasks}
    plan_id = fingerprint(plan_context)
    atomic_json(lib / "plan.json", {"status": "complete", "tasks": tasks, "context": plan_context, "identity": plan_id})
    files = {str(path.relative_to(lib / "global")): {"size": path.stat().st_size, "sha256": sha256(path)}
             for path in (lib / "global").rglob("*") if path.is_file()}
    portable = {name: {**{k: v for k, v in item.items() if k not in ("inputs", "inputs_tsv")},
                       "inputs": [{k: v for k, v in row.items() if k != "path"} for row in item["inputs"]]}
                for name, item in groups.items()}
    context = {"groups": portable, "global_files": files, "plan_identity": plan_id,
               "source_identity": fingerprint({"plan_identity": plan_id, "prepared_identities": receipt_ids})}
    catalog = {"schema": "metag-partitioned-inputs-v1", "status": "complete", "context": context,
               "identity": fingerprint(context), "groups": groups, "taxonomy_dir": str(tax),
               "seqid_map": str(lib / "global" / "seqid2taxid.map")}
    atomic_json(lib / "catalog.json", catalog)
    if native:
        target = tools / "kraken" / "bin"
        target.mkdir(parents=True)
        for name in ("estimate_capacity", "build_db", "classify", "k2mask"):
            shutil.copy2(Path(native) / name, target / name)
        atomic_json(target.parent / "install.json", {"source_tag": "2.17.2", "source_commit": SOURCE_COMMIT,
                    "files": {path.name: sha256(path) for path in target.iterdir()}})
    cfg = {"root": str(root), "library_dir": str(lib), "metadata_dir": str(meta), "reference_raw": str(raw),
           "work_dir": str(work), "tools_dir": str(tools), "storage_layout": "separated",
           "min_free_gib": 0, "group_build_threads": 2, "group_build_memory_gib": 192,
           "database_layout": "single", "classify_threads": 2, "confidence": 0.1, "minimum_hit_groups": 2}
    return cfg, catalog, sequences, shared


class GroupLibraryTests(unittest.TestCase):
    def test_catalog_validation_does_not_hash_prepared_payloads(self):
        with tempfile.TemporaryDirectory() as temporary:
            cfg, catalog, _, _ = fixture(temporary)
            with patch("group_library.sha256", side_effect=AssertionError("Unexpected SHA256 pass")):
                seen, sources, metadata = group.source_catalog(cfg)
            self.assertEqual(seen["identity"], catalog["identity"])
            self.assertEqual(sum(row["records"] for row in sources), 5)
            self.assertEqual(len(metadata), 4)
            source = Path(sources[0]["path"])
            source.write_bytes(source.read_bytes() + b"A")
            with self.assertRaisesRegex(RuntimeError, "size differs"):
                group.source_catalog(cfg)

    def test_reference_table_uses_frozen_taxonomy_and_complete_source_set(self):
        with tempfile.TemporaryDirectory() as temporary:
            cfg, catalog, _, _ = fixture(temporary)
            # Metadata taxonomy is deliberately unusable: the builder must use
            # global/taxonomy belonging to the already prepared library.
            (Path(cfg["metadata_dir"]) / "taxonomy" / "names.dmp").write_text("")
            out = Path(temporary) / "references.tsv"
            self.assertEqual(group._reference_table(cfg, catalog, out), 3)
            with out.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["group"] for row in rows}, {"human", "plant", "bacteria"})
            self.assertTrue(all(Path(row["path"]).is_absolute() for row in rows))
            self.assertEqual({row["taxid"] for row in rows}, {"9606", "3702", "562"})

    def test_stream_separates_files_without_final_newline(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            files = []
            for i in range(2):
                path = root / f"seq{i}.fna"
                path.write_text(f">seq{i}\nACGT")
                files.append({"path": str(path), "observation": group._observation(path)})
            output = root / "copied"
            group._run_stream([sys.executable, "-c", "import sys; sys.stdout.buffer.write(sys.stdin.buffer.read())"],
                              files, root / "log", output)
            self.assertEqual(output.read_text(), ">seq0\nACGT\n>seq1\nACGT\n")

    def test_stream_failure_is_reported_and_no_stage_is_published(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.write_bytes(b"A" * (1 << 20))
            with self.assertRaisesRegex(RuntimeError, "Native group command"):
                group._run_stream([sys.executable, "-c", "import sys; sys.stderr.write('fixture failure'); sys.exit(3)"],
                    [{"path": str(source), "observation": group._observation(source)}], root / "log", root / "out")
            self.assertIn("fixture failure", (root / "log").read_text())

    def test_protected_output_is_rejected_before_native_access(self):
        with tempfile.TemporaryDirectory() as temporary:
            cfg, _, _, _ = fixture(temporary)
            with self.assertRaisesRegex(RuntimeError, "must be separate"):
                group.build_group(cfg, out=Path(cfg["library_dir"]) / "index" / "oops")

    @unittest.skipUnless(os.environ.get("METAG_TEST_KRAKEN_BIN"), "Set METAG_TEST_KRAKEN_BIN for real pinned native integration")
    def test_failed_build_resumes_estimate_and_replaces_only_owned_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            cfg, _, _, _ = fixture(temporary, os.environ["METAG_TEST_KRAKEN_BIN"])
            real_stream = group._run_stream
            invocations = []
            def stop_before_build(command, *args):
                invocations.append(Path(command[0]).name)
                if Path(command[0]).name == "build_db":
                    raise RuntimeError("Fixture interrupted before build")
                return real_stream(command, *args)
            with patch("group_library._run_stream", side_effect=stop_before_build):
                with self.assertRaisesRegex(RuntimeError, "Fixture interrupted"):
                    group.build_group(cfg, max_gib=48 * 1024 / GIB)
            self.assertFalse(group.group_directory(cfg).exists())
            self.assertEqual(invocations, ["estimate_capacity", "build_db"])
            def resume_only_build(command, *args):
                self.assertEqual(Path(command[0]).name, "build_db", "Saved estimate must be reused")
                return real_stream(command, *args)
            with patch("group_library._run_stream", side_effect=resume_only_build):
                first = group.build_group(cfg, max_gib=48 * 1024 / GIB)
            with self.assertRaisesRegex(RuntimeError, "settings/inputs changed"):
                group.build_group(cfg, max_gib=64 * 1024 / GIB)
            self.assertEqual(group.compact_snapshot(cfg)["identity"], first["identity"])
            second = group.build_group(cfg, max_gib=64 * 1024 / GIB, replace=True)
            self.assertNotEqual(first["identity"], second["identity"])
            self.assertFalse(group.group_directory(cfg).with_name(".group.previous").exists())
            # --replace can repair this workflow's damaged output; it must
            # still never remove an unrelated directory at an arbitrary --out.
            (Path(second["db_dir"]) / "hash.k2d").write_bytes(b"broken")
            repaired = group.build_group(cfg, max_gib=64 * 1024 / GIB, replace=True)
            self.assertEqual(repaired["downsampling"]["capacity_fraction"], second["downsampling"]["capacity_fraction"])
            foreign = Path(temporary) / "foreign"
            foreign.mkdir()
            (foreign / "precious").write_text("keep")
            with self.assertRaises(RuntimeError):
                group.build_group(cfg, out=foreign, max_gib=64 * 1024 / GIB, replace=True)
            self.assertEqual((foreign / "precious").read_text(), "keep")

    @unittest.skipUnless(os.environ.get("METAG_TEST_KRAKEN_BIN"), "Set METAG_TEST_KRAKEN_BIN for real pinned native integration")
    def test_real_native_estimate_cap_competition_and_snapshot(self):
        with tempfile.TemporaryDirectory() as temporary:
            cfg, catalog, seq, shared = fixture(temporary, os.environ["METAG_TEST_KRAKEN_BIN"])
            cap_bytes = 48 * 1024
            # Any attempt to hash a large payload is a test failure. Small
            # native executables and reference metadata may still be checked.
            real_sha = group.sha256
            def guarded_sha(path):
                self.assertNotIn(Path(path).name, {"hash.k2d", "sequences.fna", "seqid2taxid.map"})
                return real_sha(path)
            with patch("group_library.sha256", side_effect=guarded_sha):
                snapshot = group.build_group(cfg, max_gib=cap_bytes / GIB)
            self.assertLessEqual((Path(snapshot["db_dir"]) / "hash.k2d").stat().st_size, cap_bytes)
            self.assertTrue(snapshot["downsampling"]["enabled"])
            self.assertGreater(snapshot["downsampling"]["native_options"]["minimum_acceptable_hash_value"], 0)
            self.assertEqual(snapshot["context"]["processed_records"], 5)
            self.assertEqual(snapshot["context"]["processed_bases"], 76000)
            self.assertEqual(snapshot["preparation_identity"], catalog["identity"])
            self.assertIn("NATIVE_METRICS", (Path(snapshot["db_dir"]) / "build.log").read_text())
            with patch("group_library._run_stream", side_effect=AssertionError("Completed build should reuse")):
                self.assertEqual(group.build_group(cfg, max_gib=cap_bytes / GIB)["identity"], snapshot["identity"])

            # Both ends come from one real fragment: 150 bp FR reads, insert 350.
            reads = [(9606, seq[9606][500:850]), (3702, seq[3702][500:850]),
                     (562, seq[562][500:850]), (2759, shared[500:850]), (0, "N" * 350)]
            reverse = lambda value: value.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            sample = {"sample_id": "fixture/group", "sample": "group", "dataset": "fixture"}
            for mate in (1, 2):
                path = Path(temporary) / f"reads.R{mate}.fastq"
                with path.open("w") as handle:
                    for i, (_, sequence) in enumerate(reads):
                        dna = sequence[:150] if mate == 1 else reverse(sequence[-150:])
                        handle.write(f"@fragment{i}/{mate}\n{dna}\n+\n" + "I" * 150 + "\n")
                sample[f"r{mate}"] = str(path)
            query_dir = Path(temporary) / "query"
            query_dir.mkdir()
            query = query_many(cfg, [sample], snapshot, query_dir)[0]
            assigned = struct.unpack("<" + "Q" * len(reads), Path(query["taxids"]).read_bytes())
            self.assertEqual(assigned, tuple(row[0] for row in reads))
            self.assertEqual(query["memory_mode"], "mmap")
            self.assertNotIn("--metag-hashes", query["command"])

            uncapped = group.build_group(cfg, out=Path(cfg["library_dir"]) / "group-uncapped", max_gib=0.01)
            self.assertFalse(uncapped["downsampling"]["enabled"])
            self.assertEqual(uncapped["downsampling"]["native_options"]["minimum_acceptable_hash_value"], 0)
            self.assertEqual(uncapped["context"]["processed_bases"], 76000)
            self.assertEqual(uncapped["downsampling"]["native_options"]["selected_hash_fraction"], 1)
            # Same-length corruption is detected by observation, without SHA.
            options = Path(snapshot["db_dir"]) / "opts.k2d"
            options.write_bytes(options.read_bytes())
            with self.assertRaisesRegex(RuntimeError, "changed"):
                group.compact_snapshot(cfg)


if __name__ == "__main__":
    unittest.main(verbosity=2)
