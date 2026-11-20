#!/usr/bin/env python3
"""Actual native group/species routing with byte- and ordinal-level assertions.

Run with the already installed pinned Kraken2 executable directory:
  METAG_TEST_KRAKEN_BIN=/path/to/kraken/bin python3 f/test_hierarchy_native.py

Only the full-catalog provider and unrelated release code fingerprint are test
fixtures. Native estimation, capped/group and uncapped/category builds, every
native classification, FASTQ parsing, compression and routing execute normally.
"""
from collections import Counter
import csv
import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from common import atomic_json, fingerprint
from fastq_stream import new_digest, read_id, update_digest
import group_library
import hierarchical_split as hierarchy
from partitioned_build import record
from reference import GIB, PARAMETERS, _native_index_structure, _native_opts
from test_group_library import fixture
from taxonomy import Taxonomy


def read_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_detailed_categories(cfg, catalog, compact):
    """Create real category-native files and their actual completion identities."""
    root = Path(cfg["library_dir"]) / "index"
    root.mkdir()
    capacity = compact["context"]["resource_estimate"]["capacity_before_cap"]
    identities = {}
    for category, prepared in catalog["groups"].items():
        folder = root / category
        folder.mkdir()
        command = [Path(cfg["tools_dir"]) / "kraken/bin/build_db",
                   "-k", "35", "-l", "31", "-S", PARAMETERS["seed"], "-p", "2",
                   "-H", folder / "hash.k2d", "-t", folder / "taxo.k2d", "-o", folder / "opts.k2d",
                   "-n", catalog["taxonomy_dir"], "-m", catalog["seqid_map"],
                   "-c", str(capacity), "-B", "16384", "-b", "8192", "-r", "0"]
        sources = [{"path": row["path"], "observation": group_library._observation(row["path"])}
                   for row in prepared["inputs"]]
        group_library._run_stream(command, sources, folder / "build.log", folder / "build.stdout")
        _native_index_structure(folder, expected_capacity=capacity)
        options = _native_opts(folder / "opts.k2d")
        context = {"category": category, "db_dir": str(folder), "build_identity": "native-test-" + category,
                   "preparation_identity": catalog["identity"], "parameters": PARAMETERS,
                   "files": {name: record(folder / name) for name in group_library.NATIVE_FILES},
                   "resource": {"capacity": capacity}, "native_options": options}
        identities[category] = fingerprint(context)
        atomic_json(folder / "database.json", {"schema": "metag-partitioned-index-v1", "status": "complete",
                    "identity": identities[category], "context": context})
    # Only the publication catalog provider is replaced: the single-index
    # validator still reads every category's real database.json and native file.
    full = {"identity": "fixture-full-catalog", "status": "complete", "layout": "partitioned",
            "db_dir": str(root), "taxonomy_dir": catalog["taxonomy_dir"],
            "taxonomy_files": compact["taxonomy_files"],
            "context": {"categories": list(catalog["groups"])}}
    return full, identities


def committed_flat_fixture(root):
    """A minimal real, size-committed original flat split; no classifier needed."""
    cfg, catalog, sequences, shared = fixture(root)
    cfg.update(output_dir=str(Path(root) / "species"), split_gzip_level=3,
               evaluation_pairs_per_stratum=100, evaluation_target_pairs_per_stratum=100,
               taxonomy_dir=str(Path(root) / "wrong-taxonomy"), minimum_base_quality=0)
    sample = {"sample_id": "fixture/original", "sample": "original", "dataset": "fixture"}
    tax = Taxonomy(catalog["taxonomy_dir"])
    truth = [9606, 3702, 2, 2759, 0, 3702]
    original, assignments, counts, statuses, groups = [], [], Counter(), Counter(), Counter()
    digest = new_digest()
    destination = Path(cfg["output_dir"]) / sample["sample_id"]
    qc = destination / "qc"
    qc.mkdir(parents=True)
    bin_contents = {}
    for ordinal, tid in enumerate(truth, 1):
        sequence = (shared if tid == 2759 else sequences.get(tid, sequences[562]))[300:450]
        pair = []
        for mate in (1, 2):
            pair.append(tuple(line.encode() for line in (
                f"@ERR14374050.{ordinal}.{mate} original comment\r\n", sequence + "\r\n",
                "+ untouched plus line\r\n", "I" * 150 + "\r\n")))
        a, b = pair
        original.append((a, b))
        update_digest(digest, a, b)
        info = tax.describe(tid)
        if not tid:
            bin_name = "unclassified"
        elif info["group"] == "unresolved":
            bin_name = "taxon_unresolved"
        elif not info["species_taxid"]:
            bin_name = info["group"] + "/taxon_unresolved"
        else:
            bin_name = f"{info['group']}/{info['species_taxid']}_fixture"
        counts[(bin_name, tid, info["species_taxid"], info["group"], info["status"], info["name"])] += 1
        statuses[info["status"]] += 1
        groups[info["group"]] += 1
        assignments.append({"ordinal": ordinal, "r1_read_id": read_id(a[0]).decode(),
                            "r2_read_id": read_id(b[0]).decode(), "assigned_taxid": tid,
                            # A prior host rescue must retain the final assigned
                            # taxid, rather than reverting to this earlier call.
                            "kraken_taxid": 0 if ordinal == 1 else tid})
        for mate, record_bytes in ((1, a), (2, b)):
            name = f"{bin_name}.R{mate}.fastq.gz"
            bin_contents.setdefault(name, bytearray()).extend(b"".join(record_bytes))
    for mate in (1, 2):
        path = Path(root) / f"original.R{mate}.fastq"
        path.write_bytes(b"".join(b"".join(pair[mate - 1]) for pair in original))
        sample[f"r{mate}"] = str(path)
    for name, data in bin_contents.items():
        path = destination / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(gzip.compress(bytes(data), mtime=0))
    count_rows = [dict(bin=key[0], taxid=key[1], species_taxid=key[2], group=key[3], status=key[4],
                       name=key[5], fragments=amount, reads=2 * amount) for key, amount in counts.items()]
    hierarchy.write_table(qc / "counts.tsv", ("bin", "taxid", "species_taxid", "group", "status", "name", "fragments", "reads"), count_rows)
    by_bin = Counter()
    for row in count_rows:
        by_bin[row["bin"]] += row["fragments"]
    hierarchy.write_table(qc / "bin-counts.tsv", ("bin", "fragments"),
                          [dict(bin=name, fragments=n) for name, n in by_bin.items()])
    with gzip.open(qc / "read-assignments.tsv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, tuple(assignments[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(assignments)
    summary = {"sample_id": sample["sample_id"], "paired": True, "fragments": len(truth), "reads": 2 * len(truth),
               "count_conservation": True, "requested_limit": 0, "validated_scope": "entire_input",
               "statuses": dict(statuses), "groups": dict(groups), "input_digest": digest.hexdigest()}
    atomic_json(qc / "summary.json", summary)
    full = {"identity": "committed-original-full-library", "taxonomy_dir": catalog["taxonomy_dir"]}
    context = {"database_identity": full["identity"], "confidence": 0.25,
               "minimum_hit_groups": 3, "minimum_base_quality": 9}
    atomic_json(qc / "provenance.json", {"context": context, "identity": fingerprint(context)})
    outputs = {str(path.relative_to(destination)): {"state": {"size": path.stat().st_size}}
               for path in destination.rglob("*") if path.is_file()}
    atomic_json(qc / "complete.json", {"status": "complete", "outputs": outputs})
    return cfg, sample, full, truth, original, destination


@unittest.skipUnless(os.environ.get("METAG_TEST_KRAKEN_BIN"),
                     "Set METAG_TEST_KRAKEN_BIN for actual pinned native hierarchy integration")
class NativeHierarchyTests(unittest.TestCase):
    def test_group_then_selected_species_preserves_bytes_ordinals_and_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            cfg, catalog, sequences, shared = fixture(temporary, os.environ["METAG_TEST_KRAKEN_BIN"])
            cfg.update(output_dir=str(Path(temporary) / "species"), split_gzip_level=3,
                       evaluation_pairs_per_stratum=100, evaluation_target_pairs_per_stratum=100,
                       materialize_species_min_fragments=20, species_groups=["human", "plant"],
                       group_confidence=0.0, group_minimum_hit_groups=2,
                       # Deliberately unavailable: labels must come from the
                       # frozen DB taxonomy and the coarse completion record.
                       taxonomy_dir=str(Path(temporary) / "wrong-taxonomy"))
            compact = group_library.build_group(cfg, max_gib=48 * 1024 / GIB)
            full, category_identities = build_detailed_categories(cfg, catalog, compact)
            truth = [9606, 3702, 2759, 562, 9606, 0, 3702, 2759, 562]
            groups = {9606: "human", 3702: "plant", 562: "bacteria", 2759: "ambiguous", 0: "unclassified"}
            original = []
            reverse = lambda value: value.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            for ordinal, taxid in enumerate(truth, 1):
                reference = shared if taxid == 2759 else "N" * 24000 if not taxid else sequences[taxid]
                offset = 450 + ordinal * 50
                fragment = reference[offset:offset + 350]
                pair = []
                newline = "\r\n" if ordinal % 2 else "\n"
                for mate in (1, 2):
                    sequence = fragment[:150] if mate == 1 else reverse(fragment[-150:])
                    # Every record deliberately repeats the same read ID. The
                    # comments, plus-line text and CRLF bytes are also retained.
                    lines = [f"@duplicate/{mate} fixture={ordinal}{newline}", sequence + newline,
                             f"+ original={ordinal}{newline}", "I" * 150 + newline]
                    pair.append(tuple(line.encode() for line in lines))
                original.append(tuple(pair))
            sample = {"sample_id": "fixture/interleaved", "sample": "interleaved", "dataset": "fixture"}
            for mate in (1, 2):
                path = Path(temporary) / f"interleaved.R{mate}.fastq"
                path.write_bytes(b"".join(b"".join(pair[mate - 1]) for pair in original))
                sample[f"r{mate}"] = str(path)
            digest = new_digest()
            for pair in original:
                update_digest(digest, *pair)
            with patch("hierarchical_split.code_identity", return_value="native-fixture-code"):
                root = hierarchy.split_group(cfg, sample)
            coarse = root / "group"
            summary = json.loads((coarse / "qc/summary.json").read_text())
            self.assertEqual(summary["group_counts"], {"human": 2, "plant": 2, "bacteria": 2,
                                                       "ambiguous": 2, "unclassified": 1})
            self.assertEqual([summary[k] for k in ("input_pairs", "written_pairs", "assignment_rows", "expected_pairs")], [9] * 4)
            self.assertEqual(summary["input_digest"], digest.hexdigest())
            self.assertEqual(summary["taxonomy_dir"], catalog["taxonomy_dir"])
            self.assertEqual(summary["request"]["method"], "compact_competitive_kraken")
            rows = read_rows(coarse / "qc/assignments.tsv.gz")
            self.assertEqual([int(row["ordinal"]) for row in rows], list(range(1, 10)))
            self.assertEqual([int(row["taxid"]) for row in rows], truth)
            self.assertEqual([row["group"] for row in rows], [groups[taxid] for taxid in truth])
            before = {}
            for group in set(groups.values()):
                selected = [i for i, tid in enumerate(truth) if groups[tid] == group]
                indexes = read_rows(coarse / group / "ordinals.tsv.gz")
                self.assertEqual([int(row["ordinal"]) for row in indexes], [i + 1 for i in selected])
                for mate in (1, 2):
                    path = coarse / group / f"R{mate}.fastq.gz"
                    expected = b"".join(b"".join(original[i][mate - 1]) for i in selected)
                    before[str(path)] = gzip.decompress(path.read_bytes())
                    self.assertEqual(before[str(path)], expected)
            self.assertIn(str(Path(compact["db_dir"]) / "hash.k2d"), summary["queries"][0]["command"])
            self.assertNotIn("--metag-hashes", summary["queries"][0]["command"])

            with patch("hierarchical_split.code_identity", return_value="native-fixture-code"), \
                 patch("partitioned_build.database_snapshot", return_value=full):
                # Exercise the real category identity consumed by the ordinary
                # native query validator before running selected species.
                for category in ("plants", "vertebrates"):
                    detailed = hierarchy.category_snapshot(cfg, category)
                    self.assertEqual(detailed["identity"], category_identities[category])
                hierarchy.split_species(cfg, sample)
            detail = json.loads((coarse / "species.qc/summary.json").read_text())
            self.assertEqual(detail["processed_groups"], ["human", "plant"])
            self.assertEqual(detail["scope"], "selected_coarse_groups_only")
            self.assertEqual([detail[k] for k in ("input_pairs", "written_pairs", "assignment_rows")], [4] * 3)
            self.assertEqual(detail["expected_pairs"], 9)
            self.assertEqual(detail["unprocessed_pairs"], 5)
            species_rows = read_rows(coarse / "species.qc/assignments.tsv.gz")
            self.assertEqual([int(row["ordinal"]) for row in species_rows], [1, 5, 2, 7])
            self.assertEqual([int(row["species_taxid"]) for row in species_rows], [9606, 9606, 3702, 3702])
            self.assertTrue(all(row["status"] == "assigned_species" for row in species_rows))
            # Small plant bins are grouped under 'other' without erasing their
            # species labels; human remains explicitly materialized.
            for group, leaf in (("human", "9606"), ("plant", "other")):
                for mate in (1, 2):
                    parent = coarse / group / f"R{mate}.fastq.gz"
                    child = coarse / group / "species" / leaf / f"R{mate}.fastq.gz"
                    self.assertEqual(gzip.decompress(child.read_bytes()), before[str(parent)])
            for path, data in before.items():
                self.assertEqual(gzip.decompress(Path(path).read_bytes()), data)
            for group in ("ambiguous", "unclassified", "bacteria"):
                self.assertFalse((coarse / group / "species").exists())
            self.assertEqual([int(row["ordinal"]) for row in read_rows(coarse / "species.qc/sample.tsv")], [1, 2, 5, 7])
            for query, category in zip(detail["queries"], ("vertebrates", "plants"), strict=True):
                self.assertIn(str(Path(cfg["library_dir"]) / "index" / category / "hash.k2d"), query["command"])
                self.assertNotIn("--metag-hashes", query["command"])
                self.assertEqual(query["memory_mode"], "mmap")


if __name__ == "__main__":
    unittest.main(verbosity=2)
