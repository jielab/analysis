#!/usr/bin/env python3
"""Native split plus the automatic full evaluation callback, on tiny references.

Required optional integration environment:
  METAG_TEST_KRAKEN_BIN=/path/to/pinned/kraken/bin
  METAG_TEST_KRAKEN_SOURCE=/path/to/pinned/kraken/src
  MINIMAP2_TEST_BINARY=/path/to/minimap2
Optional: METAG_TEST_EVALUATION_ARTIFACTS=/directory/to/retain/evidence

No classifier, mapper, evaluation function or result is mocked. CLI config and
sample discovery, the unrelated release fingerprint, and the prepared-catalog
loader during fixture publication are substituted. All native file validators,
real full-union/catalog queries, alignment, controls, PCA and reports execute.
"""
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import zipfile

from common import atomic_json
from common import fingerprint
from fastq_stream import fragments
import group_library
import hierarchical_split as hierarchy
import native_union
import partitioned_build
from reference import GIB
from test_group_library import fixture
from test_hierarchy_native import build_detailed_categories, read_rows
from test_stage_alignment import minimap_binary


def install_test_union(cfg, source, destination):
    """Build the actual adapter against pinned source, with no network/install."""
    source, destination = Path(source), Path(destination)
    destination.mkdir()
    for path in source.iterdir():
        if path.is_file() and (path.name == "Makefile" or path.suffix in (".cc", ".h", ".o")):
            shutil.copy2(path, destination / path.name)
    (destination / "metag_classify.cc").write_text(
        native_union.patched_classifier((destination / "classify.cc").read_text()))
    shutil.copy2(Path(native_union.__file__).with_name("metag_union.h"), destination / "metag_union.h")
    with (destination / "Makefile").open("a") as makefile:
        makefile.write(native_union.TARGET)
    with (destination / "compile.log").open("w") as log:
        subprocess.run(["make", "-j2", "metag_classify"], cwd=destination,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    target = native_union.executable_path(cfg)
    target.parent.mkdir(parents=True)
    shutil.copy2(destination / "metag_classify", target)
    return target


def normalize_fixture_category_padding(cfg, catalog):
    """Match production metag_build's explicitly zeroed IndexOptions padding.

The fixture uses official build_db directly, whose unused C++ struct padding
is unspecified. Production category construction already zero-initializes it.
Only those 11 unused bytes change; all defined fields and taxonomy must agree.
"""
    options = struct.Struct("<QQQQ?7xQiii4x")
    first_fields = first_taxonomy = None
    for category in catalog["groups"]:
        folder = Path(cfg["library_dir"]) / "index" / category
        path = folder / "opts.k2d"
        raw = path.read_bytes()
        fields = options.unpack(raw)
        taxonomy = (folder / "taxo.k2d").read_bytes()
        if first_fields is None:
            first_fields, first_taxonomy = fields, taxonomy
        elif fields != first_fields or taxonomy != first_taxonomy:
            raise AssertionError("Fixture native category fields/taxonomy disagree")
        path.write_bytes(options.pack(*fields))
        marker = json.loads((folder / "database.json").read_text())
        marker["context"]["files"]["opts.k2d"] = partitioned_build.record(path)
        marker["identity"] = fingerprint(marker["context"])
        atomic_json(folder / "database.json", marker)


def integration_available():
    source = os.environ.get("METAG_TEST_KRAKEN_SOURCE")
    return bool(os.environ.get("METAG_TEST_KRAKEN_BIN") and source and
                (Path(source) / "classify.cc").is_file() and minimap_binary())


@unittest.skipUnless(integration_available(), "Set native Kraken source/bin and minimap2 integration paths")
class NativeAutomaticEvaluationTests(unittest.TestCase):
    def test_automatic_group_evaluation_controls_pca_alignment_and_zip(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cfg, catalog, sequences, shared = fixture(root, os.environ["METAG_TEST_KRAKEN_BIN"])
            cfg.update(output_dir=str(root / "species"), split_gzip_level=3, minimum_base_quality=0,
                       taxonomy_dir=catalog["taxonomy_dir"], group_confidence=0.0, group_minimum_hit_groups=2,
                       evaluation_pairs_per_stratum=100, evaluation_target_pairs_per_stratum=100,
                       evaluation_agreement_pairs_per_stratum=100, evaluation_anchor_references=1,
                       evaluation_anchor_bases=8000, evaluation_anchor_pairs=8,
                       evaluation_seed=20261115, species_groups=["human", "plant"],
                       stage_alignment={"minimap2": minimap_binary(), "threads": 2,
                                        "max_reference_gib": .001, "max_references": 6,
                                        "offtarget_per_stratum": 100})
            sample = {"sample_id": "fixture/evaluate", "sample": "evaluate", "dataset": "fixture"}
            sample_root = Path(cfg["output_dir"]) / sample["sample_id"]
            qc = sample_root / "group/qc"
            artifact_root = os.environ.get("METAG_TEST_EVALUATION_ARTIFACTS")
            try:
                compact = group_library.build_group(cfg, max_gib=48 * 1024 / GIB)
                build_detailed_categories(cfg, catalog, compact)
                normalize_fixture_category_padding(cfg, catalog)
                # The fixture already has sealed prepared metadata and native
                # category receipts; publish and validate the real union catalog.
                with patch("partitioned_build._preparation", return_value=catalog):
                    full = partitioned_build.publish(cfg)
                self.assertEqual(partitioned_build.database_snapshot(cfg)["identity"], full["identity"])
                union_binary = install_test_union(cfg, os.environ["METAG_TEST_KRAKEN_SOURCE"], root / "union-source")
                truth = [9606, 3702, 562, 9606, 3702, 562, 2759, 0]
                reverse = lambda sequence: sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                for mate in (1, 2):
                    path = root / f"evaluate.R{mate}.fastq"
                    with path.open("w") as handle:
                        for ordinal, tid in enumerate(truth, 1):
                            reference = shared if tid == 2759 else "N" * 24000 if not tid else sequences[tid]
                            start = 500 + ordinal * 60
                            fragment = reference[start:start + 350]
                            read = fragment[:150] if mate == 1 else reverse(fragment[-150:])
                            handle.write(f"@evaluation{ordinal}/{mate}\n{read}\n+\n" + "I" * 150 + "\n")
                    sample[f"r{mate}"] = str(path)
                truth_file = root / "truth.tsv"
                hierarchy.write_table(truth_file, ("ordinal", "read_id", "true_taxid", "true_group"),
                    [dict(ordinal=i, read_id="evaluation" + str(i), true_taxid=tid,
                          true_group={9606:"human",3702:"plant",562:"bacteria",2759:"ambiguous",0:"unknown"}[tid])
                     for i, tid in enumerate(truth, 1)])
                # Call main, not evaluate directly: this verifies the explicit
                # split/evaluate sequence in one worker and cohort publication.
                with patch("hierarchical_split.load_config", return_value=cfg), \
                     patch("hierarchical_split.discover_samples", return_value=[sample]), \
                     patch("hierarchical_split.code_identity", return_value="native-evaluation-fixture"):
                    hierarchy.main(["group,evaluate", "--evaluation", "full", "--stage", "group", "--allow-local", "--threads", "2", "--sample", sample["sample_id"],
                                    "--truth", str(truth_file)])
                summary = json.loads((qc / "summary.json").read_text())
                evaluation = json.loads((qc / "evaluate.json").read_text())
                report = json.loads((qc / "stage-evaluation.json").read_text())
                self.assertEqual(evaluation["status"], "complete", evaluation)
                self.assertEqual(summary["evaluation_status"], "complete")
                self.assertEqual(report["status"], "complete", report)
                self.assertEqual(summary["input_pairs"], 8)
                self.assertEqual(report["accounting"]["status"], "complete")
                self.assertEqual(report["sample_qc"]["status"], "complete")

                agreement = evaluation["agreement"]
                self.assertEqual(agreement["status"], "complete", agreement)
                self.assertEqual(agreement["pairs"], 8)
                self.assertEqual(agreement["full_database"], full["identity"])
                self.assertEqual(agreement["query"]["native_executable"], str(union_binary))
                self.assertIn("--metag-hashes", agreement["query"]["command"])
                self.assertEqual(agreement["query"]["memory_mode"], "mmap")

                alignment = evaluation["alignment"]
                self.assertEqual(alignment["status"], "complete", alignment)
                self.assertEqual(alignment["tested_pairs"], 8)
                self.assertEqual(alignment["mapper_version"], "2.28-r1209")
                self.assertGreater(len(alignment["commands"]), 0)
                self.assertFalse(alignment["accuracy_estimated"])

                controls = report["controls"]
                self.assertEqual(report["anchors"]["anchor_pairs"], 24)
                self.assertEqual(report["anchors"]["sampling_scope"], "bounded_reference_prefix_windows")
                self.assertEqual(report["controls_execution"]["status"], "complete")
                self.assertEqual(report["controls_execution"]["pairs"], 24)
                self.assertEqual(controls["status"], "complete", controls)
                self.assertEqual(controls["truth_rows"], 24)
                self.assertEqual(controls["assignment_rows"], 24)
                self.assertGreater(controls["correct_pairs"], 0)
                self.assertTrue(0 <= controls["recall"] <= 1)
                self.assertEqual(controls["control_kind"], "reference_prefix_positive_control")
                self.assertEqual(len(read_rows(qc / "controls.assignments.tsv.gz")), 24)
                self.assertEqual(sum(1 for _ in fragments(qc / "controls.R1.fastq.gz", qc / "controls.R2.fastq.gz")), 24)
                for query in report["controls_execution"]["queries"]:
                    self.assertIn(str(Path(compact["db_dir"]) / "hash.k2d"), query["command"])

                self.assertEqual(report["composition"]["status"], "complete")
                self.assertEqual(report["composition"]["fit_scope"], "reference_anchors_only")
                self.assertEqual(report["composition"]["reference_anchor_points"], 24)
                self.assertTrue((qc / "composition.svg").read_text().lstrip().startswith("<svg"))
                self.assertEqual(report["truth"]["known_label_rows"], 6)
                self.assertEqual(report["truth"]["unknown_truth_rows"], 2)
                self.assertEqual(report["truth"]["correct_pairs"], 6)
                html = (qc / "evaluation.html").read_text()
                self.assertIn("composition.svg", html)
                self.assertIn("diagnostic.zip", html)
                self.assertNotIn("generated_not_run", html)
                link = Path(cfg["work_dir"]) / "qc/evaluate" / sample["sample_id"] / "group"
                self.assertTrue(link.is_symlink())
                self.assertEqual(link.resolve(), qc.resolve())
                self.assertTrue((Path(cfg["work_dir"]) / "qc/hierarchy/report.html").is_file())

                with zipfile.ZipFile(qc / "diagnostic.zip") as archive:
                    names = set(archive.namelist())
                    self.assertTrue({"sample.tsv", "sample.R1.fastq.gz", "sample.R2.fastq.gz",
                                     "composition.svg", "composition.tsv", "evaluation.html", "stage-evaluation.json",
                                     "controls-evaluation.json", "alignment/summary.json", "agreement/summary.json"} <= names)
                    self.assertFalse(any("/species/" in name or "hash.k2d" in name for name in names))
                    bundled = json.loads(archive.read("stage-evaluation.json"))
                    self.assertEqual(bundled["controls"]["status"], "complete")
                    self.assertEqual(bundled["status"], "complete")
                    self.assertEqual(bundled, report, "Bundled and live evaluation reports must be identical")
                    bundled_summary = json.loads(archive.read("summary.json"))
                    self.assertEqual(bundled_summary["evaluation_status"], summary["evaluation_status"],
                                     "The diagnostic ZIP must contain final stage evaluation status")
                    self.assertNotIn("diagnostic-bundle.json", names, "ZIP byte-count metadata must not be self-referential")
                bundle = json.loads((qc / "diagnostic-bundle.json").read_text())
                self.assertEqual(bundle["bytes"], (qc / "diagnostic.zip").stat().st_size)

                # Exercise the optional rescue branch using the same actual
                # native union. Unknown parents here are shared DNA and all-N;
                # they must remain unresolved instead of gaining false species.
                coarse_bytes = {str(path): gzip.decompress(path.read_bytes())
                                for group in ("human", "plant", "bacteria", "ambiguous", "unclassified")
                                for path in (sample_root / "group" / group).glob("R?.fastq.gz")}
                with patch("hierarchical_split.load_config", return_value=cfg), \
                     patch("hierarchical_split.discover_samples", return_value=[sample]), \
                     patch("hierarchical_split.code_identity", return_value="native-evaluation-fixture"):
                    hierarchy.main(["species,evaluate", "--evaluation", "full", "--stage", "species", "--allow-local", "--threads", "2", "--sample", sample["sample_id"],
                                    "--groups", "human,plant", "--include-uncertain", "--truth", str(truth_file)])
                species_qc = sample_root / "group/species.qc"
                species_summary = json.loads((species_qc / "summary.json").read_text())
                species_evaluation = json.loads((species_qc / "evaluate.json").read_text())
                species_report = json.loads((species_qc / "stage-evaluation.json").read_text())
                self.assertEqual(species_summary["processed_groups"], ["human", "plant", "ambiguous", "unclassified"])
                self.assertEqual(species_summary["selected_groups"], ["human", "plant"])
                self.assertTrue(species_summary["request"]["include_uncertain"])
                self.assertEqual(species_summary["scope"], "selected_groups_and_uncertain")
                self.assertEqual([species_summary[k] for k in ("input_pairs", "written_pairs", "assignment_rows")], [6] * 3)
                self.assertEqual(species_summary["unprocessed_pairs"], 2)
                self.assertEqual(species_summary["expected_pairs"], 8)
                self.assertEqual(species_summary["input_pairs"] + species_summary["unprocessed_pairs"], 8)
                species_rows = read_rows(species_qc / "assignments.tsv.gz")
                ordinals = [int(row["ordinal"]) for row in species_rows]
                self.assertEqual(ordinals, [1, 4, 2, 5, 7, 8])
                self.assertEqual(len(set(ordinals)), 6)
                self.assertEqual({int(row["ordinal"]) for row in species_rows if row["source_group"] == "ambiguous"}, {7})
                self.assertEqual({int(row["ordinal"]) for row in species_rows if row["source_group"] == "unclassified"}, {8})
                for row in species_rows:
                    if row["source_group"] in ("ambiguous", "unclassified"):
                        self.assertEqual(row["species_taxid"], "0")
                        self.assertEqual(row["group"], row["source_group"])
                for group in ("ambiguous", "unclassified"):
                    for mate in (1, 2):
                        parent = sample_root / "group" / group / f"R{mate}.fastq.gz"
                        child = parent.parent / "species/unresolved" / parent.name
                        self.assertEqual(gzip.decompress(child.read_bytes()), coarse_bytes[str(parent)])
                self.assertFalse((sample_root / "group/bacteria/species").exists())
                for path, data in coarse_bytes.items():
                    self.assertEqual(gzip.decompress(Path(path).read_bytes()), data)
                for query in species_summary["queries"][2:]:
                    self.assertIn("--metag-hashes", query["command"])
                    self.assertEqual(query["native_executable"], str(union_binary))
                    self.assertEqual(query["fragments"], 1)
                self.assertEqual(species_summary["evaluation_status"], "complete")
                self.assertEqual(species_evaluation["status"], "complete", species_evaluation)
                self.assertEqual(species_report["status"], "complete", species_report)
                self.assertEqual(species_report["accounting"]["status"], "complete")
                self.assertTrue(all(check["passed"] for check in species_report["accounting"]["checks"]))
                self.assertEqual(species_evaluation["agreement"]["pairs"], 6)
                self.assertEqual(species_evaluation["alignment"]["tested_pairs"], 6)
                self.assertEqual(species_report["controls_execution"]["pairs"], 16)
                self.assertEqual(species_report["controls"]["control_scope"], "conditional_on_known_source_group")
                self.assertEqual(species_report["controls"]["correct_pairs"], 16)
                self.assertEqual(species_report["truth"]["missing_predictions"], 2)
                self.assertEqual(species_report["truth"]["selected_group_scope"]["true_pairs"], 4)
                self.assertEqual(species_report["truth"]["selected_group_scope"]["correct_pairs"], 4)
                self.assertEqual(species_report["truth"]["selected_group_scope"]["recall"], 1.0)
                self.assertEqual(species_report["composition"]["status"], "complete")
                with zipfile.ZipFile(species_qc / "diagnostic.zip") as archive:
                    self.assertEqual(json.loads(archive.read("stage-evaluation.json")), species_report)
                    self.assertEqual(json.loads(archive.read("summary.json")), species_summary)
                    self.assertEqual(json.loads(archive.read("controls-evaluation.json"))["control_scope"],
                                     "conditional_on_known_source_group")
                species_bundle = json.loads((species_qc / "diagnostic-bundle.json").read_text())
                self.assertEqual(species_bundle["bytes"], (species_qc / "diagnostic.zip").stat().st_size)
                atomic_json(qc / "native-integration.json", {"status": "passed", "input_pairs": 8,
                            "control_pairs": 24, "union_source_commit": "8c190b1b668825935dbf6dee5f969227dc8269bb",
                            "minimap2_version": alignment["mapper_version"], "controls_recall": controls["recall"],
                            "species_processed_pairs": 6, "species_unprocessed_pairs": 2,
                            "species_uncertain_parents": ["ambiguous", "unclassified"],
                            "species_control_pairs": 16, "species_callback_status": species_evaluation["status"],
                            "fixture_boundaries": ["CLI config/sample discovery", "release code fingerprint",
                                                   "prepared catalog loader during real native category publication",
                                                   "unused IndexOptions padding zeroed, as production metag_build does"]})
            finally:
                if artifact_root and qc.is_dir():
                    destination = Path(artifact_root)
                    destination.mkdir(parents=True, exist_ok=True)
                    shutil.copytree(qc, destination / "qc", dirs_exist_ok=True)
                    species_qc = sample_root / "group/species.qc"
                    if species_qc.is_dir():
                        shutil.copytree(species_qc, destination / "species.qc", dirs_exist_ok=True)
                    if (root / "union-source/compile.log").is_file():
                        shutil.copy2(root / "union-source/compile.log", destination / "union-compile.log")


if __name__ == "__main__":
    unittest.main(verbosity=2)
