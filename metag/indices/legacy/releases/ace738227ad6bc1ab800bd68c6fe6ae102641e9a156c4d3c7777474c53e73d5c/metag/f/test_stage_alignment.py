#!/usr/bin/env python3
"""Meaningful small-fixture checks for sampled stage alignment evaluation.

Native mapper tests run when minimap2 is installed or the local integration-test
binary is present. They use synthetic genomes/reads, not the user's real reads.
"""
import csv
import gzip
import json
import os
from pathlib import Path
import random
import shutil
import tempfile
import unittest

from stage_alignment import evaluate_alignment


def table(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def minimap_binary():
    value = os.environ.get("MINIMAP2_TEST_BINARY") or shutil.which("minimap2")
    if value:
        return value
    base = Path(__file__).resolve().parents[3]
    candidates = sorted((base / "native-alignment-test").glob("minimap2*/minimap2"))
    return str(candidates[0]) if candidates else None


class StageAlignmentTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.sample = self.base / "sample"
        self.qc = self.sample / "group/qc"
        self.qc.mkdir(parents=True)
        self.tax = self.base / "metadata/taxonomy"
        self.tax.mkdir(parents=True)
        nodes = [(1, 1, "no rank"), (2759, 1, "superkingdom"), (33090, 2759, "kingdom"),
                 (100, 33090, "genus"), (101, 100, "species"), (102, 100, "species"),
                 (103, 100, "species"), (7742, 2759, "subphylum"), (40674, 7742, "class"),
                 (9604, 40674, "family"), (9605, 9604, "genus"), (9606, 9605, "species"),
                 (9596, 9604, "genus"), (9598, 9596, "species"), (9595, 9604, "species"),
                 (2, 1, "superkingdom"), (200, 2, "genus"), (201, 200, "species")]
        (self.tax / "nodes.dmp").write_text("".join(f"{t}\t|\t{p}\t|\t{rank}\t|\n" for t, p, rank in nodes))
        (self.tax / "names.dmp").write_text("".join(f"{t}\t|\tTest{t}\t|\t\t|\tscientific name\t|\n" for t, _, _ in nodes))
        (self.tax / "merged.dmp").write_text("")
        rng = random.Random(20261115)
        dna = lambda length: "".join(rng.choices("ACGT", k=length))
        self.sequences = {tid: dna(9000) for tid in (9606, 9598, 9595, 101, 102, 201)}
        for first, second in ((9606, 9598), (101, 102)):
            shared = dna(1000)
            for tid in (first, second):
                self.sequences[tid] = self.sequences[tid][:6000] + shared + self.sequences[tid][7000:]
        refs = self.base / "references"
        refs.mkdir()
        representatives = []
        for tid, sequence in self.sequences.items():
            path = refs / (str(tid) + ".fa")
            path.write_text(f">chromosome1 synthetic test genome {tid}\n{sequence}\n")
            group = "plant" if tid in (101, 102) else "bacteria" if tid == 201 else "human" if tid == 9606 else "mammal"
            representatives.append(dict(species_taxid=tid, name="Test" + str(tid), group=group,
                                        path=str(path), bases=len(sequence)))
        self.representatives = self.base / "representatives.tsv"
        table(self.representatives, representatives)
        def rc(sequence):
            return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        def pair(tid, start):
            return self.sequences[tid][start:start + 150], rc(self.sequences[tid][start + 300:start + 450])
        self.cases = [
            ("human", 9606, "human_correct", pair(9606, 1000)),
            ("human", 9606, "human_wrong_plant", pair(101, 1500)),
            ("human", 9606, "human_primate_tie", pair(9606, 6200)),
            ("plant", 101, "plant_correct", pair(101, 2500)),
            ("plant", 101, "plant_wrong_human", pair(9606, 2500)),
            ("plant", 101, "plant_relative_tie", pair(101, 6200)),
            ("plant", 101, "plant_other_species", pair(102, 3000)),
            ("plant", 101, "plant_cross_mates", (pair(101, 4000)[0], pair(9606, 4000)[1])),
            ("plant", 101, "plant_unmapped", (dna(150), dna(150))),
            ("bacteria", 201, "bacteria_human", pair(9606, 3500)),
            ("bacteria", 201, "bacteria_correct", pair(201, 3500)),
        ]
        self.write_sample()
        self.cfg = dict(metadata_dir=str(self.tax.parent), reference_raw=str(refs),
                        work_dir=str(self.base / "work"), tools_dir=str(self.base / "tools"),
                        storage_layout="separated",
                        stage_alignment=dict(representatives=str(self.representatives),
                                             minimap2=minimap_binary() or str(self.base / "not-installed"),
                                             max_reference_gib=.001, threads=2, offtarget_per_stratum=100))

    def tearDown(self):
        self.temporary.cleanup()

    def write_sample(self):
        sizes = {group: sum(case[0] == group for case in self.cases) for group in {case[0] for case in self.cases}}
        rows = []
        for index, (group, tid, name, _) in enumerate(self.cases):
            rows.append(dict(ordinal=100 + 7 * index, read_id=name, group=group, taxid=tid,
                             species_taxid=tid, status="assigned_species", bin=group + "/species_other",
                             source_group=group, stratum=group, stratum_pairs=1000, sampled_pairs=sizes[group]))
        table(self.qc / "sample.tsv", rows)
        for mate in (1, 2):
            with gzip.open(self.qc / f"sample.R{mate}.fastq.gz", "wt") as handle:
                for _, _, name, sequences in self.cases:
                    handle.write(f"@{name}/{mate} original comment\n{sequences[mate - 1]}\n+\n" + "I" * 150 + "\n")

    def evidence(self, qc=None):
        with ((qc or self.qc) / "alignment/sampled-evidence.tsv").open() as handle:
            return {row["read_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def native(self):
        if not minimap_binary():
            self.skipTest("Actual minimap2 binary unavailable")

    def test_actual_mapper_group_competition_and_offtarget(self):
        self.native()
        report = evaluate_alignment(self.sample, "group", self.cfg)
        self.assertEqual(report["status"], "complete", report)
        self.assertEqual(report["tested_pairs"], len(self.cases))
        self.assertFalse(report["accuracy_estimated"])
        result = self.evidence()
        expected = dict(human_correct="group_compatible", human_wrong_plant="conflicting_competitor",
                        human_primate_tie="ambiguous", plant_correct="group_compatible",
                        plant_wrong_human="conflicting_competitor", plant_relative_tie="ambiguous",
                        plant_other_species="group_compatible", plant_cross_mates="ambiguous",
                        plant_unmapped="unmapped", bacteria_human="potential_missed_human_compatible",
                        bacteria_correct="other_panel_compatible")
        self.assertEqual({name: row["category"] for name, row in result.items()}, expected)
        self.assertEqual(result["human_correct"]["ordinal"], "100")
        self.assertEqual(result["human_correct"]["r1_read_id"], "human_correct/1")
        self.assertEqual(result["human_correct"]["concordant_pair"], "True")
        self.assertEqual(report["summary_path"], "alignment/summary.json")
        self.assertTrue((self.qc / report["summary_path"]).is_file())
        self.assertFalse(list((self.qc / "alignment").glob(".panel.*")))

    def test_actual_mapper_species_claim_differs_from_group_claim(self):
        self.native()
        species_qc = self.sample / "group/species.qc"
        shutil.copytree(self.qc, species_qc)
        report = evaluate_alignment(self.sample, "split-species", self.cfg)
        self.assertEqual(report["status"], "complete", report)
        result = self.evidence(species_qc)
        self.assertEqual(result["plant_correct"]["category"], "species_compatible")
        self.assertEqual(result["plant_other_species"]["category"], "conflicting_competitor")
        self.assertEqual(result["plant_other_species"]["claim_level"], "species")

    def test_missing_mapper_is_explicit_and_accounts_for_selected_pairs(self):
        self.cfg["stage_alignment"]["minimap2"] = str(self.base / "missing-binary")
        report = evaluate_alignment(self.sample, "group", self.cfg)
        self.assertEqual(report["status"], "unavailable")
        self.assertEqual(report["tested_pairs"], 0)
        self.assertEqual(report["untested_pairs"], len(self.cases))
        self.assertEqual(report["errors"], [])
        self.assertIn("minimap2", report["reason"])

    def test_offtarget_sampling_weight_includes_second_sampling_step(self):
        self.cfg["stage_alignment"].update(offtarget_per_stratum=1, minimap2=str(self.base / "missing-binary"))
        report = evaluate_alignment(self.sample, "group", self.cfg)
        selected = [row for row in self.evidence().values() if row["role"] == "offtarget_human_screen"]
        self.assertEqual(len(selected), 1)
        self.assertAlmostEqual(float(selected[0]["evidence_weight"]), 1000)
        self.assertAlmostEqual(float(selected[0]["joint_inclusion_probability"]), .001)
        summary = next(row for row in report["by_stratum"] if row["stratum"] == "bacteria")
        self.assertEqual(summary["sampled_pairs"], 2)
        self.assertEqual(summary["selected_for_alignment"], 1)
        first = selected[0]["ordinal"]
        evaluate_alignment(self.sample, "group", self.cfg)
        second = [row for row in self.evidence().values() if row["role"] == "offtarget_human_screen"][0]["ordinal"]
        self.assertEqual(first, second)

    def test_missing_target_reference_reports_partial_without_dropping_claim(self):
        self.native()
        self.cases.append(("plant", 103, "missing_plant_reference", self.cases[3][3]))
        self.write_sample()
        report = evaluate_alignment(self.sample, "group", self.cfg)
        self.assertEqual(report["status"], "partial", report)
        missing = self.evidence()["missing_plant_reference"]
        self.assertEqual(missing["category"], "untested")
        self.assertIn("representative_missing", missing["untested_reason"])
        self.assertEqual(report["selected_pairs"], len(self.cases))

    def test_mapper_failure_is_failed_and_other_panels_are_attempted(self):
        mapper = self.base / "bad-minimap2"
        mapper.write_text('#!/bin/sh\nif [ "$1" = "--version" ]; then echo fixture; exit 0; fi\nexit 3\n')
        mapper.chmod(0o700)
        self.cfg["stage_alignment"]["minimap2"] = str(mapper)
        report = evaluate_alignment(self.sample, "group", self.cfg)
        self.assertEqual(report["status"], "failed")
        self.assertEqual(len(report["commands"]), 2)
        self.assertEqual(len(report["errors"]), 2)
        self.assertEqual(report["untested_pairs"], len(self.cases))

    def test_sample_header_or_count_mismatch_fails(self):
        with (self.qc / "sample.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows[0]["read_id"] = "wrong-order"
        table(self.qc / "sample.tsv", rows)
        report = evaluate_alignment(self.sample, "group", self.cfg)
        self.assertEqual(report["status"], "failed")
        self.assertIn("order/ID", report["reason"])
        saved = json.loads((self.qc / "alignment/summary.json").read_text())
        self.assertEqual(saved["status"], "failed")

    def test_group_coarse_taxid_without_species_is_evaluated_at_group_level(self):
        self.native()
        with (self.qc / "sample.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            if row["group"] == "plant":
                row["taxid"], row["species_taxid"] = "33090", "0"
        table(self.qc / "sample.tsv", rows)
        report = evaluate_alignment(self.sample, "group", self.cfg)
        self.assertEqual(report["status"], "complete", report)
        result = self.evidence()
        self.assertEqual(result["plant_correct"]["category"], "group_compatible")
        self.assertEqual(result["plant_other_species"]["category"], "group_compatible")
        self.assertEqual(result["plant_correct"]["claim_level"], "group")

    def test_compact_reference_catalog_selects_one_assembly_per_species(self):
        self.native()
        catalog = self.base / "library/group"
        catalog.mkdir(parents=True)
        with self.representatives.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            row.update(relative_path=Path(row["path"]).name, source_taxid=row["species_taxid"],
                       taxid=row["species_taxid"], records="1", compressed_bytes="9000", accession="GCF_test")
        # A second assembly for a species is legitimate in the compact catalog.
        duplicate = dict(rows[0])
        duplicate["accession"] = "GCA_second"
        rows.append(duplicate)
        table(catalog / "references.tsv", rows)
        self.cfg["stage_alignment"].pop("representatives")
        self.cfg["library_dir"] = str(catalog.parent)
        report = evaluate_alignment(self.sample, "group", self.cfg)
        self.assertEqual(report["status"], "complete", report)
        self.assertEqual(report["reference_table"], str(catalog / "references.tsv"))
        self.assertEqual(self.evidence()["human_correct"]["category"], "group_compatible")


if __name__ == "__main__":
    unittest.main()
