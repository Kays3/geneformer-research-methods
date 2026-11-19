import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scrna_isp.cli import analyze
from scrna_isp.encoders import score_record


class BenchmarkTests(unittest.TestCase):
    def test_relative_shift_can_hide_loss_of_goal_similarity(self):
        import math

        def vec(d):
            a = math.radians(d)
            return [math.cos(a), math.sin(a)]

        rows = [
            {
                "cell": s,
                "donor": "r",
                "role": "reference",
                "state": s,
                "original": vec(d),
            }
            for s, d in [("high", 0), ("mid", 45), ("low", 90)]
        ]
        rows.append(
            {
                "cell": "x",
                "donor": "x",
                "role": "evaluation",
                "state": "mid",
                "original": vec(65),
                "target": vec(116),
                "controls": [vec(70)],
                "control_ids": ["b"],
                "target_id": "a",
            }
        )
        _cov, out = analyze(rows, False)
        goal = next(r for r in out if r["goal"] == "low")
        self.assertLess(goal["direct"], 0)
        self.assertGreater(goal["relative"], 0)
        self.assertEqual(goal["aligned"], 1)
        rows[-1]["donor"] = "r"
        with self.assertRaises(ValueError):
            analyze(rows, False)
        self.assertEqual(len(analyze(rows, True)[0]), 1)

    def test_fixed_set_deletion_preserves_other_values(self):
        r = {
            "cell": "c",
            "donor": "d",
            "role": "evaluation",
            "state": "mid",
            "genes": ["<cls>", "a", "b", "c", "d", "e"],
            "tokens": [0, 1, 2, 3, 4, 5],
            "values": [0, 5, 6, 7, 8, 9],
            "target_id": "set",
            "control_ids": ["null"],
            "target_genes": ["a", "b"],
            "control_sets": [["d", "e"]],
        }
        calls = []

        def encode(ids, values):
            calls.append((ids, values))
            return [1, 2]

        got = score_record(r, "scgpt", encode)
        self.assertEqual(calls[1], ([0, 3, 4, 5], [0, 7, 8, 9]))
        self.assertEqual(calls[2], ([0, 1, 2, 3], [0, 5, 6, 7]))
        self.assertEqual(got["control_sets"], [["d", "e"]])

    def test_packaged_demo(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / "demo"
            subprocess.run(
                [sys.executable, "-m", "scrna_isp.cli", "demo", str(run)], check=True
            )
            receipt = json.loads((run / "receipt.json").read_text())
            self.assertEqual(receipt["eligible_cells"], 4)
            with (run / "donor_summary.tsv").open() as f:
                self.assertEqual(len(list(csv.DictReader(f, delimiter="\t"))), 8)


class CompareTests(unittest.TestCase):
    def test_intersection_preserves_native_denominators_and_rejects_mismatched_edits(
        self,
    ):
        from unittest.mock import patch

        from scrna_isp import compare

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [
                {
                    "cell": s,
                    "donor": "r",
                    "role": "reference",
                    "state": s,
                    "original": v,
                }
                for s, v in [("low", [0, 1]), ("mid", [1, 1]), ("high", [1, 0])]
            ]
            rows += [
                {
                    "cell": x,
                    "donor": "d",
                    "role": "evaluation",
                    "state": "high",
                    "original": [1, 0],
                    "target": [1, 1],
                    "controls": [[1, 0]],
                    "target_id": "a",
                    "control_ids": ["b"],
                }
                for x in ("x", "y")
            ]

            def run(name):
                inp = root / (name + ".jsonl")
                inp.write_text("".join(json.dumps(r) + "\n" for r in rows))
                subprocess.run(
                    [
                        sys.executable,
                        "-m",
                        "scrna_isp.cli",
                        "benchmark",
                        str(inp),
                        str(root / name),
                        "--model",
                        name,
                    ],
                    check=True,
                )

            run("one")
            del rows[-1]["target"], rows[-1]["controls"]
            run("two")
            args = [
                "compare",
                str(root / "matched.tsv"),
                str(root / "one"),
                str(root / "two"),
                "--intersection",
            ]
            with patch.object(sys, "argv", args):
                compare.main()
            with (root / "matched.tsv").open() as f:
                result = list(csv.DictReader(f, delimiter="\t"))
            self.assertTrue(all(r["cells"] == "1" for r in result))
            self.assertEqual(
                {r["native_eligible_start_cells"] for r in result}, {"1", "2"}
            )
            rows[-2]["control_ids"] = ["other"]
            run("three")
            args[1] = str(root / "invalid.tsv")
            args[3] = str(root / "three")
            with patch.object(sys, "argv", args), self.assertRaises(SystemExit):
                compare.main()
