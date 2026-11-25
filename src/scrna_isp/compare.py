"""Compare state shifts only when model exports cover identical cells and edits."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from statistics import mean

from .cli import sha, write_table


def paired(folder):
    receipt = json.loads((folder / "receipt.json").read_text())
    for name in ["coverage.tsv", "shifts.tsv"]:
        if (
            name not in receipt["outputs"]
            or sha(folder / name) != receipt["outputs"][name]
        ):
            raise ValueError(f"Missing, empty or modified output: {folder / name}")
    with (folder / "coverage.tsv").open() as f:
        coverage = sorted(
            (r["cell"], r["donor"], r["state"], r["eligible"])
            for r in csv.DictReader(f, delimiter="\t")
        )
    with (folder / "shifts.tsv").open() as f:
        records = {}
        for r in csv.DictReader(f, delimiter="\t"):
            key = (
                r["cell"],
                r["donor"],
                r["start"],
                r["goal"],
                r["target_id"],
                tuple(sorted(json.loads(r["control_ids"]))),
                r["target_genes"],
                r["control_sets"],
            )
            if key in records:
                raise ValueError("Duplicate cell and transition")
            records[key] = r
    return receipt, coverage, records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument("runs", type=Path, nargs="+")
    p.add_argument(
        "--intersection",
        action="store_true",
        help="Restrict to naturally eligible cells in every model; retain native counts",
    )
    a = p.parse_args()
    if a.output.exists() or len(a.runs) < 2:
        p.error("Supply a new output and at least two model runs")
    results = []
    runs = [paired(path) for path in a.runs]
    names = [receipt["model"] for receipt, _, _ in runs]
    if len(set(names)) != len(names):
        p.error("Give each checkpoint a distinct model identifier")
    base, cohort, first = runs[0]
    common = {r[0] for r in cohort if r[3] == "1"}
    for receipt, coverage, _ in runs[1:]:
        if [r[:3] for r in coverage] != [r[:3] for r in cohort] or any(
            receipt[k] != base[k]
            for k in ("scope", "target", "reference_manifest_sha256")
        ):
            p.error(
                "Evaluation identities, reference identities, state labels, target or scope differ"
            )
        if not a.intersection and coverage != cohort:
            p.error("Eligibility differs. Register or select the shared intersection.")
        common.intersection_update(r[0] for r in coverage if r[3] == "1")
    if not common:
        p.error("No commonly eligible evaluation cells")
    keys = {k for k in first if k[0] in common}
    for receipt, coverage, all_records in runs:
        model = receipt["model"]
        records = {k: v for k, v in all_records.items() if k[0] in common}
        if records.keys() != keys:
            p.error("Paired targets or control identities differ on shared cells")
        groups = defaultdict(lambda: defaultdict(list))
        for row in records.values():
            groups[row["start"], row["goal"]][row["donor"]].append(row)
        for (start, goal), donors in sorted(groups.items()):
            direct = [
                mean(float(r["direct"]) for r in cells) for cells in donors.values()
            ]
            relative = [
                mean(float(r["relative"]) for r in cells) for cells in donors.values()
            ]
            results.append(
                {
                    "model": model,
                    "scope": receipt["scope"],
                    "native_eligible_start_cells": sum(
                        r[2] == start and r[3] == "1" for r in coverage
                    ),
                    "registered_start_cells": sum(r[2] == start for r in coverage),
                    "start": start,
                    "goal": goal,
                    "units": len(donors),
                    "cells": sum(map(len, donors.values())),
                    "direct": mean(direct),
                    "relative": mean(relative),
                    "units_relative_positive": sum(x > 0 for x in relative),
                }
            )
    write_table(a.output, results)
