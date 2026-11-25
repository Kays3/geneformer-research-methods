"""Portable count baseline and matched, three-state embedding benchmark."""

import argparse
import csv
import hashlib
import json
import math
import platform
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path
from statistics import mean

STATES = ("high", "mid", "low")
META = ("cell", "donor", "role", "state")


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def unit(v):
    norm = math.hypot(*v)
    if not math.isfinite(norm) or norm == 0:
        raise ValueError("Cosine requires a finite, nonzero embedding")
    return [x / norm for x in v]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def write_table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_counts(path, target, control):
    if target == control:
        raise ValueError("Choose distinct target and control genes")
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        names = reader.fieldnames or []
        if len(set(names)) != len(names) or not {*META, target, control} <= set(names):
            raise ValueError(
                "Expected unique gene headers, metadata, target and control"
            )
        genes = [g for g in names if g not in META]
        if target not in genes or control not in genes:
            raise ValueError("Target and control must be gene columns")
        rows = []
        for r in reader:
            if None in r:
                raise ValueError("Unexpected extra count columns")
            counts = [int(r[g]) for g in genes]
            if min(counts) < 0 or sum(counts) == 0:
                raise ValueError(
                    "Counts must be nonnegative with positive library size"
                )
            library_size = sum(counts)
            emb = [math.log1p(10000 * x / library_size) for x in counts]
            out = {k: r[k] for k in META}
            out.update(original=emb, target_id=target, control_ids=[control])
            if r["role"] == "evaluation" and all(
                counts[genes.index(g)] > 0 for g in (target, control)
            ):
                out["target"] = emb.copy()
                out["target"][genes.index(target)] = 0
                null = emb.copy()
                null[genes.index(control)] = 0
                out["controls"] = [null]
            rows.append(out)
    return rows


def nearest(values):
    winners = [s for s in STATES if values[s] == max(values.values())]
    return winners[0] if len(winners) == 1 else "tie"


def analyze(rows, descriptive):
    seen, donors, refs = set(), defaultdict(set), defaultdict(list)
    dim = None
    target_id = None
    for r in rows:
        if any(not isinstance(r.get(k), str) or not r[k].strip() for k in META):
            raise ValueError("Each row needs cell, donor or pool ID, role and state")
        if (
            r["role"] not in {"reference", "evaluation"}
            or r["state"] not in STATES
            or r["cell"] in seen
        ):
            raise ValueError("Duplicate cell or unrecognized state or role")
        seen.add(r["cell"])
        donors[r["role"]].add(r["donor"])
        original = r["original"]
        unit(original)
        if dim is None:
            dim = len(original)
        if len(original) != dim:
            raise ValueError("Inconsistent embedding dimensions")
        if r["role"] == "reference":
            refs[r["state"]].append(original)
        else:
            if not isinstance(r.get("target_id"), str) or not r["target_id"].strip():
                raise ValueError("Evaluation rows require a target identifier")
            target_id = target_id or r["target_id"]
            if target_id != r["target_id"]:
                raise ValueError("Use one registered intervention per benchmark")
    if donors["reference"] & donors["evaluation"] and not descriptive:
        raise ValueError(
            "Overlapping donors: supply independent donors or use --within-library"
        )
    if any(not refs[s] for s in STATES) or not donors["evaluation"]:
        raise ValueError("All three state references and evaluation cells are required")
    centers = {s: unit([mean(col) for col in zip(*refs[s])]) for s in STATES}

    def scores(vector):
        if len(vector) != dim:
            raise ValueError("Edited embedding dimension differs from reference")
        v = unit(vector)
        return {s: dot(v, c) for s, c in centers.items()}

    coverage, shifts = [], []
    for r in rows:
        if r["role"] != "evaluation":
            continue
        orig = scores(r["original"])
        start = r["state"]
        eligible = "target" in r
        if eligible != ("controls" in r):
            raise ValueError("Supply target and control embeddings together")
        coverage.append(
            {
                "cell": r["cell"],
                "donor": r["donor"],
                "state": start,
                "nearest": nearest(orig),
                "eligible": int(eligible),
            }
        )
        if not eligible:
            continue
        ids = r.get("control_ids", [])
        if (
            not isinstance(ids, list)
            or not ids
            or not all(isinstance(i, str) and i.strip() for i in ids)
            or len(set(ids)) != len(ids)
            or r["target_id"] in ids
            or len(ids) != len(r["controls"])
        ):
            raise ValueError(
                "Supply distinct matched control IDs and one vector per control"
            )
        target = scores(r["target"])
        controls = [scores(v) for v in r["controls"]]
        for goal in STATES:
            if start == goal:
                continue
            alternate = next(s for s in STATES if s not in (start, goal))
            adj = {s: target[s] - mean(c[s] for c in controls) for s in STATES}
            aligned = nearest(orig) == start
            members = r.get("target_genes", [r["target_id"]])
            sets = r.get("control_sets", [[c] for c in ids])
            if (
                len(sets) != len(ids)
                or not members
                or len(set(members)) != len(members)
                or any(
                    len(set(g)) != len(g)
                    or len(g) != len(members)
                    or set(g) & set(members)
                    for g in sets
                )
            ):
                raise ValueError("Invalid paired gene-set identities")
            shifts.append(
                {
                    "target_genes": json.dumps(sorted(members)),
                    "control_sets": json.dumps(sorted(sorted(g) for g in sets)),
                    "cell": r["cell"],
                    "donor": r["donor"],
                    "start": start,
                    "goal": goal,
                    "alternate": alternate,
                    "target_id": r["target_id"],
                    "control_ids": json.dumps(ids),
                    "original_nearest": nearest(orig),
                    "target_nearest": nearest(target),
                    "raw_start": target[start] - orig[start],
                    "raw_alternate": target[alternate] - orig[alternate],
                    "adjusted_start": adj[start],
                    "raw_goal": target[goal] - orig[goal],
                    "control_goal": mean(c[goal] for c in controls) - orig[goal],
                    "direct": adj[goal],
                    "relative": adj[goal] - adj[start],
                    "alternate_shift": adj[alternate],
                    "aligned": int(aligned),
                    "target_cross": int(aligned and nearest(target) == goal),
                    "control_cross": mean(
                        int(aligned and nearest(c) == goal) for c in controls
                    ),
                }
            )
    return coverage, shifts


def run(args, folder):
    if args.command == "demo":
        source = Path(__file__).with_name("example") / "counts.tsv"
        shutil.copyfile(source, folder / "counts.tsv")
        args.target, args.control = "GENE_A", "GENE_B"
    else:
        source = args.input
    if args.command == "benchmark":
        with source.open() as stream:
            rows = [json.loads(line) for line in stream if line.strip()]
    else:
        rows = load_counts(source, args.target, args.control)
    coverage, shifts = analyze(rows, args.within_library)
    write_table(folder / "coverage.tsv", coverage)
    if shifts:
        write_table(folder / "shifts.tsv", shifts)
    groups = defaultdict(list)
    for r in shifts:
        groups[r["donor"], r["start"], r["goal"]].append(r)
    summaries = [
        {
            "donor": d,
            "start": s,
            "goal": g,
            "n": len(rs),
            "direct": mean(r["direct"] for r in rs),
            "relative": mean(r["relative"] for r in rs),
            "aligned": sum(r["aligned"] for r in rs),
            "target_cross": sum(r["target_cross"] for r in rs),
            "control_cross": sum(r["control_cross"] for r in rs),
        }
        for (d, s, g), rs in sorted(groups.items())
    ]
    if summaries:
        write_table(folder / "donor_summary.tsv", summaries)
    recall = {
        s: [r["nearest"] == s for r in coverage if r["state"] == s] for s in STATES
    }
    receipt = {
        "model": args.model,
        "scope": "within-library descriptive"
        if args.within_library
        else "donor-disjoint",
        "reference_manifest_sha256": hashlib.sha256(
            json.dumps(
                sorted(
                    (r["cell"], r["donor"], r["state"])
                    for r in rows
                    if r["role"] == "reference"
                )
            ).encode()
        ).hexdigest(),
        "input_sha256": sha(source),
        "code_sha256": sha(Path(__file__)),
        "python": platform.python_version(),
        "target": next(r["target_id"] for r in rows if r["role"] == "evaluation"),
        "reference_cells": len(rows) - len(coverage),
        "evaluation_cells": len(coverage),
        "eligible_cells": sum(r["eligible"] for r in coverage),
        "per_state_recall": {s: mean(v) if v else None for s, v in recall.items()},
        "outputs": {p.name: sha(p) for p in sorted(folder.iterdir())},
    }
    (folder / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("demo", "baseline", "benchmark"):
        p = sub.add_parser(name)
        if name != "demo":
            p.add_argument("input", type=Path)
        p.add_argument("output", type=Path)
        p.add_argument(
            "--within-library",
            action="store_true",
            help="Label a shared-pool comparison descriptive",
        )
        if name == "benchmark":
            p.add_argument(
                "--model", required=True, help="Checkpoint or adapter identifier"
            )
        else:
            p.set_defaults(model="log-normalized counts, fixed pre-edit library size")
        if name == "baseline":
            p.add_argument("--target", required=True)
            p.add_argument("--control", required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new output directory")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.output.parent) as temp:
        stage = Path(temp) / "result"
        stage.mkdir()
        try:
            run(args, stage)
            stage.rename(args.output)
        except (OSError, ValueError, TypeError, KeyError) as e:
            parser.exit(2, f"scrna-isp: {e}\n")


if __name__ == "__main__":
    main()
