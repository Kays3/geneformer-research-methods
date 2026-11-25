"""Register a matched single-cell ISP comparison and summarize measured model shifts."""
import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

HERE = Path(__file__).resolve().parent
MODELS = {"v1", "v2", "scgpt"}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path, fields):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not fields <= set(reader.fieldnames or ()):
            raise ValueError(f"{path.name} needs {sorted(fields)}")
        for row in reader:
            if any(not (row.get(k) or "").strip() for k in fields):
                raise ValueError(f"{path.name}:{reader.line_num}: empty required field")
            yield row


def inspect(folder):
    cfg = json.loads((folder / "study.json").read_text())
    required = ["study", "dataset", "species", "cell_type", "count_layer", "target", "edit",
                "control_rule", "positive_control", "validation_plan", "null_and_multiplicity",
                "checkpoint_selection", "eligibility", "state_reference", "software_environment"]
    for key in required:
        if not isinstance(cfg[key], str) or not cfg[key].strip() or cfg[key].startswith("FILL"):
            raise ValueError(f"Complete study.json field {key}")
    models, axes = cfg["models"], cfg["axes"]
    if not isinstance(models, dict) or not models or not set(models) <= MODELS:
        raise ValueError("Select supported transcriptome models")
    if not axes or any(type(sign) is not int or sign not in (-1, 1) for sign in axes.values()):
        raise ValueError("Register at least one axis with direction +1 or -1")
    for model, spec in models.items():
        for field in ["checkpoint", "tokenizer", "readout"]:
            if not spec[field].strip() or spec[field].startswith("FILL"):
                raise ValueError(f"Complete {model}.{field}")
        if set(spec["margins"]) != set(axes) or any(not math.isfinite(v) or v < 0 for v in spec["margins"].values()):
            raise ValueError(f"Register a finite nonnegative margin per axis for {model}")
    splits = {}
    for row in read_table(folder / "donors.tsv", {"donor", "partition"}):
        d, part = row["donor"], row["partition"]
        if d in splits or part not in {"train", "validation", "test"}:
            raise ValueError("Duplicate donor or unknown partition")
        splits[d] = part
    cells = {}
    for row in read_table(folder / "cells.tsv", {"cell", "donor", "controls"}):
        cell, donor = row["cell"], row["donor"]
        controls = row["controls"].split(",")
        if cell in cells or donor not in splits or not all(controls) or len(set(controls)) != len(controls) or cfg["target"] in controls:
            raise ValueError("Invalid paired cell, donor or control IDs")
        cells[cell] = (donor, frozenset(controls))
    coverage = {d for d, _ in cells.values()}
    panel = {d for d, p in splits.items() if p == "test"}
    if not panel or not panel <= coverage or not cells or any(splits[d] != "test" for d in coverage):
        raise ValueError("cells.tsv must cover all and only the registered test donors")
    if type(cfg["min_shared_donors"]) is not int or not 1 <= cfg["min_shared_donors"] <= len(panel):
        raise ValueError("Register an attainable shared-donor count")
    return cfg, cells


def fingerprint(folder):
    return {name: digest(folder / name) for name in ["study.json", "cells.tsv", "donors.tsv"]}


def score_direction(folder, shifts):
    cfg, cells = inspect(folder)
    if json.loads((folder / "lock.json").read_text()) != fingerprint(folder):
        raise ValueError("Registered inputs have changed")
    models, axes = cfg["models"], cfg["axes"]
    data = defaultdict(list)
    seen = set()
    for r in read_table(shifts, {"cell", "model", "axis", "target_shift", "control_mean_shift", "controls"}):
        cell, model, axis = r["cell"], r["model"], r["axis"]
        key = cell, model, axis
        if key in seen or cell not in cells or model not in models or axis not in axes:
            raise ValueError("Unexpected or duplicate score")
        donor, controls = cells[cell]
        ids = r["controls"].split(",")
        if set(ids) != controls or len(ids) != len(controls):
            raise ValueError("Scored genes differ from registered cell controls")
        effect = float(r["target_shift"]) - float(r["control_mean_shift"])
        if not math.isfinite(effect):
            raise ValueError("Nonfinite shift")
        seen.add(key)
        data[model, donor, axis].append(effect)
    if len(seen) != len(cells) * len(models) * len(axes):
        raise ValueError("Incomplete cell/model/axis coverage")
    donors = sorted({d for d, _ in cells.values()})
    result = {}
    shared = set(donors)
    for model, spec in models.items():
        by_donor = {}
        passing = set()
        for donor in donors:
            effects = {axis: mean(data[model, donor, axis]) for axis in axes}
            passed = all(axes[a] * v > spec["margins"][a] for a, v in effects.items())
            by_donor[donor] = dict(cells=len(data[model, donor, next(iter(axes))]), effects=effects, passed=passed)
            if passed:
                passing.add(donor)
        shared &= passing
        result[model] = dict(passing=len(passing), by_donor=by_donor,
                             equal_donor_effect={a: mean(by_donor[d]["effects"][a] for d in donors) for a in axes})
    return dict(study=cfg["study"], test_cells=len(cells), test_donors=len(donors), models=result,
                passing_in_all_models=sorted(shared), directional_target_met=len(shared) >= cfg["min_shared_donors"],
                registered_inputs=fingerprint(folder), shifts_sha256=digest(shifts))


def score(folder, shifts, overexpression=None):
    deletion = score_direction(folder, shifts)
    cfg, cells = inspect(folder)
    if overexpression is None:
        if cfg.get("paired_directions") is True:
            raise ValueError("This registration requires an overexpression export")
        return deletion
    if shifts.resolve() == overexpression.resolve():
        raise ValueError("Supply separate deletion and overexpression exports")
    increase = score_direction(folder, overexpression)
    results = {}
    for model in cfg["models"]:
        donors = {}
        for donor, row in deletion["models"][model]["by_donor"].items():
            down = row["effects"]
            up = increase["models"][model]["by_donor"][donor]["effects"]
            donors[donor] = {
                "cells": row["cells"], "deletion": down, "overexpression": up,
                "deletion_minus_overexpression": {axis: down[axis]-up[axis] for axis in down},
            }
        results[model] = {
            "by_donor": donors,
            "equal_donor_effects": {
                direction: {axis: mean(row[direction][axis] for row in donors.values())
                            for axis in cfg["axes"]}
                for direction in ("deletion", "overexpression", "deletion_minus_overexpression")
            },
        }
    return {"study": cfg["study"], "test_cells": len(cells), "models": results,
            "registered_inputs": fingerprint(folder), "deletion_sha256": digest(shifts),
            "overexpression_sha256": digest(overexpression)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ["init", "check", "lock", "score"]:
        sub = commands.add_parser(name)
        sub.add_argument("folder", type=Path)
        if name == "score":
            sub.add_argument("shifts", type=Path)
            sub.add_argument("--overexpression", type=Path, help="Matched overexpression shift export")
    args = parser.parse_args()
    folder = args.folder
    if args.command == "init":
        folder.mkdir(parents=True, exist_ok=False)
        for name in ("study.json", "cells.tsv", "donors.tsv"):
            shutil.copyfile(HERE / "template" / name, folder / name)
        print("Created study inputs at", folder)
    elif args.command == "score":
        print(json.dumps(score(folder, args.shifts, args.overexpression), indent=2, allow_nan=False))
    else:
        cfg, cells = inspect(folder)
        if args.command == "lock":
            with (folder / "lock.json").open("x") as f:
                json.dump(fingerprint(folder), f, indent=2)
        print(f"{cfg['study']}: {len(cells)} paired test cells; {len(cfg['models'])} model(s)")


if __name__ == "__main__":
    try:
        main()
    except (OSError, KeyError, ValueError, TypeError) as e:
        sys.exit(str(e))
