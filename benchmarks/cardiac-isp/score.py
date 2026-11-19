"""Score one predeclared cardiac benchmark arm. Inputs and outputs stay local."""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean


def rows(path, columns):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not columns.issubset(reader.fieldnames or []):
            raise ValueError(f"{path} requires columns {sorted(columns)}")
        for row in reader:
            if any(not (row.get(k) or "").strip() for k in columns):
                raise ValueError(f"Missing value in {path} at row {reader.line_num}")
            yield row


def number(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Nonfinite input")
    return value


def cell_scores(cohort, scores, margin=0.0):
    """One partition, one checkpoint, one readout, one prespecified control draw."""
    if not math.isfinite(margin) or margin < 0:
        raise ValueError("The registered margin must be finite and nonnegative")
    owners = {}
    for row in rows(cohort, {"cell", "donor"}):
        if row["cell"] in owners:
            raise ValueError("Duplicate cohort cell")
        owners[row["cell"]] = row["donor"]
    seen = set()
    groups = defaultdict(list)
    fields = {"cell", "donor", "target_nf", "target_hcm", "control_nf", "control_hcm"}
    for r in rows(scores, fields):
        cell = r["cell"]
        if owners.get(cell) != r["donor"] or cell in seen:
            raise ValueError("Cell identity, donor, or uniqueness check failed")
        seen.add(cell)
        nf = number(r["target_nf"]) - number(r["control_nf"])
        hcm = number(r["target_hcm"]) - number(r["control_hcm"])
        groups[r["donor"]].append((nf, hcm))
    if not owners or seen != set(owners):
        raise ValueError("Scores must cover exactly the registered analysis cohort")
    donors = []
    for donor, points in sorted(groups.items()):
        nf, hcm = map(mean, zip(*points))
        donors.append(dict(donor=donor, cells=len(points), nf=nf, hcm=hcm,
                           joint=(nf > margin and hcm < -margin)))
    return dict(cells=len(seen), donors=len(donors), joint=sum(r["joint"] for r in donors),
                equal_donor_nf=mean(r["nf"] for r in donors),
                equal_donor_hcm=mean(r["hcm"] for r in donors), by_donor=donors)


def cis_scores(path):
    """One predeclared primary allele/track per locus and comparator cohort."""
    by_locus = []
    seen = set()
    for r in rows(path, {"locus", "variant", "prediction", "beta", "se"}):
        if r["locus"] in seen:
            raise ValueError("Multiple primary observations for one locus")
        seen.add(r["locus"])
        pred = number(r["prediction"])
        row = dict(locus=r["locus"], variant=r["variant"], prediction=pred, agreement=None)
        if (r["beta"] == "NA") != (r["se"] == "NA"):
            raise ValueError("Missing measurements require NA in both beta and se")
        if r["beta"] != "NA":
            beta, se = number(r["beta"]), number(r["se"])
            if se <= 0:
                raise ValueError("Expected a positive standard error")
            lo, hi = beta - 1.96 * se, beta + 1.96 * se
            row.update(beta=beta, ci95=[lo, hi], agreement=pred * beta > 0,
                       directional=pred != 0 and beta != 0)
        by_locus.append(row)
    if not by_locus:
        raise ValueError("Empty locus panel")
    informative = [r for r in by_locus if r.get("directional")]
    return dict(registered_loci=len(seen), measured_loci=sum("beta" in r for r in by_locus),
                directional_loci=len(informative), same_sign=sum(r["agreement"] for r in informative),
                by_locus=by_locus)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="mode", required=True)
    cell = subs.add_parser("cell")
    cell.add_argument("cohort", type=Path)
    cell.add_argument("scores", type=Path)
    cell.add_argument("--margin", type=float, default=0.0)
    cis = subs.add_parser("cis")
    cis.add_argument("alleles", type=Path)
    args = parser.parse_args()
    result = cell_scores(args.cohort, args.scores, args.margin) if args.mode == "cell" else cis_scores(args.alleles)
    print(json.dumps(result, indent=2, allow_nan=False))
