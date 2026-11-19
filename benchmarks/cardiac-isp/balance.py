"""Evaluate the nine-donor transfer target on aligned, control-adjusted donor means."""
import argparse
import csv
import json
import math
from pathlib import Path

MODELS = ("v1", "v2", "scgpt")


def balance(rows):
    scores = {model: {} for model in MODELS}
    for row in rows:
        model, donor = row["model"], row["donor"]
        if model not in scores or not donor or donor in scores[model]:
            raise ValueError("Expected one row per declared model and donor")
        nf, hcm = float(row["nf"]), float(row["hcm"])
        if not all(map(math.isfinite, (nf, hcm))):
            raise ValueError("Missing or nonfinite effect")
        scores[model][donor] = nf > 0 and hcm < 0
    donors = set(scores[MODELS[0]])
    if len(donors) != 9 or any(set(table) != donors for table in scores.values()):
        raise ValueError("All three checkpoints must have the same nine registered donors")
    counts = {model: sum(table.values()) for model, table in scores.items()}
    shared = sorted(d for d in donors if all(scores[m][d] for m in MODELS))
    return dict(per_model=counts, weakest_model_count=min(counts.values()),
                count_spread=max(counts.values()) - min(counts.values()),
                jointly_passing_donors=shared, joint_count=len(shared),
                target_met=len(shared) >= 8)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("donor_tsv", type=Path)
    args = parser.parse_args()
    with args.donor_tsv.open() as handle:
        print(json.dumps(balance(csv.DictReader(handle, delimiter="\t")), indent=2))
