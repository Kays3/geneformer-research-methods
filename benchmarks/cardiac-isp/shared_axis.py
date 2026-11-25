"""Descriptive agreement on at least one named axis across all three models."""
import argparse
import csv
import json
import math

from balance import MODELS


def shared_axis(rows, expected_donors):
    table = {}
    for row in rows:
        donor, model = row["donor"], row["model"]
        values = float(row["nf"]), float(row["hcm"])
        if not donor or model not in MODELS or not all(map(math.isfinite, values)):
            raise ValueError("Invalid identity or nonfinite effect")
        scores = table.setdefault(donor, {})
        if model in scores:
            raise ValueError("Duplicate donor and model")
        scores[model] = values
    if len(table) != expected_donors or any(set(s) != set(MODELS) for s in table.values()):
        raise ValueError("Expected the full matched donor panel for every model")
    nf, hcm, both, either, neither = [], [], [], [], []
    for donor, scores in sorted(table.items()):
        nf_pass = all(x > 0 for x, y in scores.values())
        hcm_pass = all(y < 0 for x, y in scores.values())
        if nf_pass:
            nf.append(donor)
        if hcm_pass:
            hcm.append(donor)
        if nf_pass and hcm_pass:
            both.append(donor)
        (either if nf_pass or hcm_pass else neither).append(donor)
    return dict(donors=len(table), all_nf=nf, all_hcm=hcm, both=both, either=either, neither=neither)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("donor_tsv")
    parser.add_argument("--expected-donors", type=int, default=9)
    args = parser.parse_args()
    with open(args.donor_tsv) as stream:
        print(json.dumps(shared_axis(csv.DictReader(stream, delimiter="\t"), args.expected_donors), indent=2))
