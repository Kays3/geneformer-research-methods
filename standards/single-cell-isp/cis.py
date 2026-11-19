"""Run the accompanying locus-level benchmark for AlphaGenome 2 output."""
import argparse
import json
from pathlib import Path
import shutil
import sys

from isp import HERE, digest, read_table

# Reuse the repository's established locus-level scoring rules.
sys.path.insert(0, str(HERE.parents[1] / "benchmarks/cardiac-isp"))
from score import cis_scores

FIELDS = {"locus", "variant", "gene", "ref", "alt"}
PARAMS = {"study", "build", "model", "condition", "track", "scorer", "window_rule",
          "comparison_study", "allele_harmonization", "selection_and_controls"}


def plan(folder):
    config = json.loads((folder / "cis.json").read_text())
    for k in PARAMS:
        if not isinstance(config.get(k), str) or not config[k].strip() or config[k].startswith("FILL"):
            raise ValueError(f"Complete cis.json field {k}")
    loci = {}
    for r in read_table(folder / "loci.tsv", FIELDS):
        if r["locus"] in loci or r["ref"] == r["alt"]:
            raise ValueError("Duplicate primary locus or identical REF and ALT")
        loci[r["locus"]] = tuple(r[k] for k in ("variant", "gene", "ref", "alt"))
    if not loci:
        raise ValueError("Empty locus ledger")
    return loci


def hashes(folder):
    return {k: digest(folder / k) for k in ["cis.json", "loci.tsv"]}


def evaluate(folder, path):
    loci = plan(folder)
    if json.loads((folder / "cis-lock.json").read_text()) != hashes(folder):
        raise ValueError("Registered cis inputs changed")
    observed = {}
    for r in read_table(path, FIELDS | {"prediction", "beta", "se"}):
        observed[r["locus"]] = tuple(r[k] for k in ("variant", "gene", "ref", "alt"))
    if observed != loci:
        raise ValueError("Scored alleles, genes or coverage differ from registered loci")
    result = cis_scores(path)
    result.update(registered_inputs=hashes(folder), scores_sha256=digest(path))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="action", required=True)
    for action in ("init", "check", "lock", "score"):
        sub = subs.add_parser(action)
        sub.add_argument("folder", type=Path)
        if action == "score":
            sub.add_argument("scores", type=Path)
    args = parser.parse_args()
    try:
        if args.action == "init":
            args.folder.mkdir(parents=True, exist_ok=False)
            for name in ("cis.json", "loci.tsv"):
                shutil.copyfile(HERE / "template" / name, args.folder / name)
        elif args.action == "score":
            print(json.dumps(evaluate(args.folder, args.scores), indent=2, allow_nan=False))
        else:
            loci = plan(args.folder)
            if args.action == "lock":
                with (args.folder / "cis-lock.json").open("x") as f:
                    json.dump(hashes(args.folder), f, indent=2)
            print(f"{len(loci)} registered primary loci")
    except (OSError, ValueError, KeyError, TypeError) as e:
        sys.exit(str(e))
