"""Request one human cis-expression score per assay track from AlphaGenome 2."""

import argparse
import json
import os
import tempfile
from importlib.metadata import version
from pathlib import Path

from .cli import sha


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path, help="New output directory")
    for name in ("chromosome", "ref", "alt", "gene", "condition"):
        p.add_argument("--" + name, required=True)
    p.add_argument(
        "--position", required=True, type=int, help="1-based GRCh38 position"
    )
    p.add_argument("--fold", default="ALL_FOLDS")
    a = p.parse_args()
    if (
        a.position < 1
        or a.ref == a.alt
        or not all(s and set(s) <= set("ACGT") for s in [a.ref, a.alt])
    ):
        p.error("Expected a positive position and distinct DNA alleles")
    if a.output.exists():
        p.error("Choose a new output directory")
    from alphagenome.data import genome
    from alphagenome.models.v2 import dna_client, variant_scorers

    client = dna_client.create(
        os.environ["ALPHAGENOME_API_KEY"],
        model_version=dna_client.ModelVersion[a.fold],
        timeout=90,
    )
    variant = genome.Variant(a.chromosome, a.position, a.ref, a.alt)
    interval = variant.reference_interval.resize(1048576)
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION"]
    tracks = variant_scorers.tidy_scores(
        client.score_variant(
            interval=interval,
            variant=variant,
            condition=a.condition,
            variant_scorers=[scorer],
        )
    )
    selected = tracks.loc[tracks.gene_id == a.gene]
    if selected.empty or not selected.track_name.is_unique:
        raise ValueError("Missing or ambiguous target-gene tracks")
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.output.parent) as tmp:
        stage = Path(tmp) / "run"
        stage.mkdir()
        path = stage / "tracks.tsv"
        selected.to_csv(path, sep="\t", index=False)
        (stage / "receipt.json").write_text(
            json.dumps(
                {
                    "build": "GRCh38",
                    "variant": str(variant),
                    "interval": str(interval),
                    "gene": a.gene,
                    "condition": a.condition,
                    "fold": a.fold,
                    "scorer": "SC_GENE_EXPRESSION",
                    "sdk": version("alphagenome"),
                    "tracks_sha256": sha(path),
                    "code_sha256": sha(Path(__file__)),
                },
                indent=2,
            )
            + "\n"
        )
        stage.rename(a.output)
