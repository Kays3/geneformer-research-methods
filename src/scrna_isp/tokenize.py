"""Create model-native token records from a pre-mapped single-cell count layer."""

import argparse
import json
import pickle
import tempfile
from pathlib import Path

from .cli import META, sha


def tokenize(model, counts, vocab, medians, max_context):
    import numpy as np

    if not counts:
        raise ValueError("No mapped expressed features in this cell")
    if model == "scgpt":
        genes = sorted(counts, key=lambda g: (-counts[g], vocab[g]))
        raw = [counts[g] for g in genes]
        bins = np.digitize(raw, np.quantile(raw, np.linspace(0, 1, 50)))
        genes = genes[: max_context - 1]
        return {
            "genes": ["<cls>", *genes],
            "tokens": [vocab["<cls>"], *[vocab[g] for g in genes]],
            "values": [0, *map(int, bins[: len(genes)])],
        }
    special = model == "v2"
    genes = sorted(counts, key=lambda g: -counts[g] / medians[g])[
        : max_context - 2 * special
    ]
    if special:
        genes = ["<cls>", *genes, "<eos>"]
    return {"genes": genes, "tokens": [int(vocab[g]) for g in genes]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("model", choices=["v1", "v2", "scgpt"])
    p.add_argument("h5ad", type=Path)
    p.add_argument(
        "manifest",
        type=Path,
        help="TSV with cell, donor, role, state, target_id, control_ids (comma-separated)",
    )
    p.add_argument(
        "vocab",
        type=Path,
        help="Trusted native token dictionary pickle or scGPT vocab.json",
    )
    p.add_argument("output", type=Path)
    p.add_argument("--layer", required=True, help="Raw counts layer name, or X")
    p.add_argument(
        "--gene-column",
        required=True,
        help="Unique mapped IDs in var; use index for var_names",
    )
    p.add_argument(
        "--match-column",
        help="Unique common gene IDs in var for cross-model edit matching (default: gene-column)",
    )
    p.add_argument(
        "--medians",
        type=Path,
        help="Geneformer checkpoint-matched positive median dictionary",
    )
    p.add_argument(
        "--context", type=int, default=1200, help="scGPT context including CLS"
    )
    args = p.parse_args()
    import csv

    import anndata as ad
    import numpy as np
    from scipy import sparse

    if args.output.exists():
        p.error("Choose a new output directory")
    if args.model == "scgpt":
        if args.context < 2:
            p.error("scGPT context must accommodate CLS and genes")
        vocab = json.loads(args.vocab.read_text())
        medians = {}
        limit = args.context
    else:
        if not args.medians:
            p.error("Geneformer requires --medians")
        vocab = pickle.loads(args.vocab.read_bytes())
        medians = pickle.loads(args.medians.read_bytes())
        limit = 2048 if args.model == "v1" else 4096
    a = ad.read_h5ad(args.h5ad)
    matrix = a.X if args.layer == "X" else a.layers[args.layer]
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    identifiers = (
        a.var_names if args.gene_column == "index" else a.var[args.gene_column]
    )
    mapped = {
        i: g
        for i, g in enumerate(identifiers)
        if isinstance(g, str)
        and g in vocab
        and not g.startswith("<")
        and (
            args.model == "scgpt"
            or (g in medians and np.isfinite(medians[g]) and medians[g] > 0)
        )
    }
    shared = (
        identifiers
        if not args.match_column
        else (a.var_names if args.match_column == "index" else a.var[args.match_column])
    )
    canonical = {
        g: shared.iloc[i] if hasattr(shared, "iloc") else shared[i]
        for i, g in mapped.items()
    }
    if not all(
        isinstance(g, str) and g and not g.startswith("<") for g in canonical.values()
    ) or len(set(canonical.values())) != len(canonical):
        p.error("Supply a unique common ID for each mapped gene")
    if len(set(mapped.values())) != len(mapped) or not a.obs_names.is_unique:
        p.error("Resolve duplicate mapped genes and cell IDs before tokenization")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.output.parent) as tmp:
        stage = Path(tmp) / "run"
        stage.mkdir()
        path = stage / "tokens.jsonl"
        seen = set()
        with args.manifest.open() as f, path.open("w") as out:
            for row in csv.DictReader(f, delimiter="\t"):
                cell = row["cell"]
                if cell in seen:
                    p.error(f"Duplicate manifest cell {cell}")
                seen.add(cell)
                x = matrix[a.obs_names.get_loc(cell) : a.obs_names.get_loc(cell) + 1]
                if sparse.issparse(x):
                    x = x.copy()
                    x.sum_duplicates()
                    ix, raw = x.indices, x.data
                else:
                    raw = np.asarray(x).ravel()
                    ix = range(len(raw))
                if (
                    not np.isfinite(raw).all()
                    or (raw < 0).any()
                    or not np.equal(raw, np.floor(raw)).all()
                ):
                    p.error("Expected nonnegative integer counts")
                counts = {
                    mapped[i]: float(n)
                    for i, n in zip(ix, raw)
                    if i in mapped and n > 0
                }
                record = {k: row[k] for k in META}
                record.update(tokenize(args.model, counts, vocab, medians, limit))
                record["genes"] = [
                    g if g.startswith("<") else canonical[g] for g in record["genes"]
                ]
                if row["role"] == "evaluation":
                    record.update(
                        target_id=row["target_id"],
                        control_ids=row["control_ids"].split(","),
                    )
                    for key in ("target_genes", "control_sets"):
                        if row.get(key):
                            record[key] = json.loads(row[key])
                out.write(json.dumps(record) + "\n")
        if not seen:
            p.error("Empty manifest")
        files = [args.h5ad, args.manifest, args.vocab] + (
            [args.medians] if args.medians else []
        )
        (stage / "tokenization.json").write_text(
            json.dumps(
                {
                    "model": args.model,
                    "layer": args.layer,
                    "gene_column": args.gene_column,
                    "match_column": args.match_column,
                    "context": limit,
                    "cells": len(seen),
                    "mapped_features": len(mapped),
                    "inputs": [{"name": f.name, "sha256": sha(f)} for f in files],
                    "output_sha256": sha(path),
                    "source_sha256": sha(Path(__file__)),
                },
                indent=2,
            )
            + "\n"
        )
        stage.rename(args.output)
