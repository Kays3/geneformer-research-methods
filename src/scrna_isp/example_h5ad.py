"""Generate a small synthetic AnnData and manifest for testing installed checkpoints."""

import argparse
import csv
import tempfile
from pathlib import Path

from .cli import META


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    if a.output.exists():
        p.error("Choose a new output")
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.sparse import csr_matrix

    with (Path(__file__).with_name("example") / "counts.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    genes = ["ENSG00000075624", "ENSG00000111640", "ENSG00000166710", "ENSG00000168610"]
    var = pd.DataFrame(
        {"ensembl_id": genes, "human_symbol": ["ACTB", "GAPDH", "B2M", "STAT3"]},
        index=genes,
    )
    counts = csr_matrix(np.array([[int(r["GENE_" + c]) for c in "ABCD"] for r in rows]))
    data = ad.AnnData(
        X=counts, obs=pd.DataFrame(index=[r["cell"] for r in rows]), var=var
    )
    data.layers["counts"] = counts.copy()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.output.parent) as tmp:
        stage = Path(tmp) / "example"
        stage.mkdir()
        data.write_h5ad(stage / "counts.h5ad")
        with (stage / "manifest.tsv").open("w", newline="") as f:
            w = csv.DictWriter(
                f, fieldnames=[*META, "target_id", "control_ids"], delimiter="\t"
            )
            w.writeheader()
            for r in rows:
                w.writerow(
                    {
                        **{k: r[k] for k in META},
                        "target_id": genes[0],
                        "control_ids": genes[1],
                    }
                )
        stage.rename(a.output)
