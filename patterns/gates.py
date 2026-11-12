"""Numerical gates for an in silico perturbation run. NumPy only.

Each function returns a dict with a boolean "passed" and the numbers behind it,
so the result can be written to a JSON file that a launcher checks before the
analysis starts. Thresholds are conventions: register yours before the run.
"""
import numpy as np


def noop_gate(embed, batch) -> dict:
    """An empty edit must change nothing.

    `embed` is a function that runs a forward pass and returns an array
    (cells x dimensions). Call it twice on the same input: two independent
    passes. Comparing one result with itself is zero by algebra and proves
    nothing about stability.
    """
    first = np.asarray(embed(batch), dtype=np.float64)
    second = np.asarray(embed(batch), dtype=np.float64)
    diff = float(np.max(np.abs(first - second)))
    return {"passed": diff == 0.0, "max_abs_diff": diff, "n_cells": int(first.shape[0])}


def output_check(shifts) -> dict:
    """Count values, exact zeros and NaN, and report the median absolute shift."""
    s = np.asarray(shifts, dtype=np.float64).ravel()
    finite = s[np.isfinite(s)]
    return {
        "n_values": int(s.size),
        "n_nan": int(np.isnan(s).sum()),
        "n_exact_zero": int((s == 0).sum()),
        "median_abs_shift": float(np.median(np.abs(finite))) if finite.size else float("nan"),
    }


def resolvability_gate(shifts, state_embeddings, n_donors, *, min_nonzero_share=0.99,
                       min_median_abs=1e-4, min_state_rel_diff=0.01, min_donors=3) -> dict:
    """Can this model and precision resolve a perturbation at all?

    Run on a small sample before a screen and again on its output. A screen
    whose shifts are all zero passes a no-op gate and fails here.
    `state_embeddings` maps a state name to its centroid vector.
    """
    chk = output_check(shifts)
    n = max(chk["n_values"], 1)
    nonzero_share = 1.0 - (chk["n_exact_zero"] + chk["n_nan"]) / n
    names = sorted(state_embeddings)
    rel = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            va = np.asarray(state_embeddings[a], dtype=np.float64)
            vb = np.asarray(state_embeddings[b], dtype=np.float64)
            rel.append(float(np.linalg.norm(va - vb) / max(np.linalg.norm(va), np.linalg.norm(vb))))
    min_rel = min(rel) if rel else float("nan")
    passed = (
        chk["n_nan"] == 0
        and nonzero_share >= min_nonzero_share
        and chk["median_abs_shift"] >= min_median_abs
        and bool(rel) and min_rel >= min_state_rel_diff
        and n_donors >= min_donors
    )
    return {"passed": bool(passed), "nonzero_share": nonzero_share, "min_state_rel_diff": min_rel,
            "n_donors": int(n_donors), **chk}


def embedding_collapse_check(embeddings, strata, *, n_draw=200, seed=0, limit=0.99999) -> dict:
    """Mean pairwise cosine over a seeded draw, stratified so one cell type cannot dominate.

    A fine-tuned model that collapsed maps every cell to nearly the same vector.
    Set `limit` far above the value of a healthy model, measured on the released
    weights with the same draw.
    """
    emb = np.asarray(embeddings, dtype=np.float64)
    strata = np.asarray(strata)
    rng = np.random.default_rng(seed)
    labels = np.unique(strata)
    per = max(n_draw // len(labels), 2)
    idx = np.concatenate([
        rng.choice(np.flatnonzero(strata == lab), size=min(per, int((strata == lab).sum())), replace=False)
        for lab in labels
    ])
    x = emb[idx]
    x = x / np.linalg.norm(x, axis=1, keepdims=True)
    cos = x @ x.T
    upper = cos[np.triu_indices(len(idx), k=1)]
    mean_cos = float(upper.mean())
    return {"passed": mean_cos < limit, "mean_pairwise_cosine": mean_cos, "n_cells": int(len(idx)), "limit": limit}


def prediction_collapse_check(eval_labels, eval_scores, *, min_auc=0.5, min_sd=1e-3) -> dict:
    """After the first epoch: a constant predictor has an AUC near 0.5 and almost no spread."""
    y = np.asarray(eval_labels).astype(bool)
    s = np.asarray(eval_scores, dtype=np.float64)
    pos, neg = s[y], s[~y]
    greater = (pos[:, None] > neg[None, :]).sum()
    equal = (pos[:, None] == neg[None, :]).sum()
    auc = float((greater + 0.5 * equal) / (len(pos) * len(neg)))
    sd = float(s.std())
    return {"passed": auc > min_auc and sd >= min_sd, "auc": auc, "score_sd": sd}
