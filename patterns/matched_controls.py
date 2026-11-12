"""Matched control genes for a perturbation panel. NumPy only.

Controls are drawn once, by a registered seed, before any perturbation output
exists. Genes are identified by Ensembl ID throughout.
"""
import numpy as np


def build_strata(genes, detect_frac, rank_pct, *, log2_tol=0.5, rank_tol=5.0):
    """Group panel genes so every member is within tolerance of every other member.

    detect_frac: share of analysis cells that hold the gene's token (0 to 1].
    rank_pct: median token rank as a percentile (0 to 100).
    Returns a list of strata, each a list of gene identifiers.
    """
    d = {g: float(np.log2(v)) for g, v in zip(genes, detect_frac)}
    r = {g: float(v) for g, v in zip(genes, rank_pct)}
    strata = []
    for g in sorted(genes, key=lambda k: (d[k], r[k])):
        if strata and all(abs(d[g] - d[m]) <= log2_tol and abs(r[g] - r[m]) <= rank_tol for m in strata[-1]):
            strata[-1].append(g)
        else:
            strata.append([g])
    return strata


def draw_controls(stratum, panel_stats, candidate_stats, *, n_controls=20, seed=0,
                  exclude=(), log2_tol=0.5, rank_tol=5.0):
    """Draw controls that match every member of a stratum.

    panel_stats and candidate_stats map gene -> (detect_frac, rank_pct).
    Returns a list of n_controls genes, or None when fewer candidates exist.
    None means "not estimable: controls". Do not widen the tolerances.
    """
    members = [(np.log2(panel_stats[g][0]), panel_stats[g][1]) for g in stratum]
    banned = set(exclude) | set(panel_stats)
    pool = sorted(
        g for g, (det, rank) in candidate_stats.items()
        if g not in banned and det > 0
        and all(abs(np.log2(det) - md) <= log2_tol and abs(rank - mr) <= rank_tol for md, mr in members)
    )
    if len(pool) < n_controls:
        return None
    rng = np.random.default_rng(seed)
    return [pool[i] for i in sorted(rng.choice(len(pool), size=n_controls, replace=False))]


def donor_adjusted_value(gene_shift, control_shifts, *, min_controls=10):
    """Per-donor control-adjusted value: the gene's mean shift minus the median control shift.

    gene_shift: the gene's mean per-cell shift in this donor, or None if the gene
    is not estimable in this donor. control_shifts: one value per control,
    with None for a control that is not estimable in this donor.
    Returns (value or None, number of qualifying controls).
    """
    usable = [c for c in control_shifts if c is not None]
    if gene_shift is None or len(usable) < min_controls:
        return None, len(usable)
    return float(gene_shift - np.median(usable)), len(usable)
