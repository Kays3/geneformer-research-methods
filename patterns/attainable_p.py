"""Smallest attainable p values of discrete exact tests, and an exact signed-rank test.

Use before a study is registered: if the smallest attainable p after correction
is above alpha, the design cannot show the effect and the row is "not estimable".
"""
from itertools import product
from math import comb


def floor_sign_test(n_donors: int) -> float:
    """Two-sided floor of a sign or signed-rank test over n donors (no ties, no zeros)."""
    return min(1.0, 2.0 / 2 ** n_donors)


def floor_two_group(n_a: int, n_b: int) -> float:
    """Two-sided floor of an exact rank-sum (or permutation) test on two groups."""
    return min(1.0, 2.0 / comb(n_a + n_b, n_a))


def floor_empirical(n_null: int) -> float:
    """Floor of an empirical p value against n_null matched sets: (k + 1) / (n_null + 1), k = 0."""
    return 1.0 / (n_null + 1)


def holm_floor(floor: float, m: int) -> float:
    """Smallest Holm-adjusted p in a family of m tests that share one floor."""
    return min(1.0, floor * m)


def signed_rank_exact(values) -> float:
    """Exact two-sided Wilcoxon signed-rank p by full enumeration of sign flips.

    Handles ties with mid-ranks and drops zeros. Enumeration costs 2**n, so keep
    n at or below about 20. Use it to check a library "exact" routine when ties
    are possible.
    """
    x = [v for v in values if v != 0]
    n = len(x)
    if n == 0:
        return 1.0
    order = sorted(range(n), key=lambda i: abs(x[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(x[order[j + 1]]) == abs(x[order[i]]):
            j += 1
        mid = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = mid
        i = j + 1
    total = sum(ranks)
    observed = sum(r for r, v in zip(ranks, x) if v > 0)
    dist = abs(observed - total / 2)
    hits = 0
    for signs in product((0, 1), repeat=n):
        w = sum(r for r, s in zip(ranks, signs) if s)
        if abs(w - total / 2) >= dist - 1e-12:
            hits += 1
    return hits / 2 ** n


if __name__ == "__main__":
    for d in (3, 6, 8, 11):
        print(f"sign test, {d} donors: floor {floor_sign_test(d):.5f}")
    print(f"two groups of 6 genes: floor {floor_two_group(6, 6):.5f}")
    print(f"20 matched sets: floor {floor_empirical(20):.4f}")
    print(f"11 donors, Holm family of 7: floor {holm_floor(floor_sign_test(11), 7):.5f}")
