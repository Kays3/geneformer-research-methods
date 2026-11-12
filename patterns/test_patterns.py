"""Checks for the pattern files. Run: python test_patterns.py"""
import numpy as np

from attainable_p import floor_empirical, floor_sign_test, floor_two_group, signed_rank_exact
from gates import (embedding_collapse_check, noop_gate, output_check, prediction_collapse_check,
                   resolvability_gate)
from matched_controls import build_strata, donor_adjusted_value, draw_controls

rng = np.random.default_rng(1)

# attainable p
assert floor_sign_test(11) == 2 / 2048
assert abs(floor_two_group(6, 6) - 2 / 924) < 1e-15
assert floor_empirical(150) == 1 / 151
assert signed_rank_exact([1, 2, 3]) == 0.25                      # three pairs cannot reach 0.05
assert signed_rank_exact([0.5 + i for i in range(11)]) == 2 / 2048  # all donors agree: the floor
assert signed_rank_exact([1, -1, 2, -2]) == 1.0                  # ties, balanced

# no-op gate: deterministic passes, noisy fails
x = rng.normal(size=(8, 16))
assert noop_gate(lambda b: b * 2.0, x)["passed"]
assert not noop_gate(lambda b: b + rng.normal(scale=1e-6, size=b.shape), x)["passed"]

# resolvability: all-zero shifts fail although a no-op gate would pass
states = {"a": np.ones(16), "b": np.ones(16) * 1.05}
assert not resolvability_gate(np.zeros(500), states, n_donors=5)["passed"]
assert resolvability_gate(rng.normal(scale=5e-3, size=500), states, n_donors=5)["passed"]
assert not resolvability_gate(rng.normal(scale=5e-3, size=500), {"a": np.ones(16), "b": np.ones(16)}, 5)["passed"]
assert output_check([0.0, 1.0, float("nan")]) == {"n_values": 3, "n_nan": 1, "n_exact_zero": 1,
                                                  "median_abs_shift": 0.5}

# collapse checks
strata = np.repeat(["t", "b", "m", "n"], 100)
healthy = rng.normal(size=(400, 32))
collapsed = np.ones((400, 32)) + rng.normal(scale=1e-7, size=(400, 32))
assert embedding_collapse_check(healthy, strata)["passed"]
assert not embedding_collapse_check(collapsed, strata)["passed"]
y = rng.integers(0, 2, size=300)
assert prediction_collapse_check(y, y + rng.normal(scale=0.5, size=300))["passed"]
assert not prediction_collapse_check(y, np.full(300, 0.37))["passed"]

# matched controls
genes = ["ENSG01", "ENSG02", "ENSG03"]
strata_ = build_strata(genes, [0.50, 0.55, 0.05], [40.0, 42.0, 80.0])
assert strata_ == [["ENSG03"], ["ENSG01", "ENSG02"]]
panel = {"ENSG01": (0.50, 40.0), "ENSG02": (0.55, 42.0)}
cands = {f"ENSC{i:03d}": (0.52, 41.0) for i in range(30)}
cands["ENSG01"] = (0.50, 40.0)
a = draw_controls(["ENSG01", "ENSG02"], panel, cands, seed=7)
assert len(a) == 20 and "ENSG01" not in a and a == draw_controls(["ENSG01", "ENSG02"], panel, cands, seed=7)
assert draw_controls(["ENSG01", "ENSG02"], panel, dict(list(cands.items())[:10]), seed=7) is None
assert donor_adjusted_value(0.010, [0.002] * 12 + [None] * 8) == (0.008, 12)
assert donor_adjusted_value(0.010, [0.002] * 9 + [None] * 11) == (None, 9)
print("all pattern checks passed")
