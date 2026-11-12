# Code patterns

Small, generic examples written for this repository. They show the shape of a check; they are not the code that produced any result quoted in `docs/`. They depend only on Python 3.9+ and NumPy (`gates.py`, `matched_controls.py`) or the standard library (`attainable_p.py`). `run_guard.sh` is a bash skeleton.

| File | What it shows |
|---|---|
| `attainable_p.py` | The smallest p a discrete exact test can reach, and an exact signed-rank test by enumeration |
| `gates.py` | No-op gate, resolvability gate, output check, embedding-collapse check, prediction-collapse check |
| `matched_controls.py` | Strata by detection and token rank, a seeded control draw, per-donor control adjustment |
| `run_guard.sh` | A resumable driver: required environment variables, hash pins, done files, a GPU-second cap, distinct exit codes |
| `test_patterns.py` | Checks for the three Python files: `python test_patterns.py` |

Thresholds in these files are the conventions described in `docs/CRITERIA.md`. Set your own and register them before you run.
