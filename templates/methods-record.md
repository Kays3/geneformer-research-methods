# Methods record: <study name>, <date>

A short record that travels with every result. A reader should be able to tell what ran without opening the code.

## The five items
| # | Item | Value |
|---|---|---|
| 1 | Code path and commit | |
| 2 | Split (unit, parts, seed, table sha256) | |
| 3 | ISP start and goal states | |
| 4 | Model, checkpoint, weights sha256, precision, library commit | |
| 5 | Fine-tuned against out-of-the-box plots, or the statement that none exists | |

## Registration
File, commit, amendments in force, sign-off.

## Data
Accession, file sha256, cohort table sha256, cells and donors that entered each test.

## Run
- Host class and GPU (no host names), start and end time (UTC).
- Command, as environment variables and a driver name. No absolute paths.
- Gates and their results: no-op, resolvability, classifier, equivalence, output check (values, exact zeros, NaN, median absolute shift).
- GPU-hours: estimate, charged, cap. Reason for any gap.
- Stops and restarts, with exit codes.

## Results
Where the result tables are, their hashes, and which script wrote each number in the report.

## Deviations and changes after the run
Each change to code or text made after the run, with its reason. State whether any script that ran or any result changed.

## What cannot be recomputed from the repository
For example per-cell outputs kept on the compute host.

## Limits
