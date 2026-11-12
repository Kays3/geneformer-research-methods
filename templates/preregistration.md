# Pre-registration: <study name>

Status: DRAFT | REGISTERED at commit <hash> on <date, UTC>.
After registration nothing above the amendment line changes. Changes are appended as dated amendments.

## 0. What existed when this was written
State which data, tables and outputs the author had seen. "No perturbation output exists" is the usual line; if anything else is true, say what.

## 1. Question and hypotheses
One line per hypothesis, each with its predicted direction as a separate clause.

## 2. The five items
1. Code: path and commit of the run package.
2. Split: unit (donor), parts or folds, seed, sha256 of the split table. For a zero-shot run: "no training split".
3. ISP states: start state, goal state, how the goal is represented.
4. Model: name, checkpoint, weights sha256, precision, library commit.
5. Fine-tuned against out-of-the-box: which plots will show both, or "no fine-tuned run; limit: ...".

## 3. Data and cohort
Source accession, file sha256, selection rule (structural only), cells per donor, exclusions, what the cohort cannot represent.

## 4. Perturbation and per-donor value
Operations (delete, overexpress), the cell set each applies to (the same token-positive cells for both), the per-cell readout, the per-donor value.

## 5. Eligibility
Minimum token-positive cells per donor, minimum donors, and the status given when a gene falls short.

## 6. Controls and nulls
Matching variables and tolerances, number of controls, seed, exclusions. The random-gene or matched-set null: frame, size, seed.

## 7. Tests
Primary test, sidedness, alpha, family and its size m, correction, smallest attainable adjusted p at the planned n, ties rule.

## 8. Gates
No-op, resolvability, classifier, cross-host or repeat equivalence, stability. For each: the threshold, when it runs, and what a failure does.

## 9. Outcome vocabulary
The fixed list of statuses and the rule that maps results to each. "No effect" is not a status.

## 10. Readings fixed in advance
| Hypothesis | Reading if met | Reading if not |
|---|---|---|

## 11. Sensitivities
Each with the statement that it never changes a status.

## 12. Stop conditions and budget
GPU-hour estimate from measured rates, hard cap, time window, exit codes, what happens on a projected overrun.

## 13. Limits stated in advance
What the study cannot show.

## 14. Sign-off
Reviewer role, date, and the hashes the sign-off covers.

---

## Amendments (append only)

### Amendment <n>, <date and time, UTC>
- What existed when this was written (outputs read, runs finished).
- What changes, and the reason.
- Whether it depends on anything already observed.
- Which feasibility checks were re-run because of it.
- Check before push: zero deleted lines in the diff; the first N bytes of this file hash to the previous version (sha256 <hash>).

### Deviation log
Dated entries. State whether each is a protocol deviation or a record of a miss (for example a cost above its prediction), and whether it bears on any reading.
