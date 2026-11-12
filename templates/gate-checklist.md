# Gate checklist

Three independent reviews before a result is merged or quoted. The author does not merge. A verdict is tied to the head commit and to a hash of the description text; a verdict on an older version does not count. No push between the gate request and the last verdict.

## Science gate (design and reading)
- [ ] The five items are stated: code and commit; split; ISP start and goal states; model and checkpoint; fine-tuned against out-of-the-box plots or the statement that none exists.
- [ ] A registration exists, was committed before the output, and amendments are append-only.
- [ ] The primary test, family and predicted direction in the report match the registration.
- [ ] The unit of inference is the donor (or the limit to one sample is stated).
- [ ] Each effect is compared with a matched null, not with zero.
- [ ] Both operations were scored on the same cells; concordance is read against a random-gene baseline.
- [ ] Gates ran and passed: no-op, resolvability, classifier, equivalence, stability. A failed gate is reported, not worked around.
- [ ] Every registered gene has a row and a status from the fixed vocabulary; "no effect" does not appear.
- [ ] Claims match what the test compares. Ambient and contamination rules were applied.
- [ ] Post hoc analyses are labelled. Limits are stated.
- [ ] A cited clause of a standard is cited for the purpose it serves (a false-positive control does not qualify a negative).

## Numbers gate
- [ ] Every number in the text, tables and figures was recomputed from the committed result tables.
- [ ] Counts in prose match the tables (donors, cells, genes, calls, files).
- [ ] p values are attainable for the stated n; adjusted values use the registered m.
- [ ] Exact zeros and NaN are counted by event.
- [ ] GPU-hours match the ledger; estimate, charge and cap are all given.
- [ ] Input and code hashes in the record match the files.
- [ ] After a fix: the old value is gone in every form, rounded forms included, in the text and in the description.
- [ ] The document was built from the current source (compare strings that exist only in the last change).

## Privacy gate (full read)
- [ ] No personal names, e-mail addresses or affiliations beyond what is intended.
- [ ] No host names, network addresses, user names or home paths; paths come from environment variables.
- [ ] No donor or patient identifiers.
- [ ] No number, gene list, figure or cell count from unpublished data in a public document.
- [ ] No credentials, tokens or key names; commit messages and titles included.
- [ ] Earlier commits on the branch were checked as well as the tree.

## Re-gate
- A prose-only fix of a few lines may be checked as a delta. A change to any script, table or number gets a full read.
- The verdict states the head commit and description hash it was taken on, re-read as the last action before writing.
