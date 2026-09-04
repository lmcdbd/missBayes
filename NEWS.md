# missBayes 0.99.1

## Bug fixes

* `BayesMissingModel()`, `BayesianMissingAggregate()` and `plotPost()` now
  handle unequal group sizes. The JAGS likelihood looped over `2 * r`
  observations, where `r` was the number of samples in the *first level* of the
  group factor, rather than over the samples actually supplied. Unbalanced
  designs therefore either dropped the trailing samples of the contrast without
  comment (when `2 * r` was smaller than the number of samples) or failed inside
  JAGS for every protein and returned an entirely `NA` results table with no
  error or warning (when it was larger). The likelihood now loops over the
  supplied samples, and the answer no longer depends on the order of
  `levels(groups)`.

* Group membership is now taken from the `groups` factor rather than from
  column position. `sigma_p2params()` and `sigma_jp2params()` previously
  partitioned the columns into consecutive blocks of `r`, so the between- and
  within-group variance priors were estimated from group *mixtures* whenever the
  samples of a group were not contiguous columns - as they are not, for
  instance, when a matrix is ordered by patient. Results no longer depend on
  column order.

* Failures inside `process_row()` are no longer silent. A single warning now
  reports how many proteins failed to fit and the first failure message; the
  per-protein `tryCatch` is kept, so one bad protein still does not abort a run.

* `BayesMissingModel()` validates its inputs up front: `groups` must have one
  entry per column of `values` and no `NA`, at least two groups must be
  populated, and every group named by a contrast must be present in `groups`.

## Notes

* A balanced design whose columns are grouped together - the regime in which
  the previous code was correct - produces numerically identical results to
  0.99.0. `inst/scripts/unbalanced_demo.R` checks this against pre-patch
  reference values, alongside the unbalanced behaviour, using the bundled
  `PXD060201_vignettes.csv`.

* `tests/testthat/` added.

* `sigma_jp2params()` takes the group factor rather than the scalar `r`; its
  signature changed from `(overall_distri, r)` to `(overall_distri, group)`.
  `sigma_p2params()` returns `group_sizes` in place of `r`. Both are internal.

* The stale copies of `model.txt` and `model_logit.txt` at the package root have
  been removed; only the `inst/extdata/` copies are installed and reached by
  `system.file()`.
