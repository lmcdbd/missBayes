test_that("a balanced, group-contiguous run is unchanged by the unbalanced-design fix", {
  # balanced-reference.rds was produced by missBayes 0.99.0, before groups were
  # derived from the group factor rather than from column position. A balanced
  # design whose columns are already grouped together is the regime in which the
  # old positional code was correct, so the fix must be an exact no-op there.
  #
  # The posterior summaries come out of an MCMC run, so they depend on the JAGS
  # version as well as on this package. Comparing to a few significant figures
  # keeps a genuine change in the model detectable while leaving room for the
  # last-bit differences a JAGS point release can introduce.
  reference <- readRDS(test_path("balanced-reference.rds"))

  result <- toy_run(make_toy(c(A = 4, B = 4)))

  expect_identical(dimnames(result), dimnames(reference))
  expect_identical(result$Convergence, reference$Convergence)

  numeric_cols <- names(reference)[vapply(reference, is.numeric, logical(1))]
  expect_equal(result[numeric_cols], reference[numeric_cols], tolerance = 1e-6)
})
