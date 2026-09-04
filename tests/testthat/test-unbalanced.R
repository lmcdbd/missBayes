test_that("an unbalanced design fits whichever group is larger", {
  # Before the fix the JAGS likelihood looped to 2 * (size of the first level).
  # With the larger group first that index ran off the end of the data and every
  # protein came back NA; with the smaller group first the trailing samples were
  # silently dropped.
  #
  # A handful of proteins are still skipped for the pre-existing reason that
  # their group mean falls in an intensity bin with too few variance estimates
  # to form a prior, so this asserts "essentially all", not "all". What must
  # hold is that no protein fails inside JAGS - that is what used to wipe out
  # the whole table.
  larger_first <- expect_no_jags_failure(toy_run(make_toy(c(A = 6, B = 3))))
  smaller_first <- expect_no_jags_failure(toy_run(make_toy(c(A = 3, B = 6))))

  expect_gt(mean(!is.na(larger_first$Median)), 0.95)
  expect_gt(mean(!is.na(smaller_first$Median)), 0.95)
})

test_that("samples past the old loop bound inform the posterior", {
  # The contrast "B - A" hands JAGS the 8 B samples followed by the 3 A samples,
  # and the old loop bound of 2 * (size of the first level) = 6 stopped before
  # any A sample was read. Shifting that group must move the answer.
  #
  # Only one protein is perturbed, so the shared hyperparameters are effectively
  # unchanged and the difference is attributable to the likelihood.
  toy <- make_toy(c(A = 3, B = 8))
  ignored <- which(toy$groups == "A")   # positions 9:11, past the old bound of 6
  # A mid-range protein, so that the shift below keeps its group mean inside the
  # intensity bins that carry a variance prior.
  a_means <- rowMeans(toy$values[, toy$groups == "A"], na.rm = TRUE)
  target <- names(which.min(abs(a_means - 17)))

  base <- toy_run(toy)

  shifted <- toy
  shifted$values[target, ignored] <- shifted$values[target, ignored] + 3
  moved <- toy_run(shifted)

  # Raising group A by 3 must drag the B - A posterior down by roughly 3.
  expect_lt(moved[target, "Median"] - base[target, "Median"], -1.5)

  others <- setdiff(rownames(toy$values), target)
  expect_lt(stats::median(abs(moved[others, "Median"] - base[others, "Median"])), 0.05)
})

test_that("posterior log fold changes track the sample means of a fully observed protein", {
  toy <- make_toy(c(A = 5, B = 3))
  complete <- rownames(toy$values)[rowSums(is.na(toy$values)) == 0]
  skip_if(length(complete) < 5, "toy dataset has too few fully observed proteins")

  res <- toy_run(toy)

  observed_lfc <- vapply(complete, function(id) {
    mean(toy$values[id, toy$groups == "B"]) - mean(toy$values[id, toy$groups == "A"])
  }, numeric(1))

  expect_equal(res[complete, "Median"], unname(observed_lfc), tolerance = 0.35)
})

test_that("failed fits are surfaced as a single warning rather than silent NAs", {
  ok <- data.frame(Median = 1)
  bad <- data.frame(Median = NA)
  attr(bad, "fit_error") <- "row 2: Index out of range taking subset of is_observed"

  expect_warning(
    missBayes:::collectFitErrors(list(ok, bad)),
    "1 of 2 proteins failed to fit"
  )
  expect_warning(
    missBayes:::collectFitErrors(list(bad, bad)),
    "All 2 of 2 proteins failed to fit"
  )

  cleaned <- suppressWarnings(missBayes:::collectFitErrors(list(ok, bad)))
  expect_null(attr(cleaned[[2]], "fit_error"))
  expect_silent(missBayes:::collectFitErrors(list(ok, ok)))
})
