test_that("results do not depend on the order of levels(groups)", {
  # The old likelihood loop ran to 2 * (size of the first group level), so the
  # answer changed - or vanished entirely - when the levels were reordered.
  toy <- make_toy(c(A = 5, B = 3))

  res_ab <- toy_run(toy)

  toy_ba <- toy
  toy_ba$groups <- factor(toy$groups, levels = c("B", "A"))
  res_ba <- toy_run(toy_ba)

  expect_equal(res_ba$Median, res_ab$Median, tolerance = 1e-6)
  expect_equal(res_ba$pInROPE, res_ab$pInROPE, tolerance = 1e-6)
  expect_identical(is.na(res_ba$Median), is.na(res_ab$Median))
})

test_that("results do not depend on the column order of values", {
  # Group membership comes from the factor, so interleaving the samples of the
  # two groups must not change anything.
  contiguous <- toy_run(make_toy(c(A = 4, B = 4)))
  interleaved <- toy_run(make_toy(c(A = 4, B = 4), interleave = TRUE))

  expect_equal(interleaved$Median, contiguous$Median, tolerance = 1e-6)
})
