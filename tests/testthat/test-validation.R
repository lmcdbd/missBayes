test_that("BayesMissingModel rejects a groups vector that does not match values", {
  toy <- make_toy(c(A = 4, B = 4))
  cm <- toy_contrast(toy$groups)

  expect_error(
    BayesMissingModel(toy$values, toy$groups[-1], cm, parallel = FALSE),
    "one entry per column"
  )
})

test_that("BayesMissingModel rejects fewer than two populated groups", {
  toy <- make_toy(c(A = 4, B = 4))
  cm <- toy_contrast(toy$groups)
  one_group <- factor(rep("A", ncol(toy$values)), levels = c("A", "B"))

  expect_error(
    BayesMissingModel(toy$values, one_group, cm, parallel = FALSE),
    "at least two non-empty groups"
  )
})

test_that("BayesMissingModel rejects a contrast naming an absent group", {
  toy <- make_toy(c(A = 4, B = 4))
  three <- factor(c(as.character(toy$groups), "C"), levels = c("A", "B", "C"))
  cm <- limma::makeContrasts(contrasts = "C - A", levels = levels(three))

  expect_error(
    BayesMissingModel(toy$values, toy$groups, cm, parallel = FALSE),
    "absent from 'groups'"
  )
})

test_that("BayesMissingModel rejects NA group assignments", {
  toy <- make_toy(c(A = 4, B = 4))
  cm <- toy_contrast(toy$groups)
  with_na <- toy$groups
  with_na[1] <- NA

  expect_error(
    BayesMissingModel(toy$values, with_na, cm, parallel = FALSE),
    "must not contain NA"
  )
})

test_that("BayesMissingModel rejects a contrast that is not a pairwise comparison", {
  # setupContrasts() resolves a contrast to one "+1" group and one "-1" group,
  # so an averaged contrast has to be caught up front rather than failing on a
  # zero-length subscript deep inside it.
  toy <- make_toy(c(A = 4, B = 4))
  three <- factor(rep(c("A", "B", "C"), length.out = ncol(toy$values)),
                  levels = c("A", "B", "C"))
  cm <- limma::makeContrasts(contrasts = "(A + B) / 2 - C", levels = levels(three))

  expect_error(
    BayesMissingModel(toy$values, three, cm, parallel = FALSE),
    "not a pairwise comparison"
  )
})

test_that("plotPost applies the same groups validation as BayesMissingModel", {
  toy <- make_toy(c(A = 4, B = 4))

  expect_error(
    plotPost(toy$values, rownames(toy$values)[1], toy$groups[-1], "B - A"),
    "one entry per column"
  )
})

test_that("valid unbalanced input passes validation and returns groups as a factor", {
  toy <- make_toy(c(A = 5, B = 3))
  cm <- toy_contrast(toy$groups)

  expect_identical(
    missBayes:::validateGroups(toy$values, as.character(toy$groups), cm),
    factor(as.character(toy$groups))
  )
})
