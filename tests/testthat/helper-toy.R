# Deterministic toy datasets shared by the tests.
#
# `make_toy()` builds a protein x sample matrix with a realistic
# intensity-dependent missingness pattern, so that zeroState()'s logistic fit
# and the per-bin variance priors are all well defined. Group sizes are
# arbitrary, so the same generator produces balanced and unbalanced designs.

make_toy <- function(sizes = c(A = 4, B = 4), n_prot = 120, seed = 42,
                     effect = 1, interleave = FALSE) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")

  groups <- factor(rep(names(sizes), times = sizes), levels = names(sizes))
  n <- sum(sizes)

  # Protein-level baselines spread across the intensity range, plus a group
  # effect on the second group and per-sample noise.
  baseline <- stats::runif(n_prot, 14, 24)
  shift <- ifelse(groups == names(sizes)[1], 0, effect)
  values <- outer(baseline, shift, "+") +
    matrix(stats::rnorm(n_prot * n, sd = 0.6), nrow = n_prot)

  # Intensity-dependent (MNAR) missingness.
  p_obs <- stats::plogis(values - 15)
  values[matrix(stats::runif(n_prot * n), nrow = n_prot) > p_obs] <- NA
  values <- values[rowSums(!is.na(values)) > 1, , drop = FALSE]

  rownames(values) <- sprintf("Protein_%03d", seq_len(nrow(values)))
  colnames(values) <- sprintf("%s_%02d", groups, seq_len(n))

  if (interleave) {
    ord <- order(sequence(sizes))   # A_01 B_05 A_02 B_06 ...
    values <- values[, ord, drop = FALSE]
    groups <- groups[ord]
  }

  list(values = values, groups = groups)
}

# The model is fitted by JAGS through rjags, which needs the JAGS binary as
# well as the R package; without it the fitting tests are skipped rather than
# reported as failures.
skip_without_jags <- function() {
  testthat::skip_if_not_installed("rjags")
  testthat::skip_if(
    inherits(try(rjags::jags.version(), silent = TRUE), "try-error"),
    "JAGS is not available"
  )
}

# A protein whose group mean falls outside the intensity range covered by the
# variance priors is skipped before JAGS is reached, and collectFitErrors()
# reports it. A failure inside JAGS - the regression this suite guards against -
# reads differently, so assert on the message rather than on silence.
expect_no_jags_failure <- function(expr) {
  seen <- character()
  result <- withCallingHandlers(expr, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  unexpected <- grep("No (alpha|beta) value available", seen,
                     invert = TRUE, value = TRUE)
  testthat::expect_identical(unexpected, character(0))
  result
}

toy_contrast <- function(groups) {
  limma::makeContrasts(contrasts = "B - A", levels = levels(groups))
}

# MCMC settings small enough to keep the suite quick but long enough for the
# posterior summaries to be stable.
toy_run <- function(toy, ...) {
  skip_without_jags()
  set.seed(1)
  suppressMessages(BayesMissingModel(
    values = toy$values, groups = toy$groups,
    comparisons = toy_contrast(toy$groups),
    parallel = FALSE, n.adapt = 500, burn.in = 500, n.iter = 4000,
    ...
  ))[[1]]
}
