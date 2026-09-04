#!/usr/bin/env Rscript
#
# Demonstration: BayesMissingModel() with unbalanced group sizes.
#
# Runs entirely from the CSV bundled with the package and prints one PASS/FAIL
# line per check. It works against both the unpatched (0.99.0) and the patched
# package: it detects which is installed and asserts the behaviour appropriate
# to it, so the checks stay runnable - and stay evidence - after the fix.
#
# Reference numbers captured from the unpatched package live in
# inst/extdata/unbalanced_demo_reference.rds. Regenerate them only against an
# unpatched build:
#
#   Rscript inst/scripts/unbalanced_demo.R --write-reference
#
# parallel = FALSE throughout: generate_overdispersed_inits() sets .RNG.name and
# .RNG.seed per chain, so the serial path is exactly reproducible, whereas
# parallel = TRUE workers self-seed and are not.

suppressPackageStartupMessages({
  library(missBayes)
  library(limma)
})

write_reference <- "--write-reference" %in% commandArgs(trailingOnly = TRUE)

# sigma_jp2params() took the scalar `r` before the fix and the group factor
# after it, which is the cleanest signal of which build we are running against.
patched <- "group" %in% names(formals(missBayes:::sigma_jp2params))

reference_path <- if (write_reference) {
  "inst/extdata/unbalanced_demo_reference.rds"
} else {
  system.file("extdata", "unbalanced_demo_reference.rds", package = "missBayes")
}

MCMC <- list(n.adapt = 1000, burn.in = 500, n.iter = 4000, n.chains = 2)
SUMMARY_COLS <- c("Median", "pLtROPE", "pInROPE", "pGtROPE")

## ---------------------------------------------------------------- reporting --

results <- character(0)
report <- function(name, ok, detail) {
  results[[length(results) + 1]] <<- if (isTRUE(ok)) "PASS" else "FAIL"
  cat(sprintf("[%s] %-34s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
}
note <- function(name, detail) cat(sprintf("[----] %-34s %s\n", name, detail))

## ---------------------------------------------------------------- input data --

# The vignette's pre-processing: LFQ intensity columns, log2, zeros to NA,
# median normalisation, drop all-NA rows.
csv_path <- system.file("extdata", "PXD060201_vignettes.csv", package = "missBayes")
pg <- read.csv(csv_path, blank.lines.skip = TRUE)
pg <- pg[, c(1, 558:646)]
pg <- pg[pg$Protein.IDs != "" & !is.na(pg$Protein.IDs), ]
rownames(pg) <- pg$Protein.IDs
pg <- pg[, -1]
log2.df <- log2(pg)
log2.df[is.infinite(as.matrix(log2.df))] <- NA
all89 <- limma::normalizeMedianValues(log2.df)
all89 <- all89[rowSums(!is.na(all89)) > 0, ]

sample_group <- ifelse(grepl("_Tu", colnames(all89)), "Tu", "con")

# The 11 tumour samples with no paired control. The vignette removes them for
# that reason; here they are what makes the design unbalanced.
unpaired_tu <- c("11598_14_Tu", "12452_13_Tu", "147_11_Tu", "1602_13_Tu",
                 "201417387_Tu", "201418285_Tu", "201526492_Tu",
                 "201610390_Tu", "5752_13_Tu", "6877_14_Tu", "9531_12_Tu")
paired <- !(sub("^LFQ.intensity.", "", colnames(all89)) %in% unpaired_tu)

# bal        - the vignette's 39 con / 39 Tu, columns interleaved by patient
# bal_sorted - the same data with the columns grouped together
# all89      - 39 con / 50 Tu; levels order decides which unbalanced regime
bal <- all89[, paired]
bal_group <- sample_group[paired]
sorted <- order(factor(bal_group, levels = c("con", "Tu")))
bal_sorted <- bal[, sorted]
bal_sorted_group <- bal_group[sorted]

cat(sprintf("missBayes %s (%s)\n", packageVersion("missBayes"),
            if (patched) "unbalanced-design fix present" else "unpatched"))
cat(sprintf("bal        %d x %d  (con %d / Tu %d)\n", nrow(bal), ncol(bal),
            sum(bal_group == "con"), sum(bal_group == "Tu")))
cat(sprintf("all89      %d x %d  (con %d / Tu %d)\n\n", nrow(all89), ncol(all89),
            sum(sample_group == "con"), sum(sample_group == "Tu")))

run <- function(values, group, levels) {
  g <- factor(group, levels = levels)
  cm <- limma::makeContrasts(Tu - con, levels = levels(g))
  set.seed(1)
  suppressWarnings(suppressMessages(do.call(
    BayesMissingModel,
    c(list(values = values, groups = g, comparisons = cm, parallel = FALSE), MCMC)
  )))[[1]]
}
fitted_frac <- function(x) mean(!is.na(x$Median))

## ------------------------------------------------------------------- the runs --

cat("running models (about two minutes)...\n\n")
res_bal_sorted <- run(bal_sorted, bal_sorted_group, c("con", "Tu"))
res_bal        <- run(bal,        bal_group,        c("con", "Tu"))
res_lt         <- run(all89,      sample_group,     c("con", "Tu"))
res_gt         <- run(all89,      sample_group,     c("Tu", "con"))

# The contrast Tu - con hands JAGS the 50 Tu samples followed by the 39 con
# samples. With levels c("con", "Tu") the old loop bound was 2 * 39 = 78, so
# positions 79:89 - the last 11 controls - were never read.
ignored_cols <- which(sample_group == "con")[29:39]
target <- "A0A0C4DH31"
perturbed <- all89
perturbed[target, ignored_cols] <- perturbed[target, ignored_cols] + 8
res_lt_pert <- run(perturbed, sample_group, c("con", "Tu"))
target_shift <- abs(res_lt_pert[target, "Median"] - res_lt[target, "Median"])

# Proteins with no missing value at all: the posterior log fold change should
# sit near the difference of the sample means.
complete <- rownames(all89)[rowSums(is.na(all89)) == 0]
observed_lfc <- rowMeans(all89[complete, sample_group == "Tu"]) -
  rowMeans(all89[complete, sample_group == "con"])
posterior_lfc <- res_lt[complete, "Median"]
sanity_median <- stats::median(abs(posterior_lfc - observed_lfc))
sanity_cor <- stats::cor(posterior_lfc, observed_lfc)

if (write_reference) {
  saveRDS(list(
    bal_sorted = res_bal_sorted[, SUMMARY_COLS],
    bal = res_bal[, SUMMARY_COLS],
    fitted_lt = fitted_frac(res_lt),
    fitted_gt = fitted_frac(res_gt),
    target_shift = target_shift,
    sanity_median = sanity_median,
    sanity_cor = sanity_cor
  ), reference_path, version = 2)
  cat("wrote", reference_path, "\n")
  quit(status = 0)
}

reference <- readRDS(reference_path)

## -------------------------------------------------------------------- checks --

# 1. Regression. A balanced design whose columns are already grouped together is
#    the regime the old positional code handled correctly, so the fix must be an
#    exact no-op there.
same <- isTRUE(all.equal(res_bal_sorted[, SUMMARY_COLS], reference$bal_sorted,
                         tolerance = 0))
report("1 balanced regression", same,
       sprintf("39 con / 39 Tu, group-contiguous: %s pre-patch values",
               if (same) "identical to" else "DIFFERS from"))

# 2. Column order. The vignette's balanced matrix interleaves the two groups by
#    patient. The old code took group membership from column position, so it read
#    a con/Tu mixture as each "group"; the fix takes it from the factor.
col_delta <- max(abs(res_bal$Median - res_bal_sorted$Median), na.rm = TRUE)
if (patched) {
  report("2 column-order invariance",
         isTRUE(all.equal(res_bal[, SUMMARY_COLS], res_bal_sorted[, SUMMARY_COLS],
                          tolerance = 0)),
         "interleaved columns give the same answer as grouped columns")
} else {
  report("2 column-order invariance", col_delta > 0.01,
         sprintf("interleaved vs grouped columns differ by up to %.3f (the bug)",
                 col_delta))
}

# 3. Ignored samples. Perturbing only the samples past the old loop bound, for a
#    single protein, must move that protein's posterior - and before the fix it
#    demonstrably did not.
if (patched) {
  report("3 ignored samples now count", target_shift > 0.5,
         sprintf("+8 log2 on the last 11 controls moves %s by %.3f",
                 target, target_shift))
} else {
  report("3 ignored samples are ignored", target_shift < 0.01,
         sprintf("+8 log2 on the last 11 controls moves %s by only %.4f",
                 target, target_shift))
}

# 4. Hard failure. With levels c("Tu", "con") the old loop bound was 2 * 50 = 100
#    against 89 samples, so JAGS raised "Index out of range" for every protein and
#    process_row()'s tryCatch turned the whole table into NA.
if (patched) {
  report("4 larger group first fits", fitted_frac(res_gt) > 0.95,
         sprintf("levels c(\"Tu\",\"con\"): %.1f%% of proteins fitted (was %.1f%%)",
                 100 * fitted_frac(res_gt), 100 * reference$fitted_gt))
} else {
  report("4 larger group first fails", fitted_frac(res_gt) == 0,
         sprintf("levels c(\"Tu\",\"con\"): %.1f%% of proteins fitted, no error raised",
                 100 * fitted_frac(res_gt)))
}

# 5. Level-order invariance. Same data, same contrast, only levels(groups)
#    reordered. This cannot pass while r leaks into the answer.
if (patched) {
  lvl_delta <- max(abs(res_lt$Median - res_gt$Median), na.rm = TRUE)
  report("5 level-order invariance", lvl_delta < 1e-6,
         sprintf("max |difference in Median| = %.2e", lvl_delta))
} else {
  note("5 level-order invariance",
       "not applicable: one level order returns no fitted proteins at all")
}

# 6. Sanity. Only 28 of the 39 controls reached the likelihood before the fix,
#    so the fitted contrast drifted away from the sample means it should track.
if (patched) {
  report("6 fully observed sanity",
         sanity_median < 0.10 && sanity_cor > 0.95,
         sprintf("%d complete proteins: median |posterior - sample| %.3f (was %.3f), r = %.3f (was %.3f)",
                 length(complete), sanity_median, reference$sanity_median,
                 sanity_cor, reference$sanity_cor))
} else {
  note("6 fully observed sanity",
       sprintf("%d complete proteins: median |posterior - sample| %.3f, r = %.3f",
               length(complete), sanity_median, sanity_cor))
}

## -------------------------------------------------------------------- summary --

cat(sprintf("\n%d checks, %d passed, %d failed\n",
            length(results), sum(results == "PASS"), sum(results == "FAIL")))
quit(status = if (any(results == "FAIL")) 1L else 0L)
