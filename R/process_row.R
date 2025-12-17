#' Summarize posterior distribution of MCMC samples
#'
#' @description
#' A utility function for summarizing MCMC posterior samples. This function
#' was originally written by John K. Kruschke for his book
#' "Doing Bayesian Data Analysis" (2nd edition).
#' @param paramSampleVec A numeric vector of MCMC samples for a single parameter.
#' @param compVal An optional value for comparison (e.g. a null value).
#' @param ROPE A numeric vector of length 2 specifying the lower and upper bounds of the
#'   Region of Practical Equivalence (ROPE). Values within this interval are considered
#'   practically equivalent to the null value (e.g., `c(-0.1, 0.1)` for log fold change).
#' @param credMass A numeric value between 0 and 1 specifying the mass of the
#'   credible interval to compute.
#'
#' @return A named numeric vector of summary statistics.
#' @references Kruschke, J. K. (2015). Doing Bayesian Data Analysis:
#' A Tutorial with R, JAGS, and Stan (2nd ed.). Academic Press.
#' @keywords internal
#'


summarizePost = function( paramSampleVec ,
                          compVal=NULL , ROPE=NULL , credMass=0.95 ) {
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  mcmcEffSz = round( effectiveSize( paramSampleVec ) , 1 )
  names(mcmcEffSz) = NULL
  hdiLim = HDIofMCMC( paramSampleVec , credMass=credMass )
  if ( !is.null(compVal) ) {
    pcgtCompVal = ( 100 * sum( paramSampleVec > compVal )
                    / length( paramSampleVec ) )
  } else {
    compVal=NA
    pcgtCompVal=NA
  }
  if ( !is.null(ROPE) ) {
    pcltRope = ( 100 * sum( paramSampleVec < ROPE[1] )
                 / length( paramSampleVec ) )
    pcgtRope = ( 100 * sum( paramSampleVec > ROPE[2] )
                 / length( paramSampleVec ) )
    pcinRope = 100-(pcltRope+pcgtRope)
  } else {
    ROPE = c(NA,NA)
    pcltRope=NA
    pcgtRope=NA
    pcinRope=NA
  }
  return( c( Mean=meanParam , Median=medianParam , Mode=modeParam ,
             ESS=mcmcEffSz ,
             HDImass=credMass , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2] ,
             CompVal=compVal , PcntGtCompVal=pcgtCompVal ,
             ROPElow=ROPE[1] , ROPEhigh=ROPE[2] ,
             PcntLtROPE=pcltRope , PcntInROPE=pcinRope , PcntGtROPE=pcgtRope ) )
}




#' Perform Bayesian differential analysis for a single protein
#'
#' The core computational function that fits the Bayesian model for a single protein.
#'
#'
#'
#' @param i Integer. The row index of the protein to analyze.
#' @param fixed_data A list containing all preprocessed data, hyperparameters,
#'   and model settings for the analysis. This is returned by `runModel()`.
#' @param model_txt Character string. The JAGS model specification for the truncated normal
#'   missingness mechanism.
#' @param model_logit_txt Character string. The JAGS model specification for the logistic
#'   regression missingness mechanism.
#' @param mcmcDiag Logical. If `TRUE`, includes additional MCMC diagnostics (ESS, MCSE, R-hat)
#'   in the output.
#' @param threshold Numeric. The missing proportion threshold that determines
#'   which missingness mechanism model to use. If the protein's missingness is < `threshold`,
#'   the logistic model is used; otherwise, the truncated normal model is used.
#' @param ROPE A numeric vector of length 2 specifying the lower and upper bounds of the
#'   Region of Practical Equivalence (ROPE). Values within this interval are considered
#'   practically equivalent to the null value (e.g., `c(-0.1, 0.1)` for log fold change).
#' @return A single-row data frame containing summary statistics for the protein's
#'   posterior log fold change distribution. Returns a row of `NA` values if the model fails to run for the protein.
#'
#' @keywords internal


process_row <- function(i, fixed_data, model_txt, model_logit_txt, mcmcDiag, threshold, ROPE ) {
  # Create default summary row (all NA)
  if (mcmcDiag){
    summary_row <- data.frame(
      HDI_Low = NA,
      HDI_High = NA,
      Median = NA,
      pLtCompVal = NA,
      pLtROPE = NA,
      pInROPE = NA,
      pGtROPE = NA,
      ESS = NA,
      MCSE = NA,
      max_rhat = NA
    )
  } else{
    summary_row <- data.frame(
      HDI_Low = NA,
      HDI_High = NA,
      Median = NA,
      pLtCompVal = NA,
      pLtROPE = NA,
      pInROPE = NA,
      pGtROPE = NA
    )
  }

  y <- unlist(fixed_data$subset_data[i, ])

  # Skip if all NA
  if (all(is.na(y))) return(summary_row)

  # Calculate missing proportion
  miss_prop <- mean(is.na(y))
  cp <- fixed_data$row_mins[i]

  # Compute group statistics
  y_numeric <- as.numeric(y)
  y1 <- y_numeric[fixed_data$group1_idx]
  y2 <- y_numeric[fixed_data$group2_idx]

  sample_mean_1 <- mean(y1, na.rm = TRUE)
  sample_mean_1 <- if (is.na(sample_mean_1)) fixed_data$default_mean else sample_mean_1
  sample_mean_2 <- mean(y2, na.rm = TRUE)
  sample_mean_2 <- if (is.na(sample_mean_2)) fixed_data$default_mean else sample_mean_2

  sample_var_1 <- stats::var(y1, na.rm = TRUE)
  sample_var_1 <- if (is.na(sample_var_1)) fixed_data$default_var else sample_var_1
  sample_var_2 <- stats::var(y2, na.rm = TRUE)
  sample_var_2 <- if (is.na(sample_var_2)) fixed_data$default_var else sample_var_2

  # Calculate alpha/beta parameters
  alpha1 <- tryCatch(f_alpha(sample_mean_1, fixed_data$alpha, fixed_data$a, fixed_data$b), error = function(e) NA)
  beta1 <- tryCatch(f_beta(sample_mean_1, fixed_data$beta, fixed_data$a, fixed_data$b), error = function(e) NA)
  alpha2 <- tryCatch(f_alpha(sample_mean_2, fixed_data$alpha, fixed_data$a, fixed_data$b), error = function(e) NA)
  beta2 <- tryCatch(f_beta(sample_mean_2, fixed_data$beta, fixed_data$a, fixed_data$b), error = function(e) NA)

  # Skip if parameters are invalid
  if (any(is.na(c(alpha1, beta1, alpha2, beta2)))) return(summary_row)

  # Generate initial values
  inits_list <- generate_overdispersed_inits(sample_mean_1, sample_mean_2,
                                             sample_var_1, sample_var_2,
                                             n_chains = fixed_data$n_chains)

  # Prepare JAGS data
  is_observed <- ifelse(is.na(y), 0, 1)
  if (miss_prop < threshold) {
    data_jags <- list(
      y = y_numeric,
      is_observed = is_observed,
      group_numeric = fixed_data$group_numeric,
      r = fixed_data$r,
      mu0 = fixed_data$mu_0,
      sigma_1 = fixed_data$sigma_1,
      alpha_p = fixed_data$alpha_p,
      beta_p = fixed_data$beta_p,
      alpha1 = alpha1, beta1 = beta1,
      alpha2 = alpha2, beta2 = beta2,
      gamma0 = fixed_data$gamma0,
      gamma1 = fixed_data$gamma1
    )
    model_to_use <- model_logit_txt
  } else {
    data_jags <- list(
      y = y_numeric,
      is_observed = is_observed,
      group_numeric = fixed_data$group_numeric,
      r = fixed_data$r,
      cp = fixed_data$row_mins[i],
      mu0 = fixed_data$mu_0,
      sigma_1 = fixed_data$sigma_1,
      alpha_p = fixed_data$alpha_p,
      beta_p = fixed_data$beta_p,
      alpha1 = alpha1, beta1 = beta1,
      alpha2 = alpha2, beta2 = beta2
    )
    model_to_use <- model_txt
  }

  # Run JAGS model
  tryCatch({
    # Use textConnection to avoid disk I/O
    jagsModel <- rjags::jags.model(
      file = textConnection(model_to_use),
      data = data_jags,
      inits = inits_list,
      n.chains = fixed_data$n_chains,
      n.adapt = fixed_data$n_adpt,
      quiet = TRUE
    )
    update(jagsModel, n.iter = fixed_data$n_burnin, progress.bar = "none")
    codaSample <- rjags::coda.samples(
      jagsModel,
      variable.names = c('mu1p', 'mu2p'),
      n.iter = fixed_data$n_iter,
      progress.bar = "none"
    )

    # Process results
    posterior_combined <- as.matrix(codaSample)
    difference <- posterior_combined[, "mu1p"] - posterior_combined[, "mu2p"]
    postSummary <- summarizePost(difference, compVal = 0, credMass = 0.95, ROPE = ROPE)

    if (mcmcDiag){
      SD <- stats::sd(difference)
      EffChnLngth <- coda::effectiveSize(difference)
      rhat_values <- coda::gelman.diag(codaSample)$psrf[, 'Point est.']
      summary_row <- data.frame(
        HDI_Low = postSummary["HDIlow"],
        HDI_High = postSummary["HDIhigh"],
        Median = postSummary["Median"],
        pLtCompVal = pmin(postSummary["PcntGtCompVal"] / 100, 1 - postSummary["PcntGtCompVal"] / 100),
        pLtROPE = postSummary["PcntLtROPE"],
        pInROPE = postSummary["PcntInROPE"],
        pGtROPE = postSummary["PcntGtROPE"],
        ESS = EffChnLngth,
        MCSE = SD / sqrt(EffChnLngth),
        max_rhat = max(rhat_values)
      )
    } else{
      summary_row <- data.frame(
        HDI_Low = postSummary["HDIlow"],
        HDI_High = postSummary["HDIhigh"],
        Median = postSummary["Median"],
        pLtCompVal = pmin(postSummary["PcntGtCompVal"] / 100, 1 - postSummary["PcntGtCompVal"] / 100),
        pLtROPE = postSummary["PcntLtROPE"],
        pInROPE = postSummary["PcntInROPE"],
        pGtROPE = postSummary["PcntGtROPE"]
      )
    }

    # Update summary row

  }, error = function(e) {
    message(sprintf("Error in row %d: %s", i, e$message))
  })

  return(summary_row)
}

#' Execute the Bayesian model for differential expression analysis
#' This function conducts the full Bayesian analysis for a specified contrast.
#' @param contrast_data A string specifying which two groups are compared.
#' @param sigma_p2p A list returned by `sigma_p2params`.
#' @param sigma_jp2p A list returned by `sigma_jp2params`.
#' @param zS A vector returned by `zeroState`.
#' @param sigma_1 A numeric returned by `sigma1`.
#' @param default_data A vector returned by `default_val`.
#' @param n.adapt Integer. Number of iterations for the adaptive phase.
#' @param burn.in An integer denoting number of burn-in steps.
#' @param n.iter Integer. Number of iterations after burn-in.
#' @param n.chains Integer. Number of MCMC chains to run.
#' @param threshold Numeric. The missing proportion threshold that determines
#'   which missingness mechanism model to use. If the protein's missingness is < `threshold`,
#'   the logistic model is used; otherwise, the truncated normal model is used.
#' @param parallel Logical. If `TRUE`, the analysis for each protein is
#'   parallelized across multiple CPU cores.
#' @param ROPE A numeric vector of length 2 specifying the lower and upper bounds of the
#'   Region of Practical Equivalence (ROPE). Values within this interval are considered
#'   practically equivalent to the null value (e.g., `c(-0.1, 0.1)` for log fold change).
#' @param mcmcDiag Logical. If `TRUE`, the output includes additional MCMC diagnostics.
#' @return A data frame with one row per protein and columns containing the
#'   summarized posterior statistics. The row names correspond to the protein IDs
#'   from the input data.
#' @seealso \code{\link{BayesMissingModel}}
#' @keywords internal


runModel <- function(contrast_data, sigma_p2p, sigma_jp2p, zS, sigma_1, default_data,
                     n.chains, burn.in, n.iter, n.adapt, parallel, mcmcDiag, threshold, ROPE){



  fixed_data <- list(
    subset_data = contrast_data$subset_data,
    group_subset = factor(contrast_data$group_subset), #Drop unused levels
    groups_to_compare = contrast_data$groups_to_compare,
    row_mins = contrast_data$cp,
    group_numeric = ifelse(contrast_data$group_subset == contrast_data$groups_to_compare[1], 1,
                           ifelse(contrast_data$group_subset == contrast_data$groups_to_compare[2], 2, NA)),

    sigma_1 = sigma_1,

    alpha_p = sigma_p2p$alpha_p,
    beta_p = sigma_p2p$beta_p,
    r = sigma_p2p$r,

    alpha = sigma_jp2p$alpha,
    beta = sigma_jp2p$beta,
    a = sigma_jp2p$a,
    b = sigma_jp2p$b,

    mu_0 = zS['mu_0'],
    gamma0 = zS['gamma0'],
    gamma1 = zS['gamma1'],

    n_chains = n.chains,
    n_burnin = burn.in,
    n_iter = n.iter,
    n_adpt = n.adapt,

    group1_idx = which(contrast_data$group_subset == contrast_data$groups_to_compare[1]),
    group2_idx = which(contrast_data$group_subset == contrast_data$groups_to_compare[2]),

    default_mean = default_data['default_mean'],
    default_var = default_data['default_var']
  )

  model_txt_path <- system.file("extdata", "model.txt", package = "missBayes")
  model_logit_txt_path <- system.file("extdata", "model_logit.txt", package = "missBayes")

  model_txt <- readLines(model_txt_path, warn = FALSE)
  model_logit_txt <- readLines(model_logit_txt_path, warn = FALSE)

  if (parallel){
    n_cores <- parallel::detectCores() - 1

    if (.Platform$OS.type == "unix") {
      # Unix-like systems (Linux/macOS)
      mcmcresults_list_all <- parallel::mclapply(
        X = 1:nrow(fixed_data$subset_data),
        FUN = function(i) {
          process_row(i, fixed_data, model_txt, model_logit_txt)
        },
        mc.cores = n_cores
      )
    } else {
      # Windows systems
      e1 <- environment()
      cl <- parallel::makeCluster(n_cores)

      # Export ALL required objects and functions
      parallel::clusterExport(cl, varlist = ls(e1),
                              envir = e1)

      parallel::clusterEvalQ(cl, {
        requireNamespace("missBayes", quietly = TRUE)
        requireNamespace("rjags", quietly = TRUE)
        requireNamespace("coda", quietly = TRUE)
      })

      mcmcresults_list_all <- parallel::parLapply(
        cl,
        1:nrow(fixed_data$subset_data),
        fun = function(i, fixed_data, model_txt, model_logit_txt, mcmcDiag, threshold, ROPE) {
          process_row(i, fixed_data, model_txt, model_logit_txt, mcmcDiag, threshold, ROPE)
        },
        fixed_data = fixed_data,
        model_txt = model_txt,
        model_logit_txt = model_logit_txt,
        mcmcDiag = mcmcDiag,
        threshold = threshold,
        ROPE = ROPE
      )
      parallel::stopCluster(cl)
    }

    # Combine results
    mcmcresults_df <- do.call(rbind, mcmcresults_list_all)
    rownames(mcmcresults_df) <- rownames(fixed_data$subset_data)
  } else{

    mcmcresults_list_all <- lapply(
      X = 1:nrow(fixed_data$subset_data),
      FUN = function(i) {
        process_row(i, fixed_data, model_txt, model_logit_txt, mcmcDiag, threshold, ROPE)
      }
    )
    mcmcresults_df <- do.call(rbind, mcmcresults_list_all)
    rownames(mcmcresults_df) <- rownames(fixed_data$subset_data)
  }

  return(mcmcresults_df)
}



