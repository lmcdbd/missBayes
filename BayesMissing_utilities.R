library(rjags)
library(runjags)
library(data.table)
library(bnlearn)
library(matrixStats)
library(coda)
library(ggplot2)
library(dplyr)
library(limma)
library(QFeatures)
source("DBDA2E-utilities.R")  # this file is from Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition, which provides some very useful functions

# 2.1 calculate mu0, gamma0 and gamma1 
zeroState <- function(overall_distri){

  overall_info <- hist(as.matrix(overall_distri), breaks = 100, plot = FALSE)
  missing_prop <- sum(is.na(overall_distri)) / (ncol(overall_distri) * nrow(overall_distri)) 
  overall_info$density <- overall_info$density * (1 - missing_prop)
  mids <- overall_info$mids
  dens <- overall_info$density
  mu_0 <- mids[which.max(dens)]
  # empty vector to store estimated missing density per bin
  missing_density <- rep(0, length(overall_info$density))

  for (i in seq_along(mids)) {
    if (mids[i] < mu_0) {
      mirror_i <- which.min(abs(mids - (2 * mu_0 - mids[i])))
      expected <- dens[mirror_i]
      observed <- dens[i]
      missing_density[i] <- max(0, expected - observed)
    }
  }

  # estimate observed probability per bin
  observed_count <- dens
  total_count <- observed_count + missing_density
  obs_prob <- observed_count / total_count
  
  # fit logistic regression
  df <- data.frame(intensity = mids, p_obs = obs_prob)
  # filter valid values
  df <- df[is.finite(df$p_obs) & df$p_obs > 0 & df$p_obs < 1, ]  
  
  logit_obs <- glm(p_obs ~ intensity, data = df, family = binomial(link = "logit"))
  summary(logit_obs)
  
  # extract coefficients
  gamma0 <- coef(logit_obs)[1]
  gamma1 <- coef(logit_obs)[2]
  zS <- c(mu_0, gamma0, gamma1) 
  names(zS) <- c('mu_0', 'gamma0', 'gamma1')
  
  return (zS)
}

# 2.2 sigma1 is the standard deviation of the protein means
sigma1 <- function(overall_distri){
  Protein_means <- rowMeans(overall_distri, na.rm = TRUE)
  sigma_1 <- sd(Protein_means, na.rm = TRUE)
  
  return (sigma_1)
}

# 2.3 sigma_p^2 is from a inverse-gamma distribution with shape of alpha and rate of beta
sigma_p2params <- function(overall_distri, group){
  
  r <- unname(table(group)[1])  
  n_groups <- nlevels(group) # Number of replicates per group
  
  group_means <- sapply(1:n_groups, function(i) {
    group_indices <- (r * (i - 1) + 1):(r * i)
    rowMeans(overall_distri[,group_indices])
  })
  
  n <- n_groups # n is the number of groups
  between_group_vars_all <- var(group_means, na.rm = TRUE)
  
  mean_var <- mean(between_group_vars_all, na.rm = TRUE)
  var_var <- var(between_group_vars_all, na.rm = TRUE)
  alpha_p <- mean_var^2 / n*var_var + 2
  beta_p <- mean_var * (alpha_p - 1)
  
  sigma_p2params <- list(alpha_p[1,1], beta_p[1,1], as.numeric(r), group_means)
  names(sigma_p2params) <- c('alpha_p', 'beta_p', 'r', 'group_means')
  
  return (sigma_p2params)
}

# 2.4 sigma_jp^2 is from a inverse-gamma distribution with shape of alpha and rate of beta which are dependent of mu_jp
sigma_jp2params <- function(overall_distri, r){
  group_means_all <- c()
  group_vars_all <- c()
  for (i in 1:(ncol(overall_distri) %/% r)) {
    group_i <- overall_distri[, (r * (i - 1) + 1):(r * i)]
    group_i_mean <- rowMeans(group_i, na.rm = TRUE)
    group_i_var   <- rowVars(as.matrix(group_i), na.rm = TRUE)
    group_means_all <- c(group_means_all, group_i_mean)
    group_vars_all  <- c(group_vars_all, group_i_var)
  }
  ##### calculate alpha and beta for each interval
  a <- as.integer(min(group_means_all, na.rm = TRUE))
  b <- as.integer(max(group_means_all, na.rm = TRUE))
  alpha <- c()
  beta <- c()
  bin_centers <- c()
  for (x in a:b) {  
    selected_vars <- group_vars_all[group_means_all >= x & group_means_all <= x + 1]
    selected_vars <- selected_vars[!is.na(selected_vars)]
    bin_centers_x <- mean(c(x,x+1))
    bin_centers <- c(bin_centers, bin_centers_x)
    if(length(selected_vars) > 1) {  
      mean_vars <- mean(selected_vars)
      var_vars <- var(selected_vars)
      
      alpha_x <- (2 + (mean_vars^2 / r * var_vars))
      beta_x <- mean_vars * (alpha_x - 1)
      
      alpha <- c(alpha, alpha_x)  
      beta <- c(beta, beta_x)     
    } else {
      alpha <- c(alpha, NA)  
      beta <- c(beta, NA)
    }
  }
  params <- c(list(alpha), list(beta), a, b)
  names(params) <- c('alpha', 'beta', 'a', 'b')
  return (params)
}

# define the function f_alpha and f_beta that in a certain interval, return alpha and beta for each mu in a certain interval 
f_alpha <- function(mu_val, alpha, a, b) {
  bin_breaks <- a:(b + 1)  # From a to b+1 so that each bin is [x, x+1]
  idx <- findInterval(mu_val, bin_breaks, rightmost.closed = TRUE)
  if (!is.na(alpha[idx])) {
    return(alpha[idx])
  } else {
    stop(paste("No alpha value available for mu =", mu_val))
  }
}

f_beta <- function(mu_val, beta, a, b) {
  bin_breaks <- a:(b + 1)  # From a to b+1 so that each bin is [x, x+1]
  idx <- findInterval(mu_val, bin_breaks, rightmost.closed = TRUE)
  if (!is.na(beta[idx])) {
    return(beta[idx])
  } else {
    stop(paste("No beta value available for mu =", mu_val))
  }
}

# 2.5 specify all comparisons, filter the dataset, only keep proteins with missing values
setupContrasts <- function(overall_distri, comparisons, contrast, group, filter4NAs){
  # Split the comparison string into two group names
  groups_to_compare <- names(comparisons[,contrast])
  # Check that both groups exist in the factor
  if (!all(groups_to_compare %in% names(comparisons[,contrast]))) {
    stop("One or both groups in the comparison do not exist in 'group' levels.")
  }
  
  # Get indices for each group
  group1_indices <- which(group == groups_to_compare[1])
  group2_indices <- which(group == groups_to_compare[2])
  # Combine indices if needed
  selected_indices <- c(group1_indices, group2_indices)
  # Subset the data
  subset_data <- overall_distri[, selected_indices]
  group_subset <- group[selected_indices]
  # only keep the rows with missing values
  if (filter4NAs){
    subset_data <- subset_data[rowSums(is.na(subset_data)) > 0 & rowSums(is.na(subset_data)) < ncol(subset_data), ]
    rows_with_na <- rownames(subset_data)
    overall_distri <- overall_distri[rows_with_na, , drop = FALSE]
  }
  # get row minimums from overall distribution only for the rows present in filtered subset
  row_mins <- apply(overall_distri, 1, min, na.rm = TRUE)
  
  output <- list(subset_data, group_subset, groups_to_compare, row_mins)
  names(output) <- c('subset_data', 'group_subset', 'groups_to_compare', 'cp')
  
  return (output)
}

# 3. model ----
# model based on cutoffs
runModel <- function(overall_distri, contrast_data, sigma_p2p, sigma_jp2p, zS, sigma_1){
  modelString <- "
    model {
      # Likelihood (model.txt)
      for (i in 1:(2 * r)) {
        is_observed[i] ~ dinterval(y[i], cp)
        y[i] ~ dnorm(group_mean[group_numeric[i]], group_tau[group_numeric[i]])
      }
    
      group_mean[1] <- mu1p
      group_mean[2] <- mu2p
      group_tau[1] <- tau1p
      group_tau[2] <- tau2p
      
      # Group means
      mu1p ~ dnorm(mu_p, tau_p)
      mu2p ~ dnorm(mu_p, tau_p)
      
      # Group-specific data precision (depends on mu)
      tau1p ~ dgamma(alpha1, beta1)
      tau2p ~ dgamma(alpha2, beta2)
      
      # Priors for mu_p and tau_p
      mu_p ~ dnorm(mu0, tau1)
      tau_p ~ dgamma(alpha_p, beta_p)
      tau1 <- 1 / (sigma_1^2)
    }

  # model based on intensity-dependent missing probabilities (model_logit.txt)
  model {
    # Likelihood
    for (i in 1:(2 * r)) {
      # Intensity model (normal)
      y[i] ~ dnorm(group_mean[group_numeric[i]], group_tau[group_numeric[i]])
      
      # Missingness model: probability that y[i] is observed
      logit(p[i]) <- gamma0 + gamma1 * y[i]
      is_observed[i] ~ dbern(p[i])
    }
    
    group_mean[1] <- mu1p
    group_mean[2] <- mu2p
    group_tau[1] <- tau1p
    group_tau[2] <- tau2p
    
    # Group means
    mu1p ~ dnorm(mu_p, tau_p)
    mu2p ~ dnorm(mu_p, tau_p)
    
    # Group-specific data precision (depends on mu)
    tau1p ~ dgamma(alpha1, beta1)
    tau2p ~ dgamma(alpha2, beta2)
    
    # Priors for mu_p and tau_p
    mu_p ~ dnorm(mu0, tau1)
    tau_p ~ dgamma(alpha_p, beta_p)
    tau1 <- 1 / (sigma_1^2)
    
  }
  "
  
  subset_data <- contrast_data$subset_data
  group_subset <- contrast_data$group_subset
  groups_to_compare <- contrast_data$groups_to_compare
  row_mins <- contrast_data$cp
  
  alpha_p <- sigma_p2p$alpha_p
  beta_p <- sigma_p2p$beta_p
  r <- sigma_p2p$r
  
  alpha <- sigma_jp2p$alpha
  beta <- sigma_jp2p$beta
  a <- sigma_jp2p$a
  b <- sigma_jp2p$b
  
  mu_0 <- zS['mu_0']
  gamma0 <- zS['gamma0']
  gamma1 <- zS['gamma1']
  
  # 4. data  ----
  # loop through entire dataset
  group_numeric <- ifelse(group_subset == groups_to_compare[1], 1,
                          ifelse(group_subset == groups_to_compare[2], 2, NA))

  mcmcresults_list <- list()
  for (i in 1: nrow(subset_data)){
    y <- unlist(subset_data[i, ])

    # initialize summary with NA
    summary_row <- data.frame(
      HDI_Low = NA,
      HDI_High = NA,
      Median = NA,
      pLtCompVal = NA, 
      pLtROPE = NA,
      pInROPE = NA,
      pGtROPE = NA,
      DE = NA
    )
    if (all(is.na(y))) {
      mcmcresults_list[[i]] <- summary_row
      next
    }
    
    is_observed <- ifelse(is.na(y), 0, 1)        # 1 = observed, 0 = missing/censored
    miss_prop <- mean(is.na(y))
    cp <- row_mins[i]   #min(y, na.rm = TRUE)      # Cutoff point 
    default_mean <- quantile(unlist(overall_distri), probs = 0.05, na.rm = TRUE)    # If a group is entirely missing, set an default mean
    sample_mean_1 <- mean(as.numeric(y[group_subset == groups_to_compare[1]]), na.rm = TRUE)
    sample_mean_1 <- if (is.na(sample_mean_1)) default_mean else sample_mean_1
    sample_mean_2 <- mean(as.numeric(y[group_subset == groups_to_compare[2]]), na.rm = TRUE)
    sample_mean_2 <- if (is.na(sample_mean_2)) default_mean else sample_mean_2
    alpha1 <- tryCatch({
      f_alpha(sample_mean_1, alpha, a, b)
    }, error = function(e) NA)
    beta1 <- tryCatch({
      f_beta(sample_mean_1, beta, a, b)
    }, error = function(e) NA)
    alpha2 <- tryCatch({
      f_alpha(sample_mean_2, alpha, a, b)
    }, error = function(e) NA)
    beta2 <- tryCatch({
      f_beta(sample_mean_2, beta, a, b)
    }, error = function(e) NA)
    
    # initialize chains
    Inits_list <- list(
      list(mu1p = sample_mean_1,
           mu2p = sample_mean_2,
           tau1p = 20,
           tau2p = 20),
      list(mu1p = sample_mean_1 - 0.1,
           mu2p = sample_mean_2 - 0.1,
           tau1p = 21,
           tau2p = 21)
    )
    if (miss_prop <= 0.1){
      file_model <- 'model.txt'
      # jags data
      data_jags <- list(
        y = y,
        is_observed = is_observed,
        group_numeric = group_numeric,
        r = r,
        cp = cp,
        mu0 = mu_0, sigma_1 = sigma_1, alpha_p = alpha_p, beta_p = beta_p,
        alpha1 = alpha1, beta1 = beta1, alpha2 = alpha2, beta2 = beta2
      )
      
      if (anyNA(data_jags)) {
        mcmcresults_list[[i]] <- summary_row
        next
      }
    } else {
      file_model <- 'model_logit.txt'
      # jags data
      data_jags <- list(
        y = y,
        is_observed = is_observed,
        group_numeric = group_numeric,
        r = r,
        mu0 = mu_0, sigma_1 = sigma_1, alpha_p = alpha_p, beta_p = beta_p,
        alpha1 = alpha1, beta1 = beta1, alpha2 = alpha2, beta2 = beta2,
        gamma0 = gamma0, gamma1 = gamma1
      )

      if (anyNA(data_jags)) {
        mcmcresults_list[[i]] <- summary_row
        next
      }
    }
    
    # run the chains
    jagsModel <- jags.model(file = file_model, data = data_jags, inits = Inits_list, n.chains = 2, n.adapt = 500)
    print(jagsModel)
    update(jagsModel, n.iter = 500)
    codaSample <- coda.samples(jagsModel, variable.names = c('mu1p', 'mu2p'),
                               n.iter = 5000)
    # get HDI, median and the location of 0
    difference <- codaSample[[1]][,'mu1p'] - codaSample[[1]][,'mu2p']
    postSummary <- summarizePost(difference, compVal = 0, credMass=0.95, ROPE = c(-0.2, 0.2))
    summary_row <- data.frame(
      HDI_Low = postSummary["HDIlow"],
      HDI_High = postSummary["HDIhigh"],
      Median = postSummary["Median"],
      pLtCompVal = pmin(postSummary["PcntGtCompVal"] / 100, 
                        1 - postSummary["PcntGtCompVal"] / 100),
      pLtROPE = postSummary["PcntLtROPE"],
      pInROPE = postSummary["PcntInROPE"],
      pGtROPE = postSummary["PcntGtROPE"],
      DE = (postSummary["HDIlow"] > 0.2) | (postSummary["HDIhigh"] < -0.2)
    )
    mcmcresults_list[[i]] <- summary_row
  }
  mcmcresults_df <- do.call(rbind, mcmcresults_list)
  rownames(mcmcresults_df) <- rownames(subset_data)
  #diagMCMC(codaObject = codaSample, parName = 'mu1p' )
  
  return (mcmcresults_df)
}

# Wrapper function run by the user:
# INPUT ARGUEMENTS:
# 'values' is an array of logged protein*sample abundance values. Row names of 'values' should be used to specify protein identifiers.
# 'groups' is a factor containing groups IDs of each sample.
# 'comparisons' is a matrix showing which groups to comapre for each contrast, output by limma::makeContrasts.
# 'filter4NAs' is a bool denoting whether to only use proteins that possess missing values.
# OUTPUT:
# 'BayesResults' is a list of tables containing model results for each contrast specified in comparisons.
BayesMissingModel <-  function(values, groups, comparisons, filter4NAs = FALSE){
    
    zS <- zeroState(values)
    s1 <- sigma1(values)
    s2P <- sigma_p2params(values, groups)
    sjp2P <- sigma_jp2params(values, s2P$r)
    f_a <- f_alpha(zS['mu_0'], sjp2P$alpha, sjp2P$a, sjp2P$b)
    f_b <- f_beta(zS['mu_0'], sjp2P$beta, sjp2P$a, sjp2P$b)
    
    bayesResults <- list()
    for (contrast in colnames(comparisons)){
      cd <- setupContrasts(values,
                           comparisons,
                           contrast, 
                           groups, 
                           filter4NAs)
      
      bayesResults[[contrast]] <- runModel(values, cd, s2P, sjp2P, zS, s1)
    }
    
    return(bayesResults)
}

# QFeatures integration, runs model fitting for all specified comparisons and 
# adds them to existing QFeatures object.
# WRITE DOCUMENTATION!!
setGeneric("BayesianMissingAggregate", function(obj,...) standardGeneric("BayesianMissingAggregate"))

setMethod(
  "BayesianMissingAggregate", "QFeatures",
  function(obj, i, fcol, comparisons, name = 'BayesModels', filter4NAs = FALSE, contrastPrefix = ""){
    if (is.null(obj[[i]])) stop("QFeatures object does not contain assay ", i)
    
    log2all.df <- assay(obj[[i]])
    rownames(log2all.df) <- rowData(obj[[i]])[,fcol]

    results <- BayesMissingModel(log2all.df, colData(obj)$condition, comparisons, filter4NAs)
    
    for (contrast in colnames(comparisons)){
      rowData(obj[[i]])[[paste0(contrastPrefix, contrast)]] <- results[[contrast]]
    }
    
    return(obj)
  }
)
