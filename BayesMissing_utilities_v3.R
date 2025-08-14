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
library(parallel)
source("P:/Trost-group/Mengchun/8) MIP/Rscript/DBDA2E-utilities.R")  # this file is from Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition, which provides some very useful functions


# 2.1 calculate mu0, gamma0 and gamma1 -------
zeroState <- function(overall_distri){
  
  overall_info <- hist(as.matrix(overall_distri), breaks = 100, plot = FALSE)
  missing_prop <- sum(is.na(overall_distri)) / (ncol(overall_distri) * nrow(overall_distri)) 
  overall_info$density <- overall_info$density * (1 - missing_prop)
  mids <- overall_info$mids
  dens <- overall_info$density
  mu_0 <- mids[which.max(dens)]
  # empty vector to store estimated missing density per bin
  missing_counts <- rep(0, length(overall_info$counts))
  for (i in seq_along(mids)) {
    if (mids[i] < mu_0) {
      mirror_i <- which.min(abs(mids - (2 * mu_0 - mids[i])))
      expected <- overall_info$counts[mirror_i]
      observed <- overall_info$counts[i]
      missing_counts[i] <- max(0, expected - observed)
    }
  }
  
  # estimate observed probability per bin
  observed_count <- overall_info$counts
  total_count <- observed_count + missing_counts
  obs_prob <- observed_count / total_count
  
  # fit logistic regression
  df <- data.frame(
    intensity = mids,
    observed = observed_count,
    missing = missing_counts,
    total = total_count,
    p_obs = obs_prob
  )
  df <- df[is.finite(df$p_obs) & df$p_obs > 0 & df$p_obs < 1 & df$total > 0, ]  
  logit_obs <- glm(cbind(observed, missing) ~ intensity, data = df, family = binomial(link = "logit"))
  
  # extract coefficients
  gamma0 <- coef(logit_obs)[1]
  gamma1 <- coef(logit_obs)[2]
  zS <- c(mu_0, gamma0, gamma1) 
  names(zS) <- c('mu_0', 'gamma0', 'gamma1')
  
  return (zS)
}

# 2.2 sigma1 is the standard deviation of the protein means -------
sigma1 <- function(overall_distri){
  Protein_means <- rowMeans(overall_distri, na.rm = TRUE)
  sigma_1 <- sd(Protein_means, na.rm = TRUE)
  
  return (sigma_1)
}

# 2.3 sigma_p^2 is from a inverse-gamma distribution with shape of alpha and rate of beta -------
sigma_p2params <- function(overall_distri, group){
  
  r <- unname(table(group)[1])  # Number of replicates per group
  n <- nlevels(group) # Number of groups
  
  group_means <- t(apply(overall_distri, 1, function(row) {
    sapply(1:n, function(i) {
      group_indices <- (r * (i - 1) + 1):(r * i)
      mean(row[group_indices], na.rm = TRUE)
    })
  }))
  between_group_vars_all <- apply(group_means, 1, var, na.rm = TRUE)
  
  mean_var <- mean(between_group_vars_all, na.rm = TRUE)
  var_var <- var(between_group_vars_all, na.rm = TRUE)
  alpha_p <- mean_var^2 / n*var_var + 2
  beta_p <- mean_var * (alpha_p - 1)
  
  sigma_p2params <- list(alpha_p, beta_p, as.numeric(r), group_means) ########
  names(sigma_p2params) <- c('alpha_p', 'beta_p', 'r', 'group_means')
  
  return (sigma_p2params)
}

# 2.4 sigma_jp^2 is from a inverse-gamma distribution with shape of alpha and rate of beta which are dependent of mu_jp -------
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
  params <- c(list(alpha), list(beta), a, b, list(group_vars_all))
  names(params) <- c('alpha', 'beta', 'a', 'b', 'group_vars_all')
  return (params)
}

# 2.5 define the function f_alpha and f_beta that in a certain interval, return alpha and beta for each mu in a certain interval-------
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

# 2.6 specify all comparisons, filter the dataset, only keep proteins with missing values -------
setupContrasts <- function(overall_distri, comparisons, contrast, group, filter4NAs){
  # Split the comparison string into two group names
  groups_to_compare <- names(comparisons[,contrast])[comparisons[,contrast] != 0]
  
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

# 2.7 generate initial list -------
generate_overdispersed_inits <- function(sample_mean_1, sample_mean_2, sample_var_1, sample_var_2, 
                                         n_chains, spread_factor = 3) {
  inits_list <- vector("list", n_chains)
  
  for (i in seq_len(n_chains)) {
    # Alternate positive and negative spread
    sign1 <- ifelse(i %% 2 == 0, 1, -1)
    sign2 <- ifelse(i %% 2 == 0, -1, 1)
    
    inits_list[[i]] <- list(
      mu1p = sample_mean_1 + sign1 * rnorm(1, sd = sqrt(sample_var_1) * spread_factor),
      mu2p = sample_mean_2 + sign2 * rnorm(1, sd = sqrt(sample_var_2) * spread_factor),
      tau1p = 1 / (sample_var_1 * runif(1, 0.5, 1.5)),  # wider range
      tau2p = 1 / (sample_var_2 * runif(1, 0.5, 1.5))
    )
  }
  
  return(inits_list)
}

# 2.8 default values ------

default_val <- function(overall_distri, group_vars_all){
  default_mean <- quantile(unlist(overall_distri), probs = 0.05, na.rm = TRUE) 
  default_var <- median(unlist(group_vars_all), na.rm = TRUE)
  default <- c(default_mean, default_var)
  names(default) <- c('default_mean', 'default_var')
  
  return(default)
}
# 2.9 process_row ------
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
  
  sample_var_1 <- var(y1, na.rm = TRUE)
  sample_var_1 <- if (is.na(sample_var_1)) fixed_data$default_var else sample_var_1
  sample_var_2 <- var(y2, na.rm = TRUE)
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
  if (miss_prop <= threshold) {
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
    jagsModel <- jags.model(
      file = textConnection(model_to_use),
      data = data_jags,
      inits = inits_list,
      n.chains = fixed_data$n_chains,
      n.adapt = fixed_data$n_adpt,
      quiet = TRUE
    )
    update(jagsModel, n.iter = fixed_data$n_burnin, progress.bar = "none")
    codaSample <- coda.samples(
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
      SD <- sd(difference)
      EffChnLngth <- effectiveSize(difference)
      rhat_values <- gelman.diag(codaSample)$psrf[, 'Point est.']
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


# 3. runModel ----
runModel <- function(contrast_data, sigma_p2p, sigma_jp2p, zS, sigma_1, default_data, 
                     n.chains, burn.in, n.iter, n.adapt, parallel, mcmcDiag, threshold, ROPE){
  e1 <- environment()
  fixed_data <- list(
    subset_data = contrast_data$subset_data,
    group_subset = contrast_data$group_subset,
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
  
  
  model_txt <- readLines('model.txt')
  model_logit_txt <- readLines('model_logit.txt')
  
  if (parallel){
    n_cores <- detectCores() - 1
    
    if (.Platform$OS.type == "unix") {
      # Unix-like systems (Linux/macOS)
      mcmcresults_list_all <- mclapply(
        X = 1:nrow(fixed_data$subset_data),
        FUN = function(i) {
          process_row(i, fixed_data, model_txt, model_logit_txt)
        },
        mc.cores = n_cores
      )
    } else {
      # Windows systems
      cl <- makeCluster(n_cores)

      # Export ALL required objects and functions
      clusterExport(cl, varlist = ls(envir = globalenv()), envir = globalenv())
      
      clusterEvalQ(cl, {
        library(rjags)
        library(coda)
      })
      
      mcmcresults_list_all <- parLapply(
        cl = cl,
        X = 1:nrow(fixed_data$subset_data),
        fun = function(i) {
          process_row(i, fixed_data, model_txt, model_logit_txt, mcmcDiag, threshold, ROPE)
        }
      )
      stopCluster(cl)
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

# 4. BayesMissingModel -----
# Wrapper function run by the user: 
# INPUT ARGUEMENTS:
# 'values' is an array of logged protein*sample abundance values. Row names of 'values' should be used to specify protein identifiers.
# 'groups' is a factor containing groups IDs of each sample.
# 'comparisons' is a matrix showing which groups to compare for each contrast, output by limma::makeContrasts.
# 'filter4NAs' is a bool denoting whether to only use proteins that possess missing values.
# 'threshold' is a numeric representing the missing proportion of the groups being compared under which the logistic model will be
# -applied, otherwise the truncated normal model will be applied.
# 'parallel' is a bool denoting whether to do parallel computing.
# 'ROPE' is a vector of two elements, specifying the region of practical equivalence with in which the log fold change is considered
# -equivalent to zero.
# 'n.adapt', 'burn.in', 'n.iter' and 'n.chains' are integers denoting number of adaptation steps, burn-in steps, iteration taken by each
#  -chain and number of chains.
# 'mcmcDiag' is a bool denoting whether to show MCMC diagnostics in the output.
# OUTPUT:
# 'BayesResults' is a list of tables containing model results for each contrast specified in comparisons.
BayesMissingModel <-  function(values, groups, comparisons, filter4NAs = FALSE, threshold = 0, parallel = TRUE, ROPE = c(-0.2,0.2), 
                               n.adapt = 500,burn.in = 500, n.iter = 10000, n.chains = 2, mcmcDiag = FALSE){
  
  zS <- zeroState(values)
  s1 <- sigma1(values)
  s2P <- sigma_p2params(values, groups)
  sjp2P <- sigma_jp2params(values, s2P$r)
  d <- default_val(values, sjp2P$group_vars_all)
  
  bayesResults <- list()
  for (contrast in colnames(comparisons)){
    cd <- setupContrasts(values,
                         comparisons,
                         contrast, 
                         groups, 
                         filter4NAs)
    
    bayesResults[[contrast]] <- runModel(cd, s2P, sjp2P, zS, s1, d, n.chains, burn.in, n.iter, n.adapt, parallel, 
                                         mcmcDiag, threshold, ROPE)
  }
  
  return(bayesResults)
}

# 5. TEST ------
# import data
pg<-fread("P:/Trost-group/Mengchun/8) MIP/CQE/HYE_3spp/report.pg_matrix.tsv")
pr<-fread("P:/Trost-group/Mengchun/8) MIP/CQE/HYE_3spp/report.pr_matrix.tsv")
t <- table(unique(pr[,c('Protein.Group','Stripped.Sequence')])$Protein.Group)
pg$Peptide.Count <- t[match(pg$Protein.Group,names(t))]
metadata <- fread("P:/Trost-group/Mengchun/8) MIP/CQE/HYE_3spp/metadata.txt", header = TRUE)
#Filter by protein groups with <2 peptides
df.prot <- filter(pg, Peptide.Count > 1)
df.prot2 <- df.prot[!grepl("Cont_",df.prot$Protein.Group),]
df.prot.all<- df.prot2[, 5:29]
protein.matrix_all <- log2(as.matrix(df.prot.all))
log2all.df <- as.data.frame(protein.matrix_all)
rownames(log2all.df) <- df.prot2$Protein.Names

overall_distri <- log2all.df
group <- as.factor(metadata$Condition)
x <- c('G70_10_20 - G70_20_10', 'G35_25_40 - G35_40_25')
comparison <- limma::makeContrasts(contrasts = x, levels = levels(group))
start.time <- Sys.time()
output <- BayesMissingModel(overall_distri, group, comparison, filter4NAs = TRUE, parallel = TRUE)
finish.time <- Sys.time()

# QFeatures ------
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

# result
TP_mcmc <- 0
TN_mcmc <- 0
FP_mcmc <- 0
FN_mcmc <- 0
for (i in 1:nrow(output[[1]])) {
  row_name <- rownames(output[[1]])[i]
  p_val <- output[[1]]$pLtCompVal[i]
  fc <- output[[1]]$Median[i]
  pLtROPE <- output[[1]]$pLtROPE[i]
  pGtROPE <- output[[1]]$pGtROPE[i]
  
  # Skip iteration if p_val or fc is NA
  if (is.na(p_val) || is.na(fc)) next
  
  if (grepl("_HUMAN", row_name) && pLtROPE < 95 && pGtROPE < 95) {
    TN_mcmc <- TN_mcmc + 1
  } else if (grepl("_ECOLI", row_name) && pGtROPE > 95 && fc > 0.5 && fc < 1.5  ) {
    TP_mcmc <- TP_mcmc + 1
  } else if (grepl("_YEAST", row_name) && pLtROPE > 95 && fc < -0.5 && fc > -1.5) {
    TP_mcmc <- TP_mcmc + 1
  } else if ((grepl("_HUMAN", row_name) && pLtROPE >= 95) | (grepl("_HUMAN", row_name) && pGtROPE >= 95)) {
    FP_mcmc <- FP_mcmc + 1
  } else if (grepl("_ECOLI", row_name)|grepl("_YEAST", row_name)){
    FN_mcmc <- FN_mcmc + 1
  }
}
cat("True Positives:", TP_mcmc, "\n") # 905
cat("True Negatives:", TN_mcmc, "\n") # 579
cat("False Positives:", FP_mcmc, "\n") # 25
cat("False Negatives:", FN_mcmc, "\n") # 882
