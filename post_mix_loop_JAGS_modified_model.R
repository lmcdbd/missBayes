library(rjags)
library(runjags)
library(data.table)
library(bnlearn)
library(matrixStats)
library(coda)
library(ggplot2)
library(dplyr)
library(limma)
source("DBDA2E-utilities.R")
# 1. read data ----
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

# 2. input group factor, save data to a list ----
overall_distri <- log2all.df
group <- as.factor(metadata$Condition)
if (length(unique(table(group))) != 1) {
  stop("Error: The number of replicates is not the same across all groups.")
}

# 2.1 find mu0, gamma0 and gamma1 by plotting a histogram

overall_info <- hist(as.matrix(overall_distri), breaks = 100, plot = FALSE)
missing_prop <- sum(is.na(overall_distri)) / (ncol(overall_distri) * nrow(overall_distri)) 
overall_info$density <- overall_info$density * (1 - missing_prop)
mids <- overall_info$mids
dens <- overall_info$density
mu_0 <- mids[which.max(dens)]
# Empty vector to store estimated missing density per bin
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
df <- df[is.finite(df$p_obs) & df$p_obs > 0 & df$p_obs < 1, ]  # filter valid values

logit_obs <- glm(p_obs ~ intensity, data = df, family = binomial(link = "logit"))
summary(logit_obs)

# Extract coefficients
gamma0 <- coef(logit_obs)[1]
gamma1 <- coef(logit_obs)[2]


# 2.2 sigma1 is the standard deviation of the protein means
Protein_means <- rowMeans(overall_distri, na.rm = TRUE)
sigma_1 <- sd(Protein_means, na.rm = TRUE)

# 2.3 sigma_p^2 is from a inverse-gamma distribution with shape of alpha and rate of beta
r <- unname(table(group)[1])  # Number of replicates per group
n_groups <- nlevels(group)
n <- n_groups # n is the number of groups
between_group_vars_all <- apply(overall_distri, 1, function(row) {
  group_means <- sapply(1:n_groups, function(i) {
    group_indices <- (r * (i - 1) + 1):(r * i)
    mean(row[group_indices], na.rm = TRUE)
  })
  var(group_means, na.rm = TRUE)
})
mean_var <- mean(between_group_vars_all, na.rm = TRUE)
var_var <- var(between_group_vars_all, na.rm = TRUE)
alpha_p <- mean_var^2 / n*var_var + 2
beta_p <- mean_var * (alpha_p - 1)

# 2.4 sigma_jp^2 is from a inverse-gamma distribution with shape of alpha and rate of beta which are dependent of mu_jp
group_means_all <- c()
group_vars_all <- c()
for (i in 1:(ncol(overall_distri) %/% r)) {
  group_i <- overall_distri[, (r * (i - 1) + 1):(r * i)]
  group_i_mean <- rowMeans(group_i, na.rm = TRUE)
  group_i_var   <- rowVars(as.matrix(group_i), na.rm = TRUE)
  group_means_all <- c(group_means_all, group_i_mean)
  group_vars_all  <- c(group_vars_all, group_i_var)
}
#####calculate alpha and beta for each interval
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
# define the function f_alpha and f_beta that is mu is in a certain interval, return alpha and beta
bin_breaks <- a:(b + 1)  # From a to b+1 so that each bin is [x, x+1]
f_alpha <- function(mu_val) {
  idx <- findInterval(mu_val, bin_breaks, rightmost.closed = TRUE)
  if (!is.na(alpha[idx])) {
    return(alpha[idx])
  } else {
    stop(paste("No alpha value available for mu =", mu_val))
  }
}

f_beta <- function(mu_val) {
  idx <- findInterval(mu_val, bin_breaks, rightmost.closed = TRUE)
  if (!is.na(beta[idx])) {
    return(beta[idx])
  } else {
    stop(paste("No beta value available for mu =", mu_val))
  }
}


# 2.7 specify comparison groups, filter the dataset, only keep proteins with missing values
comparison <- "G70_20_10 - G70_10_20"
# Split the comparison string into two group names
groups_to_compare <- strsplit(comparison, " - ")[[1]]
# Check that both groups exist in the factor
if (!all(groups_to_compare %in% levels(group))) {
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
# only keep the rows with NA
subset_with_na <- subset_data[rowSums(is.na(subset_data)) > 0 & rowSums(is.na(subset_data)) < ncol(subset_data), ]
# get row minimums from overall distribution only for the rows present in filtered subset
rows_with_na <- rownames(subset_with_na)
overall_distri_with_na <- overall_distri[rows_with_na, , drop = FALSE]
row_mins <- apply(overall_distri_with_na, 1, min, na.rm = TRUE)


# 3. model ----
modelString <- "
  model {
    # Likelihood
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
"
writeLines(modelString, con = 'model.txt')

# 4. data prep ----
# loop through entire dataset
group_numeric <- ifelse(group_subset == groups_to_compare[1], 1,
                        ifelse(group_subset == groups_to_compare[2], 2, NA))
mcmcresults_list_rowmin <- list() 
for (i in 1: nrow(subset_with_na)){
  y <- unlist(subset_with_na[i, ])
  
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
    mcmcresults_list_rowmin[[i]] <- summary_row
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
    f_alpha(sample_mean_1)
  }, error = function(e) NA)
  beta1 <- tryCatch({
    f_beta(sample_mean_1)
  }, error = function(e) NA)
  alpha2 <- tryCatch({
    f_alpha(sample_mean_2)
  }, error = function(e) NA)
  beta2 <- tryCatch({
    f_beta(sample_mean_2)
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
  if (miss_prop <= 1){
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
      mcmcresults_list_rowmin[[i]] <- summary_row
      next
    }
  } else {
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
      mcmcresults_list_rowmin[[i]] <- summary_row
      next
    }
  }
  
  # run the chains
  jagsModel <- jags.model(file = file_model, data = data_jags, inits = Inits_list, n.chains = 2, n.adapt = 500)
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
  mcmcresults_list_rowmin[[i]] <- summary_row
}
mcmcresults_df_rowmin <- do.call(rbind, mcmcresults_list_rowmin)
rownames(mcmcresults_df_rowmin) <- rownames(subset_with_na)
#diagMCMC(codaObject = codaSample, parName = 'mu1p' )


# 5.1 compare limma CI width with HDI width ----
mcmcresults_df_rowmin$HDI_width <- mcmcresults_df_rowmin$HDI_High - mcmcresults_df_rowmin$HDI_Low
results_3_4_filtered$CI_width <- results_3_4_filtered$CI.R - results_3_4_filtered$CI.L
boxplot(
  list(CI = na.omit(results_3_4_filtered$CI_width),
       HDI = na.omit(mcmcresults_df_rowmin$HDI_width)),
  col = c("skyblue", "orange"),
  main = "CI vs HDI Width",
  ylab = "Interval Width"
)





# 5.2 compare with limma output ----
# calculate TPR and TNR, where TPR is yeast with -1.5 < median < -0.5, p <0.05  or ecoli with 0.5 < median < 1, p <0.05
TP <- 0
TN <- 0
FP <- 0
for (i in 1:nrow(mcmcresults_df_rowmin)) {
  row_name <- rownames(mcmcresults_df_rowmin)[i]
  p_val <- mcmcresults_df_rowmin$pLtCompVal[i]
  fc <- mcmcresults_df_rowmin$Median[i]
  DE <- mcmcresults_df_rowmin$DE[i]
  pLtROPE <- mcmcresults_df_rowmin$pLtROPE[i]
  pGtROPE <- mcmcresults_df_rowmin$pGtROPE[i]
  # Skip iteration if p_val or fc is NA
  if (is.na(p_val) || is.na(fc)) next
  
  if (grepl("_HUMAN", row_name) && DE == 'FALSE' && pLtROPE < 95 && pGtROPE < 95) {
    TN <- TN + 1
  } else if (#(grepl("_ECOLI", row_name) && DE == 'TRUE' && fc < -0.5 && fc > -1.5)
             (grepl("_ECOLI", row_name)  && fc < -0.5 && fc > -1.5 && pLtROPE > 95)
             |(grepl("_ECOLI", row_name) && fc < -0.5 && fc > -1.5 && pGtROPE > 95)) {
    TP <- TP + 1
  } else if (#(grepl("_YEAST", row_name) && DE == 'TRUE' && fc > 0.5 && fc < 1.5)
             (grepl("_YEAST", row_name) && fc > 0.5 && fc < 1.5 && pLtROPE > 95)
             |(grepl("_YEAST", row_name) && fc > 0.5 && fc < 1.5 && pGtROPE > 95)) {
    TP <- TP + 1}
  #} else if (grepl("_HUMAN", row_name) && DE == 'TRUE') {
    FP <- FP + 1
  #}
}
FP <- sum(grepl("_HUMAN", rownames(mcmcresults_df_rowmin))) - TN
cat("True Positives:", TP, "\n")
cat("True Negatives:", TN, "\n")
cat("False Positives:", FP, "\n")

# limma result

metadata <- fread("P:/Trost-group/Mengchun/8) MIP/CQE/HYE_3spp/metadata.txt", header = TRUE)

design <- model.matrix(~ 0 + Condition, data = metadata)
colnames(design) <- gsub("Condition", "", colnames(design))

fit <- lmFit(log2all.df, design)
fit <- eBayes(fit)
contrast_matrix <- makeContrasts(
  HYE_1_2 = G35_25_40 - G35_40_25,
  HYE_3_4 = G70_10_20 - G70_20_10,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
results_3_4 <- topTable(fit2, coef = "HYE_3_4", adjust.method = "BH", number = Inf, confint = TRUE)
results_3_4_filtered <- results_3_4[rownames(results_3_4) %in% rownames(overall_distri_with_na), ]

TP_limma <- 0
TN_limma <- 0
FP_limma <- 0
for (i in 1:nrow(results_3_4_filtered)) {
  row_name <- rownames(results_3_4_filtered)[i]
  p_val <- results_3_4_filtered$adj.P.Val[i]
  fc <- results_3_4_filtered$logFC[i]
  
  # Skip iteration if p_val or fc is NA
  if (is.na(p_val) || is.na(fc)) next
  
  if (grepl("_HUMAN", row_name) && p_val > 0.05) {
    TN_limma <- TN_limma + 1
  } else if (grepl("_YEAST", row_name) && p_val <= 0.05 && fc < -0.5 && fc > -1.5) {
    TP_limma <- TP_limma + 1
  } else if (grepl("_ECOLI", row_name) && p_val <= 0.05 && fc > 0.5 && fc < 1.5) {
    TP_limma <- TP_limma + 1
  } else if (grepl("_HUMAN", row_name) && p_val <= 0.05) {
    FP_limma <- FP_limma + 1
  }
}
cat("True Positives:", TP_limma, "\n")
cat("True Negatives:", TN_limma, "\n")
cat("False Positives:", FP_limma, "\n")

# perseus
impute_normal <- function(object, width=0.3, downshift=1.8, seed=100) {
  
  if (!is.matrix(object)) object <- as.matrix(object)
  mx <- max(object, na.rm=TRUE)
  mn <- min(object, na.rm=TRUE)
  if (mx - mn > 20) warning("Please make sure the values are log-transformed")
  
  set.seed(seed)
  object <- apply(object, 2, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm=TRUE)
    temp_mean <- mean(temp, na.rm=TRUE)
    shrinked_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
    temp
  })
  return(object)
}
log2all.df.perseus <- impute_normal(log2all.df, width = 0.3, downshift = 1.8, seed = 100)

fit <- lmFit(log2all.df.perseus, design)
fit <- eBayes(fit)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results_3_4_perseus <- topTable(fit2, coef = "HYE_3_4", adjust.method = "BH", number = Inf)
results_3_4_filtered_perseus <- results_3_4_perseus[rownames(results_3_4_perseus) %in% rownames(overall_distri_with_na), ]

TP_limma <- 0
TN_limma <- 0
FP_limma <- 0
for (i in 1:nrow(results_3_4_filtered_perseus)) {
  row_name <- rownames(results_3_4_filtered_perseus)[i]
  p_val <- results_3_4_filtered_perseus$adj.P.Val[i]
  fc <- results_3_4_filtered_perseus$logFC[i]
  
  # Skip iteration if p_val or fc is NA
  if (is.na(p_val) || is.na(fc)) next
  
  if (grepl("_HUMAN", row_name) && p_val > 0.05) {
    TN_limma <- TN_limma + 1
  } else if (grepl("_YEAST", row_name) && p_val <= 0.05 && fc < -0.5 && fc > -1.5) {
    TP_limma <- TP_limma + 1
  } else if (grepl("_ECOLI", row_name) && p_val <= 0.05 && fc > 0.5 && fc < 1.5) {
    TP_limma <- TP_limma + 1
  } else if (grepl("_HUMAN", row_name) && p_val <= 0.05) {
    FP_limma <- FP_limma + 1
  }
}
cat("True Positives:", TP_limma, "\n")
cat("True Negatives:", TN_limma, "\n")
cat("False Positives:", FP_limma, "\n")

# LOD
log2all.df.LOD <- log2all.df  # Make a copy

# Replace NAs with row-wise minimums
log2all.df.LOD <- t(apply(log2all.df.LOD, 1, function(row) {
  if (all(is.na(row))) return(row)  # Leave all-NA rows as-is
  row[is.na(row)] <- min(row, na.rm = TRUE)
  return(row)
}))

fit <- lmFit(log2all.df.LOD, design)
fit <- eBayes(fit)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results_3_4_LOD <- topTable(fit2, coef = "HYE_3_4", adjust.method = "BH", number = Inf)
results_3_4_filtered_LOD <- results_3_4_LOD[rownames(results_3_4_LOD) %in% rownames(overall_distri_with_na), ]

TP_limma <- 0
TN_limma <- 0
FP_limma <- 0
for (i in 1:nrow(results_3_4_filtered_LOD)) {
  row_name <- rownames(results_3_4_filtered_LOD)[i]
  p_val <- results_3_4_filtered_LOD$adj.P.Val[i]
  fc <- results_3_4_filtered_LOD$logFC[i]
  
  # Skip iteration if p_val or fc is NA
  if (is.na(p_val) || is.na(fc)) next
  
  if (grepl("_HUMAN", row_name) && p_val > 0.05) {
    TN_limma <- TN_limma + 1
  } else if (grepl("_YEAST", row_name) && p_val <= 0.05 && fc < -0.5 && fc > -1.5) {
    TP_limma <- TP_limma + 1
  } else if (grepl("_ECOLI", row_name) && p_val <= 0.05 && fc > 0.5 && fc < 1.5) {
    TP_limma <- TP_limma + 1
  } else if (grepl("_HUMAN", row_name) && p_val <= 0.05) {
    FP_limma <- FP_limma + 1
  }
}
cat("True Positives:", TP_limma, "\n")
cat("True Negatives:", TN_limma, "\n")
cat("False Positives:", FP_limma, "\n")

#kNN
# Convert to matrix
log2all.matrix <- as.matrix(log2all.df)

# Perform kNN imputation
imputed_data <- impute.knn(log2all.matrix, k = 5)$data

# Convert back to data frame
log2all.df.kNN <- as.data.frame(imputed_data)

fit <- lmFit(log2all.df.kNN, design)
fit <- eBayes(fit)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results_3_4_kNN <- topTable(fit2, coef = "HYE_3_4", adjust.method = "BH", number = Inf)
results_3_4_filtered_kNN <- results_3_4_kNN[rownames(results_3_4_kNN) %in% rownames(overall_distri_with_na), ]

TP_limma <- 0
TN_limma <- 0
FP_limma <- 0
for (i in 1:nrow(results_3_4_filtered_kNN)) {
  row_name <- rownames(results_3_4_filtered_kNN)[i]
  p_val <- results_3_4_filtered_kNN$adj.P.Val[i]
  fc <- results_3_4_filtered_kNN$logFC[i]
  
  # Skip iteration if p_val or fc is NA
  if (is.na(p_val) || is.na(fc)) next
  
  if (grepl("_HUMAN", row_name) && p_val > 0.05) {
    TN_limma <- TN_limma + 1
  } else if (grepl("_YEAST", row_name) && p_val <= 0.05 && fc < -0.5 && fc > -1.5) {
    TP_limma <- TP_limma + 1
  } else if (grepl("_ECOLI", row_name) && p_val <= 0.05 && fc > 0.5 && fc < 1.5) {
    TP_limma <- TP_limma + 1
  } else if (grepl("_HUMAN", row_name) && p_val <= 0.05) {
    FP_limma <- FP_limma + 1
  }
}
cat("True Positives:", TP_limma, "\n")
cat("True Negatives:", TN_limma, "\n")
cat("False Positives:", FP_limma, "\n")

# 6. volcano plot limma vs mcmc ----
# Add a new column for category based on rownames
mcmcresults_df_rowmin$Group <- ifelse(grepl("_HUMAN", rownames(mcmcresults_df_rowmin)), "HUMAN",
                                      ifelse(grepl("_ECOLI", rownames(mcmcresults_df_rowmin)), "ECOLI",
                                             ifelse(grepl("_YEAST", rownames(mcmcresults_df_rowmin)), "YEAST", "OTHER")))

# Plot
ggplot(mcmcresults_df_rowmin, aes(x = Median, y = -log10(pGtCompVal), color = Group)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("HUMAN" = "skyblue", "ECOLI" = "red", "YEAST" = "yellow")) +
  labs(title = "mcmc Results",
       x = "Posterior Median",
       y = "pValue") +
  ylim(0,5) +
  theme_minimal()


# 7. performance on entirely missing groups ----
first_five_na <- apply(overall_distri_filtered[, 1:5], 1, function(row) all(is.na(row)))
last_five_na <- apply(overall_distri_filtered[, 6:10], 1, function(row) all(is.na(row)))
na_rows <- which(first_five_na | last_five_na)
mcmcresults_df_rowmin_allmissing <- mcmcresults_df_rowmin[na_rows, ]
results_3_4_filtered <- results_3_4_perseus[rownames(results_3_4_perseus) %in% rownames(mcmcresults_df_rowmin_allmissing), ]


# 8. ampute a complete dataset (use ritika's data, original FC as true FC) ----
log2all.df_complete <- log2all.df[complete.cases(log2all.df), ]
hist(as.matrix(final_df), col = 'skyblue', breaks = 100)

# Flatten to long format with row index
library(tidyr)
library(dplyr)

# Add row numbers for tracking
log2all.df_complete$row_id <- seq_len(nrow(log2all.df_complete))

# Convert to long format
long_df <- log2all.df_complete %>%
  pivot_longer(-row_id, names_to = "sample", values_to = "value") %>%
  filter(!is.na(value))

# Split values into <16 and >16
left_values <- long_df %>% filter(value < 16)
right_values <- long_df %>% filter(value > 16)

# Build histograms (mirrored around 16)
breaks <- seq(floor(min(long_df$value)), ceiling(max(long_df$value)), by = 0.5)
hist_left <- hist(left_values$value, breaks = breaks, plot = FALSE)
hist_right <- hist(right_values$value, breaks = breaks, plot = FALSE)

# Mirror the left histogram to define max allowed for the right
target_counts_right <- rev(hist_left$counts[hist_left$mids < 16])

# Compare and get indices of "extra" values on the right
excess_indices <- c()
for (i in seq_along(target_counts_right)) {
  bin_mid <- hist_right$mids[hist_right$mids > 16][i]
  current_bin_vals <- right_values %>%
    filter(value >= (bin_mid - 0.25), value < (bin_mid + 0.25))
  
  excess_n <- nrow(current_bin_vals) - target_counts_right[i]
  if (excess_n > 0) {
    # Mark rows to remove (as many as needed)
    excess_rows <- head(current_bin_vals$row_id, excess_n)
    excess_indices <- c(excess_indices, excess_rows)
  }
}

# Remove entire rows containing excess right-side values
final_df <- log2all.df_complete[!log2all.df_complete$row_id %in% excess_indices, ]

# Remove helper column
final_df$row_id <- NULL

# define a logistic curve

L <- 1      # Max value
k <- 3      # Steepness
x0 <- 13.5     # Midpoint

# Plot the logistic curve
curve(L / (1 + exp(-k * (x - x0))),
      from = 10, to = 17,
      col = "blue", lwd = 2,
      ylab = "f(x)", xlab = "x",
      main = "Logistic Curve")

# ampute the final_df with the probability of the logistic regression
# Set seed for reproducibility
set.seed(123)

# Logistic function
logistic <- function(x, L = 1, k = 1.5, x0 = 14) {
  L / (1 + exp(-k * (x - x0)))
}

# Copy final_df to ampute
amputed_df <- final_df

# Apply missingness row-wise
amputed_df[] <- lapply(final_df, function(col) {
  sapply(col, function(x) {
    if (is.na(x)) return(NA)
    prob_missing <- 1 - logistic(x)
    if (runif(1) < prob_missing) NA else x
  })
})

hist(unlist(amputed_df), breaks = 50, main = "After Amputation", col = "skyblue", xlab = "x")





# 9. RMSE calculation ----
mcmcresults_df_rowmin$trueVal <- ifelse(grepl("HUMAN", rownames(mcmcresults_df_rowmin)), 0,
                                        ifelse(grepl("ECOLI", rownames(mcmcresults_df_rowmin)), 1,
                                               ifelse(grepl("YEAST", rownames(mcmcresults_df_rowmin)), -1, NA)))
rmse_mcmc <- sqrt(mean((mcmcresults_df_rowmin$trueVal - mcmcresults_df_rowmin$Median)^2, na.rm = TRUE))
#NRMSE_mcmc <- rmse_mcmc / sd(original_data, na.rm = TRUE)
results_3_4_filtered$trueVal <- ifelse(grepl("HUMAN", rownames(results_3_4_filtered)), 0,
                                        ifelse(grepl("ECOLI", rownames(results_3_4_filtered)), 1,
                                               ifelse(grepl("YEAST", rownames(results_3_4_filtered)), -1, NA)))
rmse_limma <- sqrt(mean((results_3_4_filtered$trueVal - results_3_4_filtered$logFC)^2, na.rm = TRUE))
results_3_4_filtered_kNN$trueVal <- ifelse(grepl("HUMAN", rownames(results_3_4_filtered_kNN)), 0,
                                       ifelse(grepl("ECOLI", rownames(results_3_4_filtered_kNN)), 1,
                                              ifelse(grepl("YEAST", rownames(results_3_4_filtered_kNN)), -1, NA)))
rmse_kNN <- sqrt(mean((results_3_4_filtered_kNN$trueVal - results_3_4_filtered_kNN$logFC)^2, na.rm = TRUE))


# 10. venn diagram ----
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = 308 + 10567,       # total proteins in condition 1
  area2 = 609 + 10567,       # total proteins in condition 2
  cross.area = 10567,        # shared proteins
  category = c("Condition 2", "Condition 1"),
  fill = c("lightblue", "lightgreen"),
  lty = "blank",             # no border lines
  cex = 2,                   # text size for counts
  cat.cex = 1.5,             # text size for category labels
  cat.pos = c(0, 0),     # label positions
  scaled = FALSE)

