
#' Function for computing limits of HDI's: Computes highest density interval from a sample of representative values, estimated as shortest credible interval.
#' This function was written by John K. Kruschke for his book
#' "Doing Bayesian Data Analysis" (2nd edition).
#'
#' @param sampleVec a vector of representative values from a probability distribution.
#' @param credMass is a scalar between 0 and 1, indicating the mass within the credible interval that is to be estimated.
#'
#' @return a vector containing the limits of the HDI
#' @keywords internal
#' @references
#'   Kruschke, J. K. (2015). Doing Bayesian Data Analysis:
#'   A Tutorial with R, JAGS, and Stan (2nd ed.). Academic Press.
#'

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}


#' Plot Posterior Distribution
#'
#' @description Creates a plot of posterior distribution with HDI and other statistics.
#'
#' @details This function is modified from the \code{plotPost} function written by John K. Kruschke for his book
#' "Doing Bayesian Data Analysis" (2nd edition).
#' Modifications include: Users can directly input the data, protein ID, and other information to get the posterior distribution.
#'
#'
#' @param values A numeric matrix or data.frame of log2 intensities. Rows correspond to proteins
#' (with protein IDs as rownames), and columns correspond to samples. Each row must
#' contain at least one non-NA value.
#' @param id A character string specifying the protein ID to plot.
#' @return A posterior distribution plot.
#' @author John K. Kruschke (original `plotPost` function), Mengchun Li (port and integration with `missBayes`)
#' @references
#' Kruschke, J. K. (2015). Doing Bayesian Data Analysis:
#' A Tutorial with R, JAGS, and Stan (2nd ed.). Academic Press.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_data' is a matrix and 'groups' is a factor
#' plotPost(
#'   values = my_data,
#'   id = "Protein_ID",
#'   groups = groups,
#'   contrast = "Treatment - Control",
#'   compVal = 0,
#'   ROPE = c(-0.1, 0.1)
#' )
#' }


plotPost <- function( values, id , groups, contrast, threshold = 0, cenTend=c("mode","median","mean")[1] ,
                      compVal=0, ROPE=c(-0.2,0.2),
                      n.adapt=500,burn.in = 500, n.iter = 5000, n.chains = 2,
                      credMass=0.95, HDItextPlace=0.7,
                      xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL ,
                      main=NULL , cex=NULL , cex.lab=NULL ,
                      col=NULL , border=NULL , showCurve=FALSE , breaks=NULL ,
                      ... ) {


  # calculate parameters
  zS <- zeroState(values)
  s1 <- sigma1(values)
  s2P <- sigma_p2params(values, groups)
  sjp2P <- sigma_jp2params(values, s2P$r)
  d <- default_val(values, sjp2P$group_vars_all)

  groups_to_compare <- strsplit(contrast, " - ")[[1]]
  # Check that both groups exist in the factor
  if (!all(groups_to_compare %in% levels(groups))) {
    stop("One or both groups in the comparison do not exist in 'group' levels.")
  }
  # Get indices for each group
  group1_indices <- which(groups == groups_to_compare[1])
  group2_indices <- which(groups == groups_to_compare[2])
  # Combine indices if needed
  selected_indices <- c(group1_indices, group2_indices)
  # Subset the data
  subset_data <- values[, selected_indices]
  group_subset <- groups[selected_indices]
  subset_data <- as.matrix(subset_data)
  row_mins <- apply(values, 1, min, na.rm = TRUE)





  # fixed data
  fixed_data <- list(
    subset_data = subset_data,
    group_subset = group_subset,
    groups_to_compare = groups_to_compare,
    row_mins = row_mins,
    group_numeric = ifelse(group_subset == groups_to_compare[1], 1,
                           ifelse(group_subset == groups_to_compare[2], 2, NA)),

    sigma_1 = s1,

    alpha_p = s2P$alpha_p,
    beta_p = s2P$beta_p,
    r = s2P$r,

    alpha = sjp2P$alpha,
    beta = sjp2P$beta,
    a = sjp2P$a,
    b = sjp2P$b,

    mu_0 = zS['mu_0'],
    gamma0 = zS['gamma0'],
    gamma1 = zS['gamma1'],

    n_chains = n.chains,
    n_burnin = burn.in,
    n_iter = n.iter,
    n_adpt = n.adapt,

    group1_idx = which(group_subset == groups_to_compare[1]),
    group2_idx = which(group_subset == groups_to_compare[2]),

    default_mean = d['default_mean'],
    default_var = d['default_var']
  )

  model_txt_path <- system.file("model.txt", package = "missBayes")
  model_logit_txt_path <- system.file("model_logit.txt", package = "missBayes")

  if(model_txt_path == "" || model_logit_txt_path == "") {
    model_txt <- readLines("inst/model.txt", warn = FALSE)
    model_logit_txt <- readLines("inst/model_logit.txt", warn = FALSE)
  } else {
    model_txt <- readLines(model_txt_path, warn = FALSE)
    model_logit_txt <- readLines(model_logit_txt_path, warn = FALSE)
  }

  # sampling
  y <- unlist(fixed_data$subset_data[id, ])
  # Skip if all NA
  if (all(is.na(y))) {stop("All values are missing!")}

  # Calculate missing proportion
  miss_prop <- mean(is.na(y))
  cp <- fixed_data$row_mins[id]
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
  if (any(is.na(c(alpha1, beta1, alpha2, beta2)))) {stop("Hyperparameter not avaliable.")}
  inits_list <- generate_overdispersed_inits(sample_mean_1, sample_mean_2,
                                             sample_var_1, sample_var_2,
                                             n_chains = fixed_data$n_chains)
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
      cp = fixed_data$row_mins[id],
      mu0 = fixed_data$mu_0,
      sigma_1 = fixed_data$sigma_1,
      alpha_p = fixed_data$alpha_p,
      beta_p = fixed_data$beta_p,
      alpha1 = alpha1, beta1 = beta1,
      alpha2 = alpha2, beta2 = beta2
    )
    model_to_use <- model_txt
  }

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
  })
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Param. Val."
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( compVal , ROPE , difference ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="skyblue"
  if ( is.null(border) ) border="white"

  # convert coda object to matrix:
  if ( inherits(difference, "mcmc.list") ) {
    difference = as.matrix(difference)
  }

  summaryColNames = c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh",
                      "compVal","pGtCompVal",
                      "ROPElow","ROPEhigh","pLtROPE","pInROPE","pGtROPE")
  postSummary = matrix( NA , nrow=1 , ncol=length(summaryColNames) ,
                        dimnames=list( c( xlab ) , summaryColNames ) )

  # require(coda) # for effectiveSize function
  postSummary[,"ESS"] = effectiveSize(difference)

  postSummary[,"mean"] = mean(difference)
  postSummary[,"median"] = median(difference)
  mcmcDensity = density(difference)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]

  HDI = HDIofMCMC( difference , credMass )
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]

  # Plot histogram.
  cvCol = "darkgreen"
  ropeCol = "darkred"
  if ( is.null(breaks) ) {
    if ( max(difference) > min(difference) ) {
      breaks = c( seq( from=min(difference) , to=max(difference) ,
                       by=(HDI[2]-HDI[1])/18 ) , max(difference) )
    } else {
      breaks=c(min(difference)-1.0E-6,max(difference)+1.0E-6)
      border="skyblue"
    }
  }
  if ( !showCurve ) {
    par(xpd=NA)
    histinfo = hist( difference , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , ... )
  }
  if ( showCurve ) {
    par(xpd=NA)
    histinfo = hist( difference , plot=F )
    densCurve = density( difference , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  cenTendHt = 0.9*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.55*max(histinfo$density)
  # Display central tendency:
  mn = mean(difference)
  med = median(difference)
  mcmcDensity = density(difference)
  mo = mcmcDensity$x[which.max(mcmcDensity$y)]
  if ( cenTend=="mode" ){
    text( mo , cenTendHt ,
          bquote(mode==.(signif(mo,3))) , adj=c(.5,0) , cex=cex )
  }
  if ( cenTend=="median" ){
    text( med , cenTendHt ,
          bquote(median==.(signif(med,3))) , adj=c(.5,0) , cex=cex , col=cvCol )
  }
  if ( cenTend=="mean" ){
    text( mn , cenTendHt ,
          bquote(mean==.(signif(mn,3))) , adj=c(.5,0) , cex=cex )
  }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    pGtCompVal = sum( difference > compVal ) / length( difference )
    pLtCompVal = 1 - pGtCompVal
    lines( c(compVal,compVal) , c(0.96*cvHt,0) ,
           lty="dotted" , lwd=2 , col=cvCol )
    text( compVal , cvHt ,
          bquote( .(round(100*pLtCompVal,1)) * "% < " *
                    .(signif(compVal,3)) * " < " *
                    .(round(100*pGtCompVal,1)) * "%" ) ,
          adj=c(pLtCompVal,0) , cex=0.8*cex , col=cvCol )
    postSummary[,"compVal"] = compVal
    postSummary[,"pGtCompVal"] = pGtCompVal
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    pInROPE = ( sum( difference > ROPE[1] & difference < ROPE[2] )
                / length( difference ) )
    pGtROPE = ( sum( difference >= ROPE[2] ) / length( difference ) )
    pLtROPE = ( sum( difference <= ROPE[1] ) / length( difference ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          bquote( .(round(100*pLtROPE,1)) * "% < " * .(ROPE[1]) * " < " *
                    .(round(100*pInROPE,1)) * "% < " * .(ROPE[2]) * " < " *
                    .(round(100*pGtROPE,1)) * "%" ) ,
          adj=c(pLtROPE+.5*pInROPE,0) , cex=1 , col=ropeCol )

    postSummary[,"ROPElow"]=ROPE[1]
    postSummary[,"ROPEhigh"]=ROPE[2]
    postSummary[,"pLtROPE"]=pLtROPE
    postSummary[,"pInROPE"]=pInROPE
    postSummary[,"pGtROPE"]=pGtROPE
  }
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 , lend=1 )
  text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
        adj=c(.5,-1.7) , cex=cex )
  text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=F)

}


