#' Harmonize Multi-batch Longitudinal Data
#' 
#' \code{longCombat} function will implement longitudinal ComBat harmonization for multi-batch longitudinal data. Longitudinal ComBat uses an empirical Bayes method to harmonize means and variances of the residuals across batches in a linear mixed effects model framework. Detailed methods are described in the manuscript at \url{https://www.biorxiv.org/content/10.1101/868810v4}. This is a modification of the ComBat function code from the \code{sva} package that can be found at \url{https://bioconductor.org/packages/release/bioc/html/sva.html} and \code{combat.R} that can be found at \url{https://github.com/Jfortin1/ComBatHarmonization}. Data should be in "long" format. Depends on \code{lme4} package.
#' @param idvar character string that specifies name of ID variable. ID variable can be factor, numeric, or character. 
#' @param timevar character string that specifies name of numeric variable that distinguishes within-subject repeated measures, e.g., time, age, or visit.
#' @param batchvar character string that specifies name of the batch variable. Batch variable should be a factor.
#' @param features character string that specifies names of the numeric feature variables, or the numeric indices of the corresponding columns.
#' @param formula character string representing all fixed effects on the right side of the formula for the linear mixed effects model. This should be in the notation used by \code{lme4} and include covariates, time, and any interactions. For example, \code{"age + sex + diagnosis*time"} fits model with fixed effects age, sex, diagnosis, time, and the diagnosis*time interaction. Formula should NOT include batchvar and should NOT include random effects.
#' @param ranef character string representing formula for the random effects in the notation used by \code{lme4}. For example, \code{"(1|subid)"} fits a random intercept for each unique idvar \code{subid}, and \code{"(1 + time|subid)"} fits a random intercept and random slope for each unique \code{subid}.
#' @param data name of the data frame that contains the variables above. Rows are different observations (subject/timepoints), columns are different variables.
#' @param niter number of iterations for empirical Bayes step. Usually converges quickly in less than 30 iterations. Default is 30.
#' @param method method for estimating sigma in standardization step (character string). \code{'REML'} (default, more conservative type I error control) or \code{'MSR'} (more powerful, less conservative type I error control).
#' @param verbose prints messages. Logical \code{TRUE} or \code{FALSE}. Default is \code{TRUE}.
#' @return Function outputs a list including the following:
#' \describe{
#'     \item{\code{data_combat}}{data frame with columns idvar, timevar, and ComBat-harmonized data for each feature}
#'     \item{\code{gammahat}}{data frame containing mean of standardized data for each batch (row) and feature (column)}
#'     \item{\code{delta2hat}}{data frame containing variance of standardized data for each batch (row) and feature (column)}
#'     \item{\code{gammastarhat}}{data frame containing empirical Bayes estimate of additive batch effects}
#'     \item{\code{delta2starhat}}{data frame containing empirical Bayes estimate of multiplicative batch effects}
#'     }
#'     
#' @export

longCombat <- function(idvar, timevar, batchvar, features, 
                       formula, ranef, data, niter=30, method='REML', verbose=TRUE){

  # check for missing data 
  if (sum(is.na(data)) > 0) {
    missing <- paste(names(data)[apply(data, 2, function(x) sum(is.na(x)) > 0)], collapse=', ')
    message <- paste0('Missing data in variables:\n\n', missing, '\n\nBefore running longCombat either impute the missing values, remove those rows, or remove columns with missing values if that variable is not in the model.')
    stop(message)
  }
  # make batch a factor if not already
  batch <- droplevels(as.factor(data[,batchvar]))
  # check for batches with only one observation
  if (min(table(batch)) <= 1) {
    batch_single <- paste(names(table(batch))[table(batch) <= 1], collapse=', ')
    message <- paste0('The following batches have only one observation:\n\n', batch_single, '\n\nlongCombat needs at least 2 observations per batch to harmonize variance across batch. Remove rows for these batches before running longCombat.')
    stop(message)
  }
  if (verbose) cat("[longCombat] found", nlevels(batch), 'batches\n')
  # number of batches
  m <- nlevels(batch)
  # row IDs for each batch 
  batches <- lapply(levels(batch), function(x) which(batch==x))
  # number of observations for each batch
  ni <- sapply(batches, length)
  # feature names
  if (is.numeric(features[1])) {
    featurenames <- names(data)[features]
  } else {
    featurenames <- features
  }
  # number of features
  V <- length(featurenames)
  if (verbose) cat("[longCombat] found", V, 'features\n')
  # total number of observations
  L <- nrow(data)
  if (verbose) cat("[longCombat] found", L, 'total observations\n')
  
  ##############################
  # standardize data across features
  ##############################
  if (verbose) cat('[longCombat] standardizing data across features...\n')
  # make empty data structures to store results
  sigma_estimates <- rep(NA, V)
  predicted <- matrix(nrow=L, ncol=V)
  batch_effects <- matrix(nrow=(m-1), ncol=V)
  for (v in 1:V){ # begin loop over features
    if (verbose) cat(paste0('[longCombat] fitting lme model for feature ', v, '\n'))
    # make the linear mixed effects model lmer formula
    lme_formula <- as.formula(paste0(featurenames[v], '~', formula, '+' , batchvar, '+', ranef))
    # fit lme4 model
    lme_fit <- lme4::lmer(lme_formula, data=data, REML=TRUE, control=lme4::lmerControl(optimizer='bobyqa'))
    # save sigma estimate
    if (method == 'REML'){
      corr_estimates <- as.data.frame(lme4::VarCorr(lme_fit))
      sigma_estimates[v] <- corr_estimates[corr_estimates$grp=='Residual','sdcor']
    } else if (method == 'MSR'){
      resid <- residuals(lme_fit)
      sigma_estimates[v] <- sqrt((sum((resid-mean(resid))^2)/length(resid)))
    }
    # save batch effects
    batch_effects[,v] <- lme4::fixef(lme_fit)[grep(batchvar, names(lme4::fixef(lme_fit)))]
    # save predicted values
    predicted[,v] <- fitted(lme_fit)
  } # end loop over features
  # create a L*V matrix of sigma estimates
  sigmas <- matrix(rep(sigma_estimates, each=L), nrow=L, ncol=V)
  # create a L*V matrix of batch effects
  # incorporate constraint (sum_i ni * hat{gamma}_iv = 0) 
  # to get adjusted batch effect estimates
  # calculate the gamma1 hats 
  gamma1hat <- -(ni[2:m] %*% batch_effects)/L
  # add gamma1hat to the rest of the batch effect table
  batch_effects_adjusted <- sweep(batch_effects, 2, gamma1hat, FUN='+')
  # add gamma1hat as the top row
  batch_effects_adjusted <- rbind(gamma1hat, batch_effects_adjusted)
  # expand the adjusted batch effects to all timepoints
  batch_effects_expanded <- matrix(nrow=L, ncol=V)
  for(i in 1:m){ # begin loop over batches
    batch_effects_expanded[batches[[i]],] <- matrix(
      rep(batch_effects_adjusted[i,],length(batches[[i]])),
      ncol=V, byrow=TRUE) 
  } # end loop over batches
  # standardize the data
  data_std <- (data[,featurenames] - predicted + batch_effects_expanded) / sigmas
  
  ##############################
  # method of moments to estimate hyperparameters
  ##############################
  if (verbose) cat('[longCombat] using method of moments to estimate hyperparameters\n')
  gammahat <- matrix(nrow=m, ncol=V)
  delta2hat <- matrix(nrow=m, ncol=V)
  for (i in 1:m){ # begin loop over batches
      gammahat[i,] <- colMeans(data_std[batches[[i]],])
      delta2hat[i,] <- apply(data_std[batches[[i]],], 2, var)
  } # end loop over batches
  gammabar <- rowMeans(gammahat)
  tau2bar <- apply(gammahat, 1, var)
  Dbar <- rowMeans(delta2hat)
  S2bar <- apply(delta2hat, 1, var)
  # inverse gamma parameters
  lambdabar <- (Dbar^2 + 2*S2bar) / S2bar
  thetabar <- (Dbar^3 + Dbar*S2bar) / S2bar
  
  ##############################
  # empirical Bayes to estimate batch effects
  ##############################
  if (verbose) cat('[longCombat] using empirical Bayes to estimate batch effects...\n')
  if (verbose) cat('[longCombat] initializing...\n')
  # get initial estimates
  gammastarhat0 <- matrix(nrow=m, ncol=V)
  for (v in 1:V){ # begin loop over features
    gammastarhat0[,v] <- ((ni * tau2bar * gammahat[,v]) + (delta2hat[,v] * gammabar))/((ni * tau2bar) + delta2hat[,v])
  } # end loop over features
  delta2starhat0 <- matrix(nrow=m, ncol=V)
  for (v in 1:V){ # begin loop over features
    for(i in 1:m){ # begin loop over batches
      zminusgammastarhat2 <- sum((data_std[batches[[i]],v] - gammastarhat0[i,v])^2)
      delta2starhat0[i,v] <- (thetabar[i] + 0.5*zminusgammastarhat2) / (ni[i]/2 + lambdabar[i] - 1)
    } # end loop over features
  } # end loop over batches
  # iterate
  gammastarhat <- array(dim=c(m, V, (niter+1)))
  gammastarhat[,,1] <- gammastarhat0
  delta2starhat <- array(dim=c(m, V, (niter+1)))
  delta2starhat[,,1] <- delta2starhat0
  for(b in 2:(niter+1)){ # begin loop over iterations
    if (verbose) cat(paste0('[longCombat] starting EM algorithm iteration ', (b-1), '\n')) 
    for (v in 1:V){ # begin loop over features
      gammastarhat[,v,b] <- ((ni * tau2bar * gammahat[,v]) + (delta2starhat[,v,(b-1)] * gammabar))/((ni * tau2bar) + delta2starhat[,v,(b-1)])
      for(i in 1:m){ # begin loop over batches
        zminusgammastarhat2 <- sum((data_std[batches[[i]],v] - gammastarhat[i,v,(b-1)])^2)
        delta2starhat[i,v,b] <- (thetabar[i] + 0.5*zminusgammastarhat2) / (ni[i]/2 + lambdabar[i] - 1)
      } # end loop over batches
    } # end loop over features
  } # end loop over iterations
  # save final result
  gammastarhat_final <- gammastarhat[,,niter+1]
  delta2starhat_final <- delta2starhat[,,niter+1]

  ##############################
  # adjust data for batch effects
  ##############################
  if (verbose) cat('[longCombat] adjusting data for batch effects\n')
  # repeat each row the correct number of times
  gammastarhat_expanded <- matrix(nrow=L, ncol=V)
  delta2starhat_expanded <- matrix(nrow=L, ncol=V)
  for(i in 1:m){ # loop over batches
    gammastarhat_expanded[batches[[i]],] <- matrix(
      rep(gammastarhat_final[i,],length(batches[[i]])),
      ncol=V, byrow=TRUE) 
    delta2starhat_expanded[batches[[i]],] <- matrix(
      rep(delta2starhat_final[i,],length(batches[[i]])),
      ncol=V, byrow=TRUE) 
  } # end loop over batches
  # do ComBat 
  data_combat <- (sigmas/sqrt(delta2starhat_expanded))*(data_std - gammastarhat_expanded) + predicted - batch_effects_expanded
  
  ##############################
  # label the data
  ##############################
  # add IDs, time variable, and batch variable to data_combat
  data_combat <- cbind(data[,c(idvar, timevar, batchvar)], data_combat)
  # add names
  colnames(data_combat) <- c(idvar, timevar, batchvar, paste0(featurenames, '.combat'))
  colnames(gammahat) <- featurenames
  colnames(delta2hat) <- featurenames
  colnames(gammastarhat_final) <- featurenames
  colnames(delta2starhat_final) <- featurenames
  rownames(gammahat) <- levels(batch)
  rownames(delta2hat) <- levels(batch)
  rownames(gammastarhat_final) <- levels(batch)
  rownames(delta2starhat_final) <- levels(batch)
  
  ##############################
  # return results
  ##############################
  return(list(data_combat=data_combat,
              gammahat=gammahat, 
              delta2hat=delta2hat,
              gammastarhat=gammastarhat_final, 
              delta2starhat=delta2starhat_final
              ))
}