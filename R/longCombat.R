###############################################################
# longCombat function will implement Longitudinal ComBat
# harmonization for multi-batch longitudinal data
# Author: Joanne C. Beer, joannecbeer@gmail.com
###############################################################
# as described in the manuscript at 
# https://www.biorxiv.org/content/10.1101/868810v1
###############################################################
# This is a modification of the ComBat function code 
# from the sva package that can be found at
# https://bioconductor.org/packages/release/bioc/html/sva.html 
# and the combat.R that can be found at 
# https://github.com/Jfortin1/ComBatHarmonization
###############################################################
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license. 
###############################################################

longCombat <- function(idvar, batchvar, features, 
                       formula, ranef, niter=30, data, verbose=TRUE){
  ###########################################################
  # idvar:    name of ID variable (character string)
  # batchvar: name of the batch/site/scanner variable (character string)
  # features: vector of names of the feature variables (character string)
  #           or their numeric indices of the corresponding colunms
  # formula:  character string representing everything on the right side of the formula
  #           for the model, in the notation used by lm or lme4
  #           including covariates, time, and any interactions
  #           e.g. "age + sex + diagnosis*time"
  #           fits model with main effects age, sex, diagnosis, and time
  #           and the diagnosis*time interaction
  #           should NOT include batchvar
  #           should NOT include random effects 
  # ranef:    character string representing formula for the random effects
  #           in the notation used by lme4
  #           e.g. "(1|subid)" fits a random intercept for each unique idvar "subid"
  #           e.g. "(1 + time|subid)" fits a random intercept and slope unique "subid"
  # niter:    number of iterations for empirical Bayes step
  #           usually converges quickly in less than 30 iterations
  # data:     name of the data.frame that contains the variables above
  #           rows are different subject/timepoints, columns are different variables
  # verbose:  prints messages (logical TRUE/FALSE)
  ###########################################################
  
  # make batch a factor if not already
  batch <- as.factor(data[,batchvar])
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
  # total number of observations
  L <- nrow(data)
  
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
    # save sigma estimates 
    corr_estimates <- as.data.frame(lme4::VarCorr(lme_fit))
    sigma_estimates[v] <- corr_estimates[corr_estimates$grp=='Residual','sdcor']
    # save batch effects
    batch_effects[,v] <- fixef(lme_fit)[grep(batchvar, names(fixef(lme_fit)))]
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
  # add gamma1hat to the rest of the scanner effect table
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
  # add names
  ##############################
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
              gammahat=gammahat, delta2hat=delta2hat,
              gammastarhat=gammastarhat_final, delta2starhat=delta2starhat_final,
              gammabar=gammabar, tau2bar=tau2bar
              ))
}
