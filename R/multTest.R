###############################################################
# for multiplicative batch effects in the residuals for each feature 
# after fitting linear mixed effects model
# Author: Joanne C. Beer, joannecbeer@gmail.com
###############################################################
# as described in the manuscript at 
# https://www.biorxiv.org/content/10.1101/868810v1
###############################################################
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license. 
###############################################################

multTest <- function(idvar, batchvar, features, 
                       formula, ranef, data, verbose=TRUE){
  ###########################################################
  # DATA SHOULD BE IN "LONG" FORMAT
  # PACKAGE DEPENDENCIES: lme4
  # INPUTS ##################################################
  # idvar:    name of ID variable (character string)
  # batchvar: name of the batch/site/scanner variable (character string)
  # features: vector of names of the feature variables (character string)
  #           or the numeric indices of the corresponding colunms
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
  # data:     name of the data.frame that contains the variables above
  #           rows are different subject/timepoints (long format), columns are different variables
  # verbose:  prints messages (logical TRUE/FALSE)
  # OUTPUTS #################################################
  # multEffect_ordered: table of Fligner-Killeen test results for each feature
  ###########################################################
  
  # make batch a factor if not already
  batch <- as.factor(data[,batchvar])
  if (verbose) cat("[multTest] found", nlevels(batch), 'batches\n')
  # feature names
  if (is.numeric(features[1])) {
    featurenames <- names(data)[features]
  } else {
    featurenames <- features
  }
  # number of features
  V <- length(featurenames)
  if (verbose) cat("[multTest] found", V, 'features\n')
  
  ##############################
  # fit model and do Fligner-Killeen test
  ##############################
  # make empty data structure to store results
  multEffect <- data.frame(feature=rep(NA, V), Chisq=rep(NA, V), df=rep(NA, V), p=rep(NA, V))
  # loop over features
  for (v in 1:V){ # begin loop over features
    if (verbose) cat("[multTest] testing for multiplicative batch effect for feature ", v, '\n')
    # create full model formula
    # full model includes a batch fixed effect
    full_formula <- as.formula(paste0(featurenames[v], '~', formula, '+' , batchvar, '+', ranef))
    # fit full model
    fit_full <- lme4::lmer(full_formula, data=data, REML=FALSE, control=lme4::lmerControl(optimizer='bobyqa'))
    # save residuals
    fit_residuals <- residuals(fit_full)
    # do Fligner-Killeen test
    FKtest <- fligner.test(fit_residuals ~ batch)
    # save results
    multEffect[v,] <- c(featurenames[v], FKtest$statistic, FKtest$parameter, FKtest$p.value)
  } # end loop over features
  # sort according to FK chi-squared statistic
  multEffect_ordered <- multEffect[order(as.numeric(multEffect$Chisq), decreasing=TRUE),]
  # add column names
  colnames(multEffect_ordered) <- c('Feature', 'ChiSq', 'DF', 'p-value')
  # return result
  return(multEffect_ordered)
}
