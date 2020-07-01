###############################################################
# addTest function will test (Kenward-Roger method) 
# for additive batch effects in the residuals for each feature 
# after fitting linear mixed effects model
# Author: Joanne C. Beer, joannecbeer@gmail.com
###############################################################
# as described in the manuscript at 
# https://www.biorxiv.org/content/10.1101/868810v4
###############################################################
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license. 
###############################################################

addTest <- function(idvar, batchvar, features, 
                       formula, ranef, data, verbose=TRUE){
  ###########################################################
  # DATA SHOULD BE IN "LONG" FORMAT
  # PACKAGE DEPENDENCIES: lme4, pbkrtest
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
  # addEffect_ordered: table of Kenward-Roger test results for each feature
  ###########################################################
  
  # make batch a factor if not already
  batch <- as.factor(data[,batchvar])
  if (verbose) cat("[addTest] found", nlevels(batch), 'batches\n')
  # feature names
  if (is.numeric(features[1])) {
    featurenames <- names(data)[features]
  } else {
    featurenames <- features
  }
  # number of features
  V <- length(featurenames)
  if (verbose) cat("[addTest] found", V, 'features\n')
  
  ##############################
  # fit model and do Kenward Roger test
  ##############################
  # make empty data structure to store results
  addEffect <- data.frame(feature=rep(NA, V), KRFstat=rep(NA, V), KRddf=rep(NA, V), KR.p=rep(NA, V))
  # loop over features
  for (v in 1:V){ # begin loop over features
    if (verbose) cat("[addTest] testing for additive batch effect for feature ", v, '\n')
    # create full and reduced model formulas
    # full model includes a batch fixed effect
    full_formula <- as.formula(paste0(featurenames[v], '~', formula, '+' , batchvar, '+', ranef))
    # reduced model omits a batch fixed effect
    reduced_formula <- as.formula(paste0(featurenames[v], '~', formula, '+', ranef))
    # fit full and reduced models
    fit_full <- lme4::lmer(full_formula, data=data, REML=FALSE, control=lme4::lmerControl(optimizer='bobyqa'))
    fit_reduced <- lme4::lmer(reduced_formula, data=data, REML=FALSE, control=lme4::lmerControl(optimizer='bobyqa'))
    # do KR test (slowest step)
    KR <- pbkrtest::KRmodcomp(fit_full, fit_reduced)
    # save in results
    addEffect[v,] <- c(featurenames[v], KR$stats$Fstat, KR$stats$ddf, KR$stats$p.value)
  } # end loop over features
  # sort according to KR F statistic
  addEffect_ordered <- addEffect[order(as.numeric(addEffect$KRFstat), decreasing=TRUE),]
  # add column names
  colnames(addEffect_ordered) <- c('Feature', paste0('KR F(', KR$stats$ndf, ', KRddf)'), 'KRddf', 'KR p-value')
  # return result
  return(addEffect_ordered)
}
