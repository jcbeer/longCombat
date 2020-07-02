#' Test for Multiplicative Batch Effects
#' 
#' \code{multTest} function will test (Fligner-Killeen method) for multiplicative batch effects in the residuals for each feature after fitting a linear mixed effects model. Data should be in "long" format. Depends on \code{lme4} package.
#' @param idvar name of ID variable (character string).
#' @param batchvar name of the batch/site/scanner variable (character string).
#' @param features vector of names of the feature variables (character string) or the numeric indices of the corresponding columns.
#' @param formula character string representing everything on the right side of the formula for the model, in the notation used by \code{lme4} including covariates, time, and any interactions, e.g., \code{"age + sex + diagnosis*time"} fits model with main effects age, sex, diagnosis, and time and the diagnosis*time interaction. Formula should NOT include batchvar and should NOT include random effects.
#' @param ranef character string representing formula for the random effects in the notation used by lme4, e.g., \code{"(1|subid)"} fits a random intercept for each unique idvar \code{"subid"}, and \code{"(1 + time|subid)"} fits a random intercept and slope unique \code{"subid"}.
#' @param data name of the data.frame that contains the variables above. Rows are different subject/timepoints (long format), columns are different variables.
#' @param verbose prints messages (logical \code{TRUE} or \code{FALSE}).
#' @return A data frame of Fligner-Killeen test results for each feature.

multTest <- function(idvar, batchvar, features, 
                       formula, ranef, data, verbose=TRUE){
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