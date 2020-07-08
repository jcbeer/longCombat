#' Test for Multiplicative Batch Effects
#' 
#' \code{multTest} function will test for multiplicative batch effects in the residuals for each feature after fitting a linear mixed effects model. Uses Fligner-Killeen method for significance testing. Data should be in "long" format. Depends on \code{lme4} package.
#' @param idvar character string that specifies name of ID variable. ID variable can be factor, numeric, or character. 
#' @param batchvar character string that specifies name of the batch variable. Batch variable should be a factor.
#' @param features character string that specifies names of the numeric feature variables, or the numeric indices of the corresponding columns.
#' @param formula character string representing all fixed effects on the right side of the formula for the linear mixed effects model. This should be in the notation used by \code{lme4} and include covariates, time, and any interactions. For example, \code{"age + sex + diagnosis*time"} fits model with fixed effects age, sex, diagnosis, time, and the diagnosis*time interaction. Formula should NOT include batchvar and should NOT include random effects.
#' @param ranef character string representing formula for the random effects in the notation used by \code{lme4}. For example, \code{"(1|subid)"} fits a random intercept for each unique idvar \code{subid}, and \code{"(1 + time|subid)"} fits a random intercept and random slope for each unique \code{subid}.
#' @param data name of the data frame that contains the variables above. Rows are different observations (subject/timepoints), columns are different variables.
#' @param verbose prints messages. Logical \code{TRUE} or \code{FALSE}. Default is \code{TRUE}.
#' @return A data frame of Fligner-Killeen test results for each feature.
#' 
#' @export

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