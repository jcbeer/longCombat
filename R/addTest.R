#' Test for Additive Batch Effects
#' 
#' \code{addTest} function will test for additive batch effects in the residuals for each feature after fitting a linear mixed effects model. Uses Kenward-Roger method for significance testing. Data should be in "long" format. Depends on \code{lme4} and \code{pbkrtest} packages.
#' @param idvar character string that specifies name of ID variable. ID variable can be factor, numeric, or character. 
#' @param batchvar character string that specifies name of the batch variable. Batch variable should be a factor.
#' @param features character string that specifies names of the numeric feature variables, or the numeric indices of the corresponding columns.
#' @param formula character string representing all fixed effects on the right side of the formula for the linear mixed effects model. This should be in the notation used by \code{lme4} and include covariates, time, and any interactions. For example, \code{"age + sex + diagnosis*time"} fits model with fixed effects age, sex, diagnosis, time, and the diagnosis*time interaction. Formula should NOT include batchvar and should NOT include random effects.
#' @param ranef character string representing formula for the random effects in the notation used by \code{lme4}. For example, \code{"(1|subid)"} fits a random intercept for each unique idvar \code{subid}, and \code{"(1 + time|subid)"} fits a random intercept and random slope for each unique \code{subid}.
#' @param data name of the data frame that contains the variables above. Rows are different observations (subject/timepoints), columns are different variables.
#' @param verbose prints messages. Logical \code{TRUE} or \code{FALSE}. Default is \code{TRUE}.
#' @return A data frame of Kenward-Roger test results for each feature.

addTest <- function(idvar, batchvar, features,
                    formula, ranef, data, verbose=TRUE){
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