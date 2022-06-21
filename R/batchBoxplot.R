#' Boxplot for Batch Effects
#' 
#' \code{batchBoxplot} function will plot residuals of linear mixed effects model for a single feature by batch to visualize additive and multiplicative batch effects. Data should be in "long" format. Depends on \code{lme4} package.
#' @param idvar character string that specifies name of ID variable. ID variable can be factor, numeric, or character.
#' @param batchvar character string that specifies name of the batch variable. Batch variable should be a factor.
#' @param feature character string that specifies name of the numeric feature variable, or the numeric index of the corresponding column.
#' @param formula character string representing all fixed effects on the right side of the formula for the linear mixed effects model. This should be in the notation used by \code{lme4} and include covariates, time, and any interactions. For example, \code{"age + sex + diagnosis*time"} fits model with fixed effects age, sex, diagnosis, time, and the diagnosis*time interaction. Formula should NOT include batchvar and should NOT include random effects.
#' @param ranef character string representing formula for the random effects in the notation used by \code{lme4}. For example, \code{"(1|subid)"} fits a random intercept for each unique idvar \code{subid}, and \code{"(1 + time|subid)"} fits a random intercept and random slope for each unique \code{subid}.
#' @param data name of the data frame that contains the variables above. Rows are different observations (subject/timepoints), columns are different variables.
#' @param adjustBatch should residuals be adjusted for the fixed effect of batch? Logical \code{TRUE} or \code{FALSE}. Use \code{FALSE} to illustrate additive (and multiplicative) batch effects. Use \code{TRUE} to illustrate only multiplicative batch effects. Default is \code{FALSE}.
#' @param orderby \code{'mean'} orders boxplots by increasing mean; best for illustrating additive batch effects (use with \code{adjustBatch=FALSE}). \code{'var'} orders boxplots by increasing variance; best for illustrating multiplicative batch effects. Default is \code{'mean'}.
#' @param plotMeans should batch means be plotted on top of the boxplots? Logical \code{TRUE} or \code{FALSE}. Default is \code{TRUE}.
#' @param colors vector of colors the same length and order as \code{levels(as.factor(data[,batchvar]))} that determines the colors of the boxplots (character string of color names or hexadecimal codes). Default is \code{"grey"} for all.
#' @param xlabel x-axis label (character string). Default is \code{'batch'}.
#' @param ylabel y-axis label (character string). Default is \code{'residuals'}.
#' @param ylim y-axis limits.
#' @param title main title for the plot, default is no title (character string).
#' @param verbose prints messages. Logical \code{TRUE} or \code{FALSE}. Default is \code{TRUE}.
#' @param ... other graphical parameter arguments passed to \code{\link[graphics]{par}}.
#' @return Creates a boxplot.
#' 
#' @export

batchBoxplot <- function(idvar, batchvar, feature, 
                         formula, ranef, data,
                         adjustBatch=FALSE, orderby='mean', 
                         plotMeans=TRUE, colors='grey',
                         xlabel='batch', ylabel='residuals',
                         ylim=NULL,
                         title='', 
                         verbose=TRUE, ...){
  
  # make batch a factor if not already
  data[,batchvar] <- droplevels(as.factor(data[,batchvar]))
  if (verbose) cat("[batchBoxplot] found", nlevels(data[,batchvar]), 'batches\n')
  # get feature names
  if (is.numeric(feature)) {
    featurename <- names(data)[feature]
  } else {
    featurename <- feature
  }
  # make color vector 
  if (length(colors) < nlevels(data[,batchvar])){
    colors <- rep_len(colors, length.out=nlevels(data[,batchvar]))
  }
  
  ##############################
  # fit linear mixed effect model
  ##############################
  if (verbose) cat(paste0('[longCombat] fitting lme model for feature ', feature, '\n'))
  # make the lmer formula
  if (adjustBatch==TRUE){
    lme_formula <- as.formula(paste0(featurename, '~', formula, '+' , batchvar, '+', ranef))
  } else if (adjustBatch==FALSE){
    lme_formula <- as.formula(paste0(featurename, '~', formula, '+', ranef))
  }
  # fit lme4 model
  lme_fit <- lme4::lmer(lme_formula, data=data, REML=TRUE, control=lme4::lmerControl(optimizer='bobyqa'))
  # save residuals and their means and variances
  fit_residuals <- data.frame(residuals=residuals(lme_fit), batch=data[,batchvar])
  fit_residuals_means <- aggregate(fit_residuals$residuals, by=list(fit_residuals$batch), FUN=mean)
  fit_residuals_var <- aggregate(fit_residuals$residuals, by=list(fit_residuals$batch), FUN=var)
  # order boxplots by mean or variance
  if (orderby=='mean'){
    batchorder <- with(fit_residuals, reorder(batch, residuals, mean))
    colors <- colors[order(fit_residuals_means[,2])]
  } else if (orderby=='var'){
    batchorder <- with(fit_residuals, reorder(batch, residuals, var))
    colors <- colors[order(fit_residuals_var[,2])]
  }
  
  ##############################
  # make plot
  ##############################
  par(mar=c(3, 5, 3, 1), ...)
  boxplot(fit_residuals$residuals ~ batchorder, main=title, ylab='', xlab='', 
          ylim=ylim, lty=1, col=colors, las=1, xaxt='n')
  if (plotMeans==TRUE){
    points(fit_residuals_means[,2][order(fit_residuals_means[,2])], pch=5, col='red', cex=0.6)
  }
  mtext(text=ylabel, side=2, line=3.5, cex=1.25, font=2)
  mtext(text=xlabel, side=1, line=1.5, cex=1.25, font=2)
  abline(h=0)
}