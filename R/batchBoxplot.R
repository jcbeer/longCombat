###############################################################
# batchBoxplot will plot residuals of linear mixed effect model
# for a simgle feature by batch 
# to visualize additive and multiplicative batch effects
# Author: Joanne C. Beer, joannecbeer@gmail.com
###############################################################
# as described in the manuscript at 
# https://www.biorxiv.org/content/10.1101/868810v4
###############################################################
# The present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license. 
###############################################################

batchBoxplot <- function(idvar, batchvar, feature, 
                         formula, ranef, data,
                         adjustBatch=FALSE, orderby='mean', 
                         plotMeans=TRUE, colors='grey',
                         xlabel='batch', ylabel='residuals',
                         title='', verbose=TRUE, ...){
  ###########################################################
  # DATA SHOULD BE IN "LONG" FORMAT
  # idvar:    name of ID variable (character string)
  # batchvar: name of the batch/site/scanner variable (character string)
  # feature:  name of the feature variable to plot (character string)
  #           or the numeric index of the corresponding column
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
  #           e.g. "(1 + time|subid)" fits a random intercept and slope for unique "subid"
  # data:     name of the data.frame that contains the variables above
  #           rows are different subject/timepoints (long format)
  #           columns are different variables
  # adjustBatch: should residuals be adjusted for batch? (logical TRUE/FALSE)
  #           use FALSE to illustrate additive (and multiplicative) batch effects
  #           use TRUE to illustrate only multiplicative batch effects
  # orderby:  'mean' orders boxplots by increasing mean
  #           best for illustrating additive batch effects (use with adjustBatch=FALSE)
  #           'var' orders boxplots by increasing variance
  #           best for illustrating multiplicative batch effects
  # plotMeans: should batch means be plotted on top of the boxplots (logical TRUE/FALSE)
  # colors:   vector of colors the same length and order as levels(as.factor(data$batchvar))
  #           that determines the colors of the boxplots 
  #           (character string of color names or hexadecimal codes)
  # xlabel:   x-axis label, default is 'time' (character string)
  # ylabel:   y-axis label, default is 'batch' (character string)
  # title:    main title for the plot, default is no title (character string)
  # verbose:  prints messages (logical TRUE/FALSE)
  # ...:      other graphical parameter arguments passed to par()
  ###########################################################
  
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
  # fit liner mixed effect model
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
  boxplot(fit_residuals$residuals ~ batchorder, main=title, ylab='', xlab='', lty=1, col=colors[batchorder], las=1, xaxt='n')
  if (plotMeans==TRUE){
    points(fit_residuals_means[,2][order(fit_residuals_means[,2])], pch=5, col='red', cex=0.6)
  }
  mtext(text=ylabel, side=2, line=3.5, cex=1.25, font=2)
  mtext(text=xlabel, side=1, line=1.5, cex=1.25, font=2)
  abline(h=0)
}