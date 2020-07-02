###############################################################
# batchBoxplot will plot residuals of linear mixed effect model 
# by batch to help assess additive and multiplicative batch effects
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
                         plotMeans=TRUE, orderby='mean', adjustBatch=FALSE,
                         xlabel='batch', ylabel='residuals',
                         title='', verbose=TRUE, ...){
  ###########################################################
  # DATA SHOULD BE IN "LONG" FORMAT
  # idvar:    name of ID variable (character string)
  # batchvar: name of the batch/site/scanner variable (character string)
  # feature:  name of the feature variable to plot (character string)
  #           or the numeric index of the corresponding column
  # xlabel:   x-axis label, default is 'time' (character string)
  # ylabel:   y-axis label, default is 'batch' (character string)
  # title:    main title for the plot, default is no title (character string)
  # data:     name of the data.frame that contains the variables above
  #           rows are different subject/timepoints (long format), columns are different variables
  # verbose:  prints messages (logical TRUE/FALSE)
  # ...:      other graphical parameter arguments passed to par()
  ###########################################################
  
  # make batch a factor if not already
  batch <- as.factor(data[,batchvar])
  if (verbose) cat("[batchViz] found", nlevels(batch), 'batches\n')
  # number of batches
  m <- nlevels(batch)
  # row IDs for each batch 
  batches <- lapply(levels(batch), function(x) which(batch==x))
  
  ##############################
  # order batches by earliest timepoint and duration
  ##############################
  batchtime_first <- c()
  batchtime_duration <- c()
  for(i in 1:m){ # begin loop over batches
    # get earliest batch time point
    batchtime_first[i] <- min(data[batches[[i]],timevar])
    # get duration
    batchtime_duration[i] <- max(data[batches[[i]],timevar]) - min(data[batches[[i]],timevar])
  } # end loop over batches
  
  # reorder batches 
  batch <- factor(batch, levels=levels(batch)[order(batchtime_first, batchtime_duration)])
  # new row IDs for each batch 
  batches <- lapply(levels(batch), function(x) which(batch==x))
  
  # set limits for x axis (time)
  xlimits <- c(min(data[,timevar]), max(data[,timevar]))
  # limits for y axis
  ylimits <- c(0, m)
  
  ##############################
  # make plot
  ##############################
  par(mar=c(5, 3, 4, 2), ...)
  # empty plot
  plot(0, 0, xlim=xlimits, ylim=ylimits, col='white', xlab=xlabel, ylab='', main=title, yaxt='n')
  mtext(ylabel, side=2, line=1)
  for(i in 1:m){ # begin loop over batches
    # get batch time points
    batchdata <- data[batches[[i]],timevar]
    # plot line for batch
    lines(sort(batchdata), rep(i, length(batchdata)), type='o')
  } # end loop over batches
  
}