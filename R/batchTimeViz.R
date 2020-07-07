#' Visualize Batches Over Time
#' 
#' \code{batchTimeViz} is a simple function that will visualize batches over time for multi-batch longitudinal data. Data should be in "long" format.
#' @param batchvar name of the batch/site/scanner variable (character string).
#' @param timevar name of the time variable, e.g. age or time from baseline (character string).
#' @param xlabel x-axis label, default is \code{'time'} (character string).
#' @param ylabel y-axis label, default is \code{'batch'} (character string).
#' @param title main title for the plot, default is no title (character string).
#' @param data name of the data frame that contains the variables above rows are different subject/timepoints (long format), columns are different variables.
#' @param verbose prints messages (logical \code{TRUE} or \code{FALSE}). Default is \code{TRUE}.
#' @param ... other graphical parameter arguments passed to \code{par()}.
#' @return Creates a plot.

batchTimeViz <- function(batchvar, timevar, xlabel='time', ylabel='batch', 
                     title='', data, verbose=TRUE, ...){
  
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