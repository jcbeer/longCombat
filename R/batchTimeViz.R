#' Visualize Batches Over Time
#' 
#' \code{batchTimeViz} is a simple function that will visualize batches over time for multi-batch longitudinal data. Data should be in "long" format.
#' @param batchvar character string that specifies name of the batch variable. Batch variable should be a factor.
#' @param timevar character string that specifies name of numeric variable that distinguishes within-subject repeated measures, e.g., time, age, or visit. Will be plotted along x-axis.
#' @param data name of the data frame that contains the variables above. Rows are different observations (subject/timepoints), columns are different variables.
#' @param xlabel x-axis label (character string). Default is \code{'time'}.
#' @param ylabel y-axis label (character string). Default is \code{'batch'}.
#' @param title main title for the plot (character string). Default is no title.
#' @param verbose prints messages. Logical \code{TRUE} or \code{FALSE}. Default is \code{TRUE}.
#' @param ... other graphical parameter arguments passed to \code{\link[graphics]{par}}.
#' @return Creates a plot.

batchTimeViz <- function(batchvar, timevar, data, 
                         xlabel='time', ylabel='batch', title='', 
                         verbose=TRUE, ...){
  
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