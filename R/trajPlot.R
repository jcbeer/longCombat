#' Plot Individual Subject Trajectories
#' 
#' \code{trajPlot} will plot individual subject trajectories for multi-batch longitudinal data. Each line represents a trajectory of observations of a single feature (e.g., left fusiform cortical thickness) over time for an individual subject. Plotting points are coded by batch. Data should be in "long" format.
#' @param idvar name of ID variable (character string).
#' @param timevar name of variable that distinguishes within-subject repeated measures, e.g., time, age, or visit (character string). Will be plotted along x-axis.
#' @param feature name of the feature variable (character string) to be plotted over time, or the numeric index of the corresponding column. Will be plotted along y-axis.
#' @param batchvar name of the batch/site/scanner variable (character string). Will determine the shape of plotting points.
#' @param point.shape vector that encodes the point shape for the batch variable levels. Should be in the same order as \code{levels(data[subset, batchvar])}. See the \code{pch} argument in \code{?points}. Default cycles through \code{pch=c(1:25, 0)} by batch level.
#' @param point.col vector that encodes colors of the point for each observation (character string of color names or hexadecimal codes). Length should be equal to number of unique observations included (\code{nrow(data[subset,])}) and in the same order as in the data frame. Default is \code{'black'} for all.
#' @param line.col vector that encodes colors of the line for each subject (character string of color names or hexadecimal codes). Length should be equal to number of unique subjects included (\code{length(unique(data[subset,idvar]))}) and in the same order as in the data frame. Default is \code{'black'} for all.
#' @param xlabel x-axis label, default is \code{'time'} (character string).
#' @param ylabel y-axis label, default is \code{'feature'} (character string).
#' @param title main title for the plot, default is no title (character string).
#' @param xlimits two dimensional numeric vector giving minimum and maximum x-axis values.
#' @param ylimits two dimensional numeric vector giving minimum and maximum y-axis values.
#' @param data name of the data frame that contains the variables above. Rows are different observations (subject/timepoints), columns are different variables.
#' @param subset logical or numeric 0/1 vector the same length as \code{nrow(data)} denoting which observations (i.e., rows) to plot. Default is to include all rows.
#' @param verbose prints messages (logical \code{TRUE} or \code{FALSE}). Default is \code{TRUE}.
#' @param ... other graphical parameter arguments passed to \code{par()}.
#' @return Creates a plot.

trajPlot <- function(idvar, timevar, feature, batchvar, 
                     point.shape=NULL, point.col=NULL, line.col=NULL, 
                     xlabel='time', ylabel='feature', title='', 
                     xlimits=NULL, ylimits=NULL, margins=c(5, 5, 3, 1),
                     data, subset=rep(TRUE, nrow(data)), verbose=TRUE, ...){
  
  # get data subset
  if (is.null(subset)) subset <- rep(TRUE, nrow(data))
  data.subset <- data[subset,]
  # number of unique subjects 
  subjectIDs <- unique(data.subset[,idvar])
  n <- length(subjectIDs)
  # number of observations
  L <- nrow(data.subset)
  # make batch a factor if not already
  batch <- droplevels(as.factor(data.subset[,batchvar]))
  # number of batches
  m <- nlevels(batch)
  # print messages
  if (verbose) cat("[trajPlot] found", n, 'unique subjects\n')
  if (verbose) cat("[trajPlot] found", L, 'observations\n')
  if (verbose) cat("[trajPlot] found", m, 'batches\n')

  # set limits for x-axis (time) and y-axis (feature)
  if (is.null(xlimits)){
    xlimits <- c(min(data.subset[,timevar]), max(data.subset[,timevar]))
  }
  if (is.null(ylimits)){
    ylimits <- c(min(data.subset[,feature]), max(data.subset[,feature]))
  }
  
  # make point shape and color vectors
  if (is.null(point.shape)) point.shape <- as.numeric(data.subset[,batchvar]) %% 26
  if (is.null(point.col)) point.col <- rep('black', L)
  if (is.null(line.col)) line.col <- rep('black', n)
  data.subset$point.shape <- point.shape
  data.subset$point.col <- point.col
  
  ##############################
  # make plot
  ##############################
  par(mar=margins, ...)
  # plot first subject trajectory
  # get the data
  data.subject <- data.subset[data.subset[,idvar]==subjectIDs[1],]
  # sort by time point
  data.subject <- data.subject[order(data.subject[,timevar]),]
  # plot line
  plot(data.subject[,timevar], data.subject[,feature], col=line.col[1], type='l', xlim=xlimits, ylim=ylimits, xlab=xlabel, ylab=ylabel, main=title, las=1)
  # plot points
  points(data.subject[,timevar], data.subject[,feature], pch=data.subject$point.shape, col=data.subject$point.col)
  for(i in 2:n){ # begin loop over remaining subjects
    # get the data
    data.subject <- data.subset[data.subset[,idvar]==subjectIDs[i],]
    # sort by time point
    data.subject <- data.subject[order(data.subject[,timevar]),]
    # plot line
    lines(data.subject[,timevar], data.subject[,feature], col=line.col[i])
    # plot points
    points(data.subject[,timevar], data.subject[,feature], pch=data.subject$point.shape, col=data.subject$point.col)
  } # end loop over remaining subjects
}
