#' Plot Individual Subject Trajectories
#' 
#' \code{trajPlot} will plot individual subject trajectories for multi-batch longitudinal data. Each line represents a trajectory of observations of a single feature (e.g., left fusiform cortical thickness) over time for an individual subject. Plotting points are coded by batch. Data should be in "long" format.
#' @param idvar character string that specifies name of ID variable. ID variable can be factor, numeric, or character. 
#' @param timevar character string that specifies name of numeric variable that distinguishes within-subject repeated measures, e.g., time, age, or visit. Will be plotted along x-axis.
#' @param feature character string that specifies name of the numeric feature variable to be plotted over time, or the numeric index of the corresponding column. Will be plotted along y-axis.
#' @param batchvar character string that specifies name of the batch variable. Batch variable should be a factor. Will determine the shape of plotting points.
#' @param data name of the data frame that contains the variables above. Rows are different observations (subject/timepoints), columns are different variables.
#' @param point.shape optional numeric vector that encodes the point shape for the batch variable levels. If not using default setting, length should be equal to number of unique observations included (\code{nrow(data)}) and in the same order as in the data frame. See the \code{pch} argument in \code{\link[graphics]{points}}. Default cycles through \code{pch=c(1:25, 0)} by batch level.
#' @param point.col optional vector that encodes colors of the point for each observation (character string of color names or hexadecimal codes). Length should be equal to number of unique observations included (\code{nrow(data)}) and in the same order as in the data frame. Default is \code{'black'} for all.
#' @param line.col optional vector that encodes colors of the line for each subject (character string of color names or hexadecimal codes). Length should be equal to number of unique subjects included (\code{length(unique(data[,idvar]))}) and in the same order as in the data frame. Default is \code{'black'} for all.
#' @param line.type optional numeric or character vector that encodes line type for each subject. Length should be equal to number of unique subjects included (\code{length(unique(data[,idvar]))}) and in the same order as in the data frame. See the \code{lty} argument in \code{\link[graphics]{par}}. Default is \code{'solid'} for all.
#' @param xlabel x-axis label (character string). Default is \code{'time'}.
#' @param ylabel y-axis label (character string). Default is \code{'feature'}.
#' @param title main title for the plot (character string). Default is no title.
#' @param xlimits two dimensional numeric vector giving minimum and maximum x-axis values. Default is minimum and maximum time variable values.
#' @param ylimits two dimensional numeric vector giving minimum and maximum y-axis values. Default is minimum and maximum feature variable values.
#' @param margins numeric vector of 4 values specifying margins in the order (bottom, left, top, right).
#' @param verbose prints messages. Logical \code{TRUE} or \code{FALSE}. Default is \code{TRUE}.
#' @param ... other graphical parameter arguments passed to \code{\link[graphics]{par}}.
#' @return Creates a plot.
#' 
#' @export

trajPlot <- function(idvar, timevar, feature, batchvar, data,
                     point.shape=NULL, point.col=NULL, 
                     line.col=NULL, line.type=NULL,
                     xlabel='time', ylabel='feature', title='', 
                     xlimits=NULL, ylimits=NULL, margins=c(5, 5, 3, 1),
                     verbose=TRUE, ...){
  
  # number of unique subjects 
  subjectIDs <- unique(data[,idvar])
  n <- length(subjectIDs)
  # number of observations
  L <- nrow(data)
  # make batch a factor if not already
  batch <- droplevels(as.factor(data[,batchvar]))
  # number of batches
  m <- nlevels(batch)
  # print messages
  if (verbose) cat("[trajPlot] found", n, 'unique subjects\n')
  if (verbose) cat("[trajPlot] found", L, 'observations\n')
  if (verbose) cat("[trajPlot] found", m, 'batches\n')

  # set limits for x-axis (time) and y-axis (feature)
  if (is.null(xlimits)){
    xlimits <- c(min(data[,timevar]), max(data[,timevar]))
  }
  if (is.null(ylimits)){
    ylimits <- c(min(data[,feature]), max(data[,feature]))
  }
  
  # make point shape and color vectors
  if (is.null(point.shape)) point.shape <- as.numeric(data[,batchvar]) %% 26
  if (is.null(point.col)) point.col <- rep('black', L)
  if (is.null(line.col)) line.col <- rep('black', n)
  if (is.null(line.type)) line.type <- rep('solid', n)
  data$point.shape <- point.shape
  data$point.col <- point.col
  
  ##############################
  # make plot
  ##############################
  par(mar=margins, ...)
  # plot first subject trajectory
  # get the data
  data.subject <- data[data[,idvar]==subjectIDs[1],]
  # sort by time point
  data.subject <- data.subject[order(data.subject[,timevar]),]
  # plot line
  plot(data.subject[,timevar], data.subject[,feature], col=line.col[1], lty=line.type[1], type='l', xlim=xlimits, ylim=ylimits, xlab=xlabel, ylab=ylabel, main=title, las=1)
  # plot points
  points(data.subject[,timevar], data.subject[,feature], pch=data.subject$point.shape, col=data.subject$point.col)
  for(i in 2:n){ # begin loop over remaining subjects
    # get the data
    data.subject <- data[data[,idvar]==subjectIDs[i],]
    # sort by time point
    data.subject <- data.subject[order(data.subject[,timevar]),]
    # plot line
    lines(data.subject[,timevar], data.subject[,feature], col=line.col[i], lty=line.type[i])
    # plot points
    points(data.subject[,timevar], data.subject[,feature], pch=data.subject$point.shape, col=data.subject$point.col)
  } # end loop over remaining subjects
}