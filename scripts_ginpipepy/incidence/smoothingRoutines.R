#' Smoothing over of the dated point estimates
#'
#' The function resamples a certain proportion of the given point estimates, if given with the respective weight, if not randomly.
#' For each resampling the values are collected, if no data is available for certain days, the values are interpolated.
#' The median and 5th and 95th percentiles at each time point are smoothed with a 2-sided convolution filter of the given window size.
#'
#' @param input.table a data frame containing a column "value" with the point estimates and a column "t" wiht the time points
#' @param ts a vector of time points which should be interpolated
#' @param weights (optional) a weight between 0 and 1 for each value in the table, for being sampled (default: all equal)
#' @param N_sampling (optional) number of resampling to draw the median and percentiles from (default: 1000)
#' @param p_sample (optional) proportion of point estimates to sample (default: 0.5)
#' @param width (optional) windowsize for smoothing (default: 7)
#' @return a data frame containing the smoothed values for each time point  ("t" = time points in ts, "smoothMedian" = median, "smooth5" = 5th percentile and "smooth95" = 95th percentile) 
smooth_point_estimates <- function (input.table, ts, weights=NULL, N_sampling=NULL, p_sample=NULL, width=NULL) {
  # number of resamplings
  if(is.null(N_sampling))
    N_sampling <- 1000
  # proportion of samples
  if(is.null(p_sample))
    p_sample <- 0.5
  #smoothing window of one week
  if(is.null(width))
  width <- 7

  # weights per point estimate
  if(is.null(weights))
    weights=rep(1, nrow(input.table))

  samplings <- matrix(ncol=length(ts), nrow=0)
  for(i in seq(N_sampling)) {
    # sample for each row if it is in the sampling set
    #sub_input.table <- input.table[runif(nrow(input.table))<p_sample, ]
    sub_value.table <- input.table[sample(seq(nrow(input.table)), round(nrow(input.table)*p_sample), replace = F, prob=NULL),]
    #sort by timepoint
    sub_value.table <- sub_value.table[order(sub_value.table$t),]

    #### first interpolate, collect and smooth median afterwards
    interpol <- approx(sub_value.table$t, sub_value.table$value, xout=ts)
    samplings <-rbind(samplings, interpol$y)
  }

  # taking the quantiles of all samplings
  quantiles <- apply(samplings, 2, quantile, c(0.05, 0.5, 0.95), na.rm=T)
  smoothed_quantiles <- t(apply(quantiles, 1, filter, filter=rep(1,width)/width, sides=2))
  df <- data.frame(t=ts, smoothMedian=smoothed_quantiles[2,], smooth5=smoothed_quantiles[1,], smooth95=smoothed_quantiles[3,])
  return(df)
}

