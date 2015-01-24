summary <- function(x) UseMethod("summary")
#' Get a summary of the online bayesian changepoint output returned from the simulation
#' 
#' Returns a data table of relevatnt information
#' 
#' @param x The output of the obcp function. A class of "obcp_output"
#' @return a data frame of relevant info. 
#' @keywords summary, obcp
#' 
#' @export
summary.obcp_output <- function(x) {
  diff_maxes <- diff(x$maxes)
  cutoff <- mean(range(maxes))
  locations <- which(abs(diff_maxes) > cutoff)
  x$locations <- locations
  return (data.frame(x))
}

