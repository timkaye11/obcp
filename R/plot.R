plot <- function(x) UseMethod("plot")
#' Plot the changepoints from the algorithm
#' 
#' Plots the changepoints on top of the original data using the ggplot2 library. 
#' 
#' @param x the output from the obcp function
#' @keywords plot
#' @export
#' @examples
#' output <- obcp(data)
#' plot(output)
plot.obcp_output <- function(x) {
  data <- x$data
  if (!is.data.frame(data)) { stop("data must be data frame") }
  if (length(dim(data)) > 2) { stop("data must be 2 dimensional") }
  col <- which(sapply(data, is.numeric) == TRUE)
  
  p <- ggplot2::ggplot(data, aes_string(names(data)[col], names(data)[-col])) + geom_line()
  return (p)
}