## Online Bayesian Change-point Detection in R

This algorithm is based off the paper written by ![Adams and Mackay](https://hips.seas.harvard.edu/files/adams-changepoint-tr-2007.pdf). 

Essentially, given time-series data, this package detects changepoints (i.e instances where the mean or the variance changed) using the Bayesian framework. 

## Quick Example

```r
devtools::install_github("timkaye11", "onlineBayesChangepoint")
library(onlineBayesChangepoint)
library(quantmod)

# from the quantmod package to load stock data
getSymbols("GOOG", src="yahoo")
data <- GOOG$GOOG.Adjusted
plot.ts(data)

# from the plot, the mean is around 500, and there are some definite changepoints
c_points <- obcp(data, mu0=500, kappa0=0, alpha0=7, beta0=4)

# to get a sense of r_t, the run length probabilities
plot.ts(c_points$maxes)

# because there are splits around 38th, 100th, 130th, 150th, and a few 
# around the 200th day in the  plot of the r_t, we see that 
# we correctly detected the changepoints. 

```

## TODO
 - be able to support data streams, online data coming in
 - functionality for other distributions, besides normal-inverse-gamma
 