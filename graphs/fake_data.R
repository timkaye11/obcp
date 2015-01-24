# Function to make fake data, for testing purposes
# The value of lambda somewhat controls the number of breakpoints that are detected
fake_data <- function(num_iters, harzard_func=hazard_func, lam=200) {
  # parameters for normal-inverse-gamma distribution
  mu0 <- 0
  kappa0 <- 1
  alpha0 <- 1
  beta0 <- 1
  # holds the data
  X <- rep(1, num_iters)
  # keeps track of the changepoints
  c_points <- c()
  
  new_ivar <- rgamma(1, alpha0) / beta0
  new_mean <- (kappa0*new_ivar)^(-0.5)*rnorm(1) + mu0
  new_run <- 0
  for (i in 1:num_iters) {
    p <- hazard_func(i, lam) 
    if (runif(1) < p) {
      print(paste("\nChangepoint",i))
      new_ivar <- rgamma(1, alpha0) * beta0
      new_mean <- (kappa0 * new_ivar)^(-0.5)*rnorm(1) + mu0
      new_run <- 0
      c_points <- c(c_points, i)
    } else {
      new_run <- new_run + 1
    }
    # draw data based on the distribution
    X[i] <- new_ivar^(-0.5) * rnorm(1) + new_mean
  }
  # keep track of the changepoints
  attr(X, "num_cp") <- length(c_points)
  attr(X, "changepoints") <- c_points
  return (X)
}

