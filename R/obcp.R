obcp_output <- structure(list(), class = "obcp_output")
#'
#' Performs the algorithm outlined in the paper by Adams and MacKay. 
#' 
#' Input must be a vector. The function outputs relevant data from the simulation, including the r_t table
#'
#' @param data A vector of data to be analyzed. Must be a one-dimensional time series
#' @param mu0 The prior for mu. The default is set to 0, but it may not be suitable for all data
#' @param kappa0 The prior for kappa. The default is set to 1. 
#' @param alpha0 The prior for alpha in the inverse-gamma distribution
#' @param beta0 The prior for beta in the inverse-gamma distribution
#' 
#' @keywords online, bayesian, changepoint, normal-inverse-gamma
#' @return a S3 object containing the r_t, maxes, and other information. 
#' @export
#' @examples
#' output <- obcp(data, mu=0, kappa=2, alpha=2, beta=4)
#' summary(output)
obcp <- function(data, mu0, kappa0, alpha0, beta0) {
  # for referencing data size
  t <- length(data)
  
  # parameters for the normal-inverse-gamma distribution
  mu <- mu0
  kappa <- kappa0
  alpha <- alpha0
  beta <- beta0
  
  # matrix containing the run length probabilities (r_t)
  R <- matrix(0, nrow=t+1, ncol=t+1)
    
  # at t=1, the run length is definitely 1. Initialize table 
  R[1,1] <- 1
  
  # keep track of the maximum probabilities
  maxes <- rep(0, t+1)
  
  for (i in 1:t) {  
    # the predictive posterior distribution for each new data observation under each of the params
    pred_post <- dt_ls(data[i], 2*alpha, mu, beta*(kappa+1)/(alpha*kappa))
    
    # calculate the hazard function for this interval
    H <- hazard_func(i, 100)
    
    # calculate the growth probability (step 4 in Algorithm 1)
    R[2:(i+1), i+1] <- R[1:i, i] * pred_post * (1-H)
    
    # calculate the changepoint probability (step 5 in Algorithm 1)
    # this is the probability that there was an observed changepoint
    R[1, i+1] <- sum( R[1:i, i] * pred_post * H )
    
    # normalize run length probability - improves accuracy
    R[,i+1] = R[,i+1] / sum(R[,i+1])
    
    
    # update the sufficient statistics (parameters) - Step 8 in Alg 1
    mu_new = c(mu0, (kappa*mu + data[i]) / (kappa + 1))
    kappa_new <- c(kappa0, kappa + 1)
    alpha_new <- c(alpha0, alpha + 0.5)
    beta_new <- c(beta0, beta + (kappa *(data[i]-mu)^2) / (2*(kappa+1)))
    
    mu <- mu_new
    kappa <- kappa_new
    alpha <- alpha_new
    beta <- beta_new
    
    # store the max value
    maxes[i] <- which(R[,i] == max(R[,i]))
  }
  #print(R)
  coeffs <- list(mu=mu, kappa=kappa, alpha=alpha, beta=beta)
  output <- list(data=data, maxes=maxes, coeffs=coeffs, r_t=R, t=t)
  return (structure(output, class="obcp_output"))
}
