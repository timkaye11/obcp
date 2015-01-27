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
  output <- list(data=data, maxes=maxes, coeffs=coeffs, r_t=R, t=t, last_col=R[,t+1])
  return (structure(output, class="obcp_output"))
}

#' Add new data to the online bayesian changepoint detection model. Used to add new incoming data to the 
#' original model. 
#' 
#' Given a new data input, the posterior estimate is calculated, as well as the probability of a changepoint. The 
#' posterior is updated and the max value is recorded in the obcp object. 
#' 
#' @param obcp The output of the original obcp algorithm
#' @param data the data to be added to the model 
#' @keywords online, data, add data
#' @export
#' @examples
#' model <- obcp(data, 0, 1, 1, 1)
#' new_data <- c(24.4, 25.5, 26.6)
#' sapply(new_data, function(x) { add_data(model, x) }) 
#' model
#' 
add_data <- function(obcp, data) {
  if (length(data) > 1) { stop("data must be 1 observation")}
  if (length(obcp$last_col) < 20) { stop("it is recommended that you have at least 20 obs")}
  
  
  coeffs <- obcp$coeffs
  mu <- coeffs$mu
  kappa <- coeffs$kappa
  alpha <- coeffs$alpha
  beta <- coeffs$beta
  # calculate posterior
  pred_post <- dt_ls(data, 2*alpha, mu, beta*(kappa+1)/(alpha*kappa))
  
  # for changepoint detection
  last_col <- obcp$last_col
  t <- length(last_col)
  probs <- rep(0, t+1)
  
  # growth probability
  probs[2:(t+1)] <- last_col * pred_post * (1-H)
  # changepoint probability
  probs[1] <- sum( last_col * pred_post * H)

  # update the sufficient statistics (parameters) - Step 8 in Alg 1
  mu_new = c(mu0, (kappa*mu + data) / (kappa + 1))
  kappa_new <- c(kappa0, kappa + 1)
  alpha_new <- c(alpha0, alpha + 0.5)
  beta_new <- c(beta0, beta + (kappa *(data-mu)^2) / (2*(kappa+1)))
  
  obcp$coeffs$mu <<- mu_new
  obcp$coeffs$kappa <<- kappa_new
  obcp$coeffs$alpha <<- alpha_new
  obcp$coeffs$beta <<- beta_new
  obcp$maxes <<- c(obcp$maxes, which(probs == max(probs)))
  obcp$last_col <<- probs
  obcp$data <<- append(obcp$data, data)
  
  return (obcp)
}


#' Find where the changepoints occured. 
#' 
#' Detects the location of the changepoints from the data points that the model has seen so far
#' 
#' @param obcp the online bayesian changepoint detection model
#' @param num the number of changepoints to detect. sorted by p-value
#' @keywords find, changepoints
#' @export
#' @examples
#' model <- obcp(data, 0, 1, 1, 1)
#' # finds the top three changepoints
#' find_changepoints(model, 3)
find_changepoints <- function(obcp, num) {
  maxes <- obcp$maxes
  if (num > length(maxes)) stop("the number of changepoints cannot be longer than the data observed")
  diffed <- diff(maxes)
  sort_diff <- sort(diffed)
  locs <- match(sort_diff[1:num], diffed)
  return (locs)
}
  