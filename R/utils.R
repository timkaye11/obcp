# location-scale t-distrubution functions
# not provided in the standard t-distribution in R
# mu is the location, a is the scale
dt_ls <- function(x, df, mu, a) 1/a * dt((x - mu)/a, df)
pt_ls <- function(x, df, mu, a) pt((x - mu)/a, df)
qt_ls <- function(prob, df, mu, a) qt(prob, df)*a + mu
rt_ls <- function(n, df, mu, a) rt(n,df)*a + mu

hazard_func <- function(r, lam) { 
  odds <- 1/lam * rep(1, r) 
  #return (sum(odds[which(odds < 1)]))
  return (1/lam)
}
