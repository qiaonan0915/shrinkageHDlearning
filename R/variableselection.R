#library(MASS)
require(dlsa)
#library(doParallel)

linear_fit <- function(X, Y, init.beta = NULL, sd_random. = NULL, sigma_e. = NULL){ #, sigma_ind. = NULL
 #sd_random. <- sd_random; sigma_e. <- sigma_e
  N = nrow(X)
  p = ncol(X)
  Z_k <- rep(1, nrow(X))
  H_inv <-  diag(1,nrow = nrow(X)) - Z_k%*%t(Z_k)*(sd_random./sigma_e.)/(1 + nrow(X)* (sd_random./sigma_e.) )
  hess <- t(X) %*%H_inv%*% X/sigma_e. ##这样算也对
 # gamma_inv <- solve((sd_random./sigma_e.)* Z_k%*%t(Z_k) + diag(1,nrow = nrow(X)))
  # hess <- sigma_ind. * t(X) %*% gamma_inv %*% X ##这样算也对,但需要sigma_ind.
  return(list(theta = init.beta, Sig_inv = hess))
}


variable_selection <- function(X, Y, XY, K,init.beta, sd_random., sigma_e.){
  reults_dlsa <- dlsa.fit(Y~., XY, K, lasso_fun = 1, fit_fun = linear_fit, init.beta = beta_lmm, sd_random. = sd_random, sigma_e. = sigma_e )#$result$theta
  beta_dlsa <- reults_dlsa$lsa_result$beta.bic
  return(beta_dlsa)
}
