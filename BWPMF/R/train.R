init_model <- function(a1, a2, b2, c1, c2, d2, k, history) {
  set_K(k)
  prior <- new(Prior, a1, a2, b2, c1, c2, d2)
  new(Model, prior, k, count_cookie_history(history), count_hostname_history(history))
}

#'@title Train a Poisson Matrix Factorization Model
#'@param src path. Please see details for more information.
#'@param prior list. A named numeric vector of prior.
#'@param output path to a directory. If the path does not exist, R will try
#'  to create the directory. The meta data and model will be written to the 
#'  directory.
#'@details
#'The Poisson Matrix Factorization(PMF) model assumes that the counting response
#'\eqn{y_{u,i}} corresponding to the user $u$ and item $i$ is poisson distributed 
#'with \eqn{Ey_{u,i} = \sum_{k=1}^K {\theta_{u,k} \beta_{i,k}}}. 
#'
#'
#'@export
train_pmf <- function(src, prior, output) {
  
}