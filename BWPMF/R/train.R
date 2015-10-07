init_model <- function(a1, a2, b2, c1, c2, d2, k, history) {
  prior <- new(Prior, a1, a2, b2, c1, c2, d2)
  new(Model, prior, k, count_cookie_history(history), count_hostname_history(history))
}