init_model <- function(a1, a2, b2, c1, c2, d2, k, history = NULL) {
  set_K(k)
  prior <- new(Prior, a1, a2, b2, c1, c2, d2)
  if (is.null(history)) {
    new(Model, prior, k, 0, 0)
  } else {
    new(Model, prior, k, count_cookie_history(history), count_hostname_history(history))
  }
}