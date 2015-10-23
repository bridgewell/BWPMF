library(BWPMF)

check <- function(a) {
  a2 <- test_phi_on_disk(tmp.path <- tempfile(), a)
  stopifnot(abs(a - a2$retval1) < 1e-5)
  stopifnot(abs(a - a2$retval2 + 1) < 1e-5)
  stopifnot(abs(a2$retval1 - a2$retval2 + 1) < 1e-5)
}

check(a <- matrix(rnorm(50), 5, 10))
check(a <- matrix(rnorm(50), 10, 10))
check(a <- matrix(rnorm(50), 50, 10))

