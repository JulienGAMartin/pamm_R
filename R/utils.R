ci.p <- function(x, confidence = 0.95, alpha = 1 - confidence, na.rm = TRUE, ...) {

  est <- mean(x, na.rm = na.rm)
  if(nobs.p(x) <= 1) {
    stderr <- NA
    ci.low <- NA
    ci.high <- NA
  } else {
    stderr <- sd(x, na.rm = na.rm)/sqrt(nobs.p(x))
    ci.low <- est + qt(alpha/2, nobs.p(x) - 1) * stderr
    ci.high <- est - qt(alpha/2, nobs.p(x) - 1) * stderr
  }
  retval <- c(Estimate = est, `CI lower` = ci.low, `CI upper` = ci.high, `Std. Error` = stderr)

  retval
}

nobs.p <- function(object) {
  sum(!is.na(object))
}
