logCI <-
function(x1, x2, est, se, conf) {
  D12 <- sum(x1 > 0 & x2 > 0)
  t <- est - D12
  z <- qnorm((1-conf)/2, lower.tail=F)
  K <- exp(z * sqrt(log(1 + se^2 / t^2)))
  CI <- c(D12 + t / K, D12 + t * K)
  return(CI)
}
