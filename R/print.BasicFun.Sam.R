print.BasicFun.Sam <-
function(x, ...) {
  cat("(1)  BASIC DATA INFORMATION:", "\n\n")
  cat("                      (Number of samples from community 1)   t1  = ", x$t1, "\n")
  cat("                      (Number of samples from community 2)   t2  = ", x$t2, "\n")
  cat("               (Number of observed species in community 1)   D1  = ", x$D1, "\n")
  cat("               (Number of observed species in community 2)   D2  = ", x$D2, "\n")
  cat("     (Number of observed shared species in two communities)  D12 = ", x$D12, "\n")
  cat("              (Bootstrap replications for s.e. estimate)    ", x$B,   "\n\n")
  cat("     Some Statistics:", "\n")
  cat("         --------------------------------------------------------------", "\n")
  cat("         Q[11] =", x$Q11, "; ", "Q[1+]", x$Q1.plus, ";", "Q[+1] =", x$Qplus.1, "; ", "Q[2+] =", x$Q2.plus, "; ", "Q[+2] =", x$Qplus.2, "\n")
  cat("         --------------------------------------------------------------", "\n")
}
