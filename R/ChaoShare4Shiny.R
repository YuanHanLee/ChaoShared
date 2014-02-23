BasicFun.Sam2 <- function(y1, y2, B) {
  t1 <- y1[1]
  t2 <- y2[1]
  x1 <- y1[-1]
  x2 <- y2[-1]
  D1 <- sum(x1 > 0)
  D2 <- sum(x2 > 0)
  D12 <- sum(x1 > 0 & x2 > 0)  
  
  Q11 <- sum(x1 == 1 & x2 == 1)
  Q1.plus <- sum(x1 == 1 & x2 >= 1)
  Qplus.1 <- sum(x2 == 1 & x1 >= 1)
  Q2.plus <- sum(x1 == 2 & x2 >= 1)
  Qplus.2 <- sum(x2 == 2 & x1 >= 1)
  out1 <- data.frame(Value=paste(c("t1 =", "t2 =", "D1 =", "D2 =", "D12 =", "B ="),
                                 c(t1, t2, D1, D2, D12, B)))
  rownames(out1) <- c("Number of samples from community 1", 
                      "Number of samples from community 2",
                      "Number of observed species in community 1",
                      "Number of observed species in community 2",
                      "Number of observed shared species in two communities",
                      "Bootstrap replications for s.e. estimate")
  
  out2 <- matrix(c(Q11, Q1.plus, Qplus.1, Q2.plus, Qplus.2), nrow=1)
  rownames(out2) <- "Value"
  colnames(out2) <- c("Q[11]", "Q[1+]", "Q[+1]", "Q[2+]", "Q[+2]")
  out2 <- as.data.frame(out2)
  
  out <- list(Basic=out1, Rara=out2)
  return(out)
}

BasicFun2 <- function(x1, x2, B) {
  n1 <- sum(x1)
  n2 <- sum(x2)
  D1 <- sum(x1 > 0)
  D2 <- sum(x2 > 0)
  D12 <- sum(x1 > 0 & x2 > 0)
  x1_share <- x1[which(x1 > 0 & x2 > 0)]
  x2_share <- x2[which(x1 > 0 & x2 > 0)]
  f11 <- sum(x1_share == 1 & x2_share == 1)
  f1.plus <- sum(x1_share == 1 & x2_share <= 10)
  fplus.1 <- sum(x2_share == 1 & x1_share <= 10)
  f2.plus <- sum(x1_share == 2 & x2_share <= 10)
  fplus.2 <- sum(x2_share == 2 & x1_share <= 10)
  D12_rare <- sum(x1_share <= 10 & x2_share <= 10)
  
  pos <- (x1 > 0 & x2 > 0) & (x1 > 10 | x2 > 10)
  n1_rare <- n1 - sum(x1[pos])
  n2_rare <- n2 - sum(x2[pos])
  
  pos_r <- (x1_share <= 10 & x2_share <= 10)
  pos1_r <- (x1_share == 1 & x2_share <= 10)
  pos2_r <- (x2_share == 1 & x1_share <= 10)
  
  tmp <- sum(x1_share[pos_r] * x2_share[pos_r])
  C12_rare <- 1 - (sum(x2_share[pos1_r]) + sum(x1_share[pos2_r]) - f11) / tmp
  #   C12_rare <- round(C12_rare, 4)
  
  T10 <- sum(x1_share[x1_share <= 10 & x2_share <= 10])
  T01 <- sum(x2_share[x1_share <= 10 & x2_share <= 10])
  T11 <- tmp
  T21 <- sum(x1_share[pos_r] * (x1_share - 1)[pos_r] * x2_share[pos_r])
  T12 <- sum(x1_share[pos_r] * (x2_share - 1)[pos_r] * x2_share[pos_r])
  
  T22 <- sum(x1_share[pos_r] * x2_share[pos_r] * 
               (x1_share - 1)[pos_r] * (x2_share - 1)[pos_r])
  
  S12_0 <- D12_rare / C12_rare
  CCV_1 <- S12_0 * n1_rare * T21 / (n1_rare - 1) / T10 / T11 - 1
  CCV_2 <- S12_0 * n2_rare * T12 / (n2_rare - 1) / T01 / T11 - 1
  CCV_12 <- n1_rare * n2_rare * S12_0^2 * T22 / 
    ((n1_rare - 1) * (n2_rare - 1) * T10 * T01 * T11) - 
    S12_0 * T11 / T10 / T01 - CCV_1 - CCV_2

  out1 <- data.frame(Value=paste(c("n1 =", "n2 =", "D1 =", "D2 =", "D12 =", "B"),
                                 c(n1, n2, D1, D2, D12, B)))
  rownames(out1) <- c("Number of observed individuals in community 1", 
                     "Number of observed individuals in community 2",
                     "Number of observed species in community 1",
                     "Number of observed species in community 2",
                     "Number of observed shared species", 
                     "Bootstrap replications for s.e. estimate")
  
  out2 <- matrix(c(f11, f1.plus, fplus.1, f2.plus, fplus.2), nrow=1)
  rownames(out2) <- "Value"
  colnames(out2) <- c("f[11]", "f[1+]", "f[+1]", "f[2+]", "f[+2]")
  out2 <- as.data.frame(out2)
  
  out3 <- paste(c("n1_rare =", "n2_rare =", "D12_rare =", "C12_rare =", "CCV_1 =", "CCV_2 =", "CCV_12 ="),
                c(n1_rare, n2_rare, D12_rare, C12_rare, CCV_1, CCV_2, CCV_12))
  out3 <- data.frame(Value=out3)
  rownames(out3) <- c("Number of observed individuals for rare species in community 1", 
                      "Number of observed individuals for rare species in community 2", 
                      "Number of observed shared species for rare species", 
                      "Estimated sample coverage for rare species",
                      "Estimated CCV_1", 
                      "Estimated CCV_2", 
                      "Estimated CCV_12")
  out <- list(Basic=out1, Rare1=out2, Rare2=out3)
  return(out)
}

ChaoShared.Basic <- function(data, datatype = c("abundance", "incidence"), nboot = 200) {
  
  if (nboot < 1)
    nboot <- 1
  if (nboot == 1)
    cat("Warning: When \"nboot\" =" ,nboot, ", the bootstrap s.e. and confidence interval can't be calculated.", 
        "\n\n")  
  
  datatype <- match.arg(datatype)
  if (datatype == "abundance") {
    x1 <- data[, 1]
    x2 <- data[, 2]
    Basic <- BasicFun2(x1, x2, nboot)
  }
  if (datatype == "incidence") {
    y1 <- data[, 1]
    y2 <- data[, 2]
    Basic <- BasicFun.Sam2(y1, y2, B=nboot)
  }
  out <- Basic
  return(out)
}