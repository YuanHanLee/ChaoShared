ChaoShared <-
function(data, datatype = c("abundance", "incidence"), 
                       se = TRUE, nboot = 200, conf = 0.95) {
  
  method <- "all"
  if (se == TRUE) {
    if (nboot < 1)
      nboot <- 1
    if (nboot == 1)
      cat("Warning: When \"nboot\" =" ,nboot, ", the bootstrap s.e. and confidence interval can't be calculated.", 
          "\n\n")  
  }
  
  if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) {
    cat("Warning: \"conf\"(confidence level) must be a numerical value between 0 and 1, e.g. 0.95.",
        "\n")
    cat("          We use \"conf\" = 0.95 to calculate!", 
        "\n\n")
    conf <- 0.95
  }
  
  datatype <- match.arg(datatype)
  if (datatype == "abundance") {
    x1 <- data[, 1]
    x2 <- data[, 2]
    Basic <- BasicFun(x1, x2, nboot)
    #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
    output <- ChaoShared.Ind(x1, x2, method, nboot, conf, se)
  }
  if (datatype == "incidence") {
    y1 <- data[, 1]
    y2 <- data[, 2]
    Basic <- BasicFun.Sam(y1, y2, B=nboot)
    #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
    output <- ChaoShared.Sam(y1, y2, method, conf, se)
  }
  out <- list(BASIC_DATA_INFORMATION=Basic, 
              ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES=output)
  return(out)
}
