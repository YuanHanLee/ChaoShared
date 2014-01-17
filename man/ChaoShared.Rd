\name{ChaoShared}
\alias{ChaoShared}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ChaoShared(data, datatype = c("abundance", "incidence"), se = TRUE, nboot = 200, conf = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{datatype}{
%%     ~~Describe \code{datatype} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
}
  \item{nboot}{
%%     ~~Describe \code{nboot} here~~
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, datatype = c("abundance", "incidence"), se = TRUE, 
    nboot = 200, conf = 0.95) 
{
    method <- "all"
    if (se == TRUE) {
        if (nboot < 1) 
            nboot <- 1
        if (nboot == 1) 
            cat("Warning: When \"nboot\" =", nboot, ", the bootstrap s.e. and confidence interval can't be calculated.", 
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
        output <- ChaoShared.Ind(x1, x2, method, nboot, conf, 
            se)
    }
    if (datatype == "incidence") {
        y1 <- data[, 1]
        y2 <- data[, 2]
        Basic <- BasicFun.Sam(y1, y2, B = nboot)
        output <- ChaoShared.Sam(y1, y2, method, conf, se)
    }
    out <- list(BASIC_DATA_INFORMATION = Basic, ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES = output)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
