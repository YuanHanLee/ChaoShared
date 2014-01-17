\name{VarEstFun.Sam}
\alias{VarEstFun.Sam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VarEstFun.Sam(y1, y2, diffFun, FunName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
%%     ~~Describe \code{y1} here~~
}
  \item{y2}{
%%     ~~Describe \code{y2} here~~
}
  \item{diffFun}{
%%     ~~Describe \code{diffFun} here~~
}
  \item{FunName}{
%%     ~~Describe \code{FunName} here~~
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
function (y1, y2, diffFun, FunName) 
{
    n1 <- y1[1]
    n2 <- y2[1]
    x1 <- y1[-1]
    x2 <- y2[-1]
    D12 <- sum(x1 > 0 & x2 > 0)
    f11 <- sum(x1 == 1 & x2 == 1)
    f22 <- sum(x1 == 2 & x2 == 2)
    f12 <- sum(x1 == 1 & x2 == 2)
    f21 <- sum(x1 == 2 & x2 == 1)
    f1p <- sum(x1 == 1 & x2 >= 1)
    fp1 <- sum(x1 >= 1 & x2 == 1)
    f2p <- sum(x1 == 2 & x2 >= 1)
    fp2 <- sum(x1 >= 1 & x2 == 2)
    K1 <- (n1 - 1)/n1
    K2 <- (n2 - 1)/n2
    ksi <- c(D12, f11, fp1, f1p, fp2, f2p, f22, K1, K2)
    S12 <- FunName(y1, y2)
    CovFun <- function(i, j) {
        if (i == j) {
            cov <- ksi[i] - ksi[i]^2/S12
        }
        else if ((i == 1 & j == 2) | (j == 1 & i == 2)) {
            cov <- f11 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 1 & j == 3) | (j == 1 & i == 3)) {
            cov <- fp1 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 1 & j == 4) | (j == 1 & i == 4)) {
            cov <- f1p - ksi[i] * ksi[j]/S12
        }
        else if ((i == 1 & j == 5) | (j == 1 & i == 5)) {
            cov <- fp2 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 1 & j == 6) | (j == 1 & i == 6)) {
            cov <- f2p - ksi[i] * ksi[j]/S12
        }
        else if ((i == 2 & j == 3) | (j == 2 & i == 3)) {
            cov <- f11 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 2 & j == 4) | (j == 2 & i == 4)) {
            cov <- f11 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 3 & j == 4) | (j == 3 & i == 4)) {
            cov <- f11 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 3 & j == 6) | (j == 3 & i == 6)) {
            cov <- f21 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 4 & j == 5) | (j == 4 & i == 5)) {
            cov <- f12 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 5 & j == 6) | (j == 5 & i == 6)) {
            cov <- f22 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 1 & j == 7) | (j == 1 & i == 7)) {
            cov <- f22 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 5 & j == 7) | (j == 5 & i == 7)) {
            cov <- f22 - ksi[i] * ksi[j]/S12
        }
        else if ((i == 6 & j == 7) | (j == 6 & i == 7)) {
            cov <- f22 - ksi[i] * ksi[j]/S12
        }
        else {
            cov <- 0 - ksi[i] * ksi[j]/S12
        }
        return(cov)
    }
    diff <- diffFun
    i <- rep(c(1:7), 7)
    j <- rep(c(1:7), each = 7)
    var <- sum(mapply(function(i, j) diff(ksi, i) * diff(ksi, 
        j) * CovFun(i, j), i, j))
    se <- sqrt(var)
    return(se)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
