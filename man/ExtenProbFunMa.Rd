\name{ExtenProbFunMa}
\alias{ExtenProbFunMa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ExtenProbFunMa(z1, z2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z1}{
%%     ~~Describe \code{z1} here~~
}
  \item{z2}{
%%     ~~Describe \code{z2} here~~
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
function (z1, z2) 
{
    x1 <- z1
    x2 <- z2
    n1 <- sum(x1)
    n2 <- sum(x2)
    D1 <- sum(x1 > 0)
    D2 <- sum(x2 > 0)
    D12 <- sum(x1 > 0 & x2 > 0)
    f11 <- sum(x1 == 1 & x2 == 1)
    f22 <- sum(x1 == 2 & x2 == 2)
    f1p <- sum(x1 == 1 & x2 >= 1)
    fp1 <- sum(x1 >= 1 & x2 == 1)
    f2p <- sum(x1 == 2 & x2 >= 1)
    fp2 <- sum(x1 >= 1 & x2 == 2)
    K1 <- (n1 - 1)/n1
    K2 <- (n2 - 1)/n2
    f0p <- ifelse(f2p == 0, f0p <- ceiling(K1 * f1p * (f1p - 
        1)/2/(f2p + 1)), f0p <- ceiling(K1 * f1p^2/2/f2p))
    fp0 <- ifelse(fp2 == 0, fp0 <- ceiling(K2 * fp1 * (fp1 - 
        1)/2/(fp2 + 1)), fp0 <- ceiling(K2 * fp1^2/2/fp2))
    f00 <- ifelse(f22 == 0, f00 <- ceiling(K1 * K2 * f11 * (f11 - 
        1)/4/(f22 + 1)), f00 <- ceiling(K1 * K2 * f11^2/4/f22))
    ga1 <- D1 - D12
    ga2 <- D2 - D12
    tmp1 <- Cf0Fun(x1)
    Chat1 <- tmp1[1]
    f0_1 <- tmp1[2]
    more1 <- max(f0_1, (max(ga1, fp0) - ga1 + f0p + f00))
    add1 <- max(f0_1 - max(ga1, fp0) + ga1 - f0p - f00, 0)
    lambda1 <- (1 - Chat1)/sum(x1/n1 * (1 - x1/n1)^n1)
    pi1 <- x1/n1 * (1 - lambda1 * (1 - x1/n1)^n1)
    p0_1 <- (1 - Chat1)/more1
    tmp2 <- Cf0Fun(x2)
    Chat2 <- tmp2[1]
    f0_2 <- tmp2[2]
    more2 <- max(f0_2, (max(ga2, f0p) - ga2 + fp0 + f00))
    add2 <- max(f0_2 - max(ga2, f0p) + ga2 - fp0 - f00, 0)
    lambda2 <- (1 - Chat2)/sum(x2/n2 * (1 - x2/n2)^n2)
    pi2 <- x2/n2 * (1 - lambda2 * (1 - x2/n2)^n2)
    p0_2 <- (1 - Chat2)/more2
    prob1 <- c(pi1[pi1 > 0], rep(p0_1, max(fp0 - ga1, 0)), rep(p0_1, 
        f0p), rep(0, max(ga2 - f0p, 0)), rep(p0_1, f00), rep(p0_1, 
        add1), rep(0, add2))
    prob2 <- c(pi2[1:D12], rep(p0_2, fp0), rep(0, max(ga1 - fp0, 
        0)), pi2[(D12 + ga1 + 1):(D12 + ga1 + ga2)], rep(p0_2, 
        max(f0p - ga2, 0)), rep(p0_2, f00), rep(0, add1), rep(p0_2, 
        add2))
    out <- list(prob1 = prob1, prob2 = prob2)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
