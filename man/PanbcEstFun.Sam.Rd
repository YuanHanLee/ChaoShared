\name{PanbcEstFun.Sam}
\alias{PanbcEstFun.Sam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PanbcEstFun.Sam(y1, y2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
%%     ~~Describe \code{y1} here~~
}
  \item{y2}{
%%     ~~Describe \code{y2} here~~
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
function (y1, y2) 
{
    n1 <- y1[1]
    n2 <- y2[1]
    x1 <- y1[-1]
    x2 <- y2[-1]
    D12 <- sum(x1 > 0 & x2 > 0)
    f11 <- sum(x1 == 1 & x2 == 1)
    f22 <- sum(x1 == 2 & x2 == 2)
    f1p <- sum(x1 == 1 & x2 >= 1)
    fp1 <- sum(x1 >= 1 & x2 == 1)
    f2p <- sum(x1 == 2 & x2 >= 1)
    fp2 <- sum(x1 >= 1 & x2 == 2)
    K1 <- (n1 - 1)/n1
    K2 <- (n2 - 1)/n2
    est <- D12 + K1 * f1p * (f1p - 1)/2/(f2p + 1) + K2 * fp1 * 
        (fp1 - 1)/2/(fp2 + 1) + K1 * K2 * f11 * (f11 - 1)/4/(f22 + 
        1)
    return(est)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
