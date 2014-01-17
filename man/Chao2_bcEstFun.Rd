\name{Chao2_bcEstFun}
\alias{Chao2_bcEstFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Chao2_bcEstFun(y1, y2)
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
    t1 <- y1[1]
    t2 <- y2[1]
    c1 <- (t1 - 1)/t1
    c2 <- (t2 - 1)/t2
    x1 <- y1[-1]
    x2 <- y2[-1]
    D12 <- sum(x1 > 0 & x2 > 0)
    Q11 <- sum(x1 == 1 & x2 == 1)
    Q1.plus <- sum(x1 == 1 & x2 >= 1)
    Qplus.1 <- sum(x2 == 1 & x1 >= 1)
    Q2.plus <- sum(x1 == 2 & x2 >= 1)
    Qplus.2 <- sum(x2 == 2 & x1 >= 1)
    est <- D12 + Q11 * c1 * c2 * Q1.plus * Qplus.1/(4 * (Q2.plus + 
        1) * (Qplus.2 + 1)) + c1 * Q1.plus * (Q1.plus - 1)/(2 * 
        (Q2.plus + 1)) + c2 * Qplus.1 * (Qplus.1 - 1)/(2 * (Qplus.2 + 
        1))
    return(est)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
