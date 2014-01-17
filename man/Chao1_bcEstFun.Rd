\name{Chao1_bcEstFun}
\alias{Chao1_bcEstFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Chao1_bcEstFun(x1, x2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
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
function (x1, x2) 
{
    D12 <- sum(x1 > 0 & x2 > 0)
    x1_share <- x1[which(x1 > 0 & x2 > 0)]
    x2_share <- x2[which(x1 > 0 & x2 > 0)]
    f11 <- sum(x1_share == 1 & x2_share == 1)
    f1.plus <- sum(x1_share == 1 & x2_share >= 1)
    fplus.1 <- sum(x2_share == 1 & x1_share >= 1)
    f2.plus <- sum(x1_share == 2 & x2_share >= 1)
    fplus.2 <- sum(x2_share == 2 & x1_share >= 1)
    est <- D12 + f11 * f1.plus * fplus.1/(4 * (f2.plus + 1) * 
        (fplus.2 + 1)) + f1.plus * (f1.plus - 1)/(2 * (f2.plus + 
        1)) + fplus.1 * (fplus.1 - 1)/(2 * (fplus.2 + 1))
    return(est)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
