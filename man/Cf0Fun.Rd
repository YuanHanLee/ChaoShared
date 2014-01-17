\name{Cf0Fun}
\alias{Cf0Fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Cf0Fun(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x) 
{
    n <- sum(x)
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if (f2 > 0) {
        C <- 1 - f1/n * ((n - 1) * f1/((n - 1) * f1 + 2 * f2))
        f0 <- (n - 1)/n * f1^2/(2 * f2)
    }
    else if (f2 == 0 & f1 != 0) {
        C <- 1 - f1/n * ((n - 1) * (f1 - 1)/((n - 1) * (f1 - 
            1) + 2))
        f0 <- (n - 1)/n * f1 * (f1 - 1)/2
    }
    else {
        C <- 1
        f0 <- (n - 1)/n * f1 * (f1 - 1)/2
    }
    f0 <- ceiling(f0)
    return(c(C, f0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
