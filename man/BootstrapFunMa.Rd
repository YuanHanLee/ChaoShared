\name{BootstrapFunMa}
\alias{BootstrapFunMa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BootstrapFunMa(x1, x2, B, FunName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
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
function (x1, x2, B, FunName) 
{
    n1 <- sum(x1)
    n2 <- sum(x2)
    z <- SortDataFun(x1, x2)
    z1 <- z[, 1]
    z2 <- z[, 2]
    newprob <- ExtenProbFunMa(z1, z2)
    p1 <- newprob$prob1
    p2 <- newprob$prob2
    set.seed(123)
    X1 <- rmultinom(B, n1, p1)
    set.seed(123)
    X2 <- rmultinom(B, n2, p2)
    X <- rbind(X1, X2)
    se <- sd(apply(X, 2, function(x) {
        y <- matrix(x, ncol = 2)
        y1 <- y[, 1]
        y2 <- y[, 2]
        FunName(y1, y2)
    }), na.rm = T)
    return(se)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
