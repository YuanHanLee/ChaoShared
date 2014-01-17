\name{SortDataFun}
\alias{SortDataFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SortDataFun(x1, x2)
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
    D1 <- sum(x1 > 0)
    D2 <- sum(x2 > 0)
    D12 <- sum(x1 > 0 & x2 > 0)
    sit <- order(x1 > 0 & x2 > 0, decreasing = T)
    common <- cbind(x1[sit[1:D12]], x2[sit[1:D12]])
    sit1 <- order(x1 > 0 & x2 == 0, decreasing = T)
    set.seed(123)
    sit1 <- sample(sit1[1:(D1 - D12)], replace = F)
    special1 <- cbind(x1[sit1[1:(D1 - D12)]], x2[sit1[1:(D1 - 
        D12)]])
    sit2 <- order(x1 == 0 & x2 > 0, decreasing = T)
    set.seed(123)
    sit2 <- sample(sit2[1:(D2 - D12)], replace = F)
    special2 <- cbind(x1[sit2[1:(D2 - D12)]], x2[sit2[1:(D2 - 
        D12)]])
    z <- rbind(common, special1, special2)
    return(z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
