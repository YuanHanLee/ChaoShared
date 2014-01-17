\name{HomoEstFun}
\alias{HomoEstFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HomoEstFun(x1, x2)
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
    x1_share <- x1[which(x1 > 0 & x2 > 0)]
    x2_share <- x2[which(x1 > 0 & x2 > 0)]
    f11 <- sum(x1_share == 1 & x2_share == 1)
    D12_rare <- sum(x1_share <= 10 & x2_share <= 10)
    pos_r <- (x1_share <= 10 & x2_share <= 10)
    pos1_r <- (x1_share == 1 & x2_share <= 10)
    pos2_r <- (x2_share == 1 & x1_share <= 10)
    tmp <- sum(x1_share[pos_r] * x2_share[pos_r])
    if (tmp == 0) 
        tmp <- 1
    C12_rare <- 1 - (sum(x2_share[pos1_r]) + sum(x1_share[pos2_r]) - 
        f11)/tmp
    if (C12_rare == 0 || C12_rare > 1) 
        C12_rare <- 1
    est <- D12 - D12_rare + D12_rare/C12_rare
    return(est)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
