\name{BasicFun.Sam}
\alias{BasicFun.Sam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BasicFun.Sam(y1, y2, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
%%     ~~Describe \code{y1} here~~
}
  \item{y2}{
%%     ~~Describe \code{y2} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
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
function (y1, y2, B) 
{
    t1 <- y1[1]
    t2 <- y2[1]
    x1 <- y1[-1]
    x2 <- y2[-1]
    D1 <- sum(x1 > 0)
    D2 <- sum(x2 > 0)
    D12 <- sum(x1 > 0 & x2 > 0)
    Q11 <- sum(x1 == 1 & x2 == 1)
    Q1.plus <- sum(x1 == 1 & x2 >= 1)
    Qplus.1 <- sum(x2 == 1 & x1 >= 1)
    Q2.plus <- sum(x1 == 2 & x2 >= 1)
    Qplus.2 <- sum(x2 == 2 & x1 >= 1)
    out <- list(t1 = t1, t2 = t2, D1 = D1, D2 = D2, D12 = D12, 
        B = B, Q11 = Q11, Q1.plus = Q1.plus, Q1.plus = Q1.plus, 
        Qplus.1 = Qplus.1, Q2.plus = Q2.plus, Qplus.2 = Qplus.2)
    class(out) <- "BasicFun.Sam"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
