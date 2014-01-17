\name{Chao2_sharedFun}
\alias{Chao2_sharedFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Chao2_sharedFun(y1, y2, conf = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
%%     ~~Describe \code{y1} here~~
}
  \item{y2}{
%%     ~~Describe \code{y2} here~~
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
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
function (y1, y2, conf = 0.95) 
{
    x1 <- y1[-1]
    x2 <- y2[-1]
    f2p <- sum(x1 == 2 & x2 >= 1)
    fp2 <- sum(x1 >= 1 & x2 == 2)
    est <- Chao2_sharedEstFun(y1, y2)
    if (f2p == 0 || fp2 == 0) {
        se <- VarEstFun.Sam(y1, y2, diffFun = diff_Chao2bc, FunName = Chao2_bcEstFun)
    }
    else {
        se <- VarEstFun.Sam(y1, y2, diffFun = diff_Chao2, FunName = Chao2_sharedEstFun)
    }
    CI <- logCI(y1, y2, est, se, conf)
    out <- matrix(c(est, se, CI), nrow = 1)
    rownames(out) <- c("Chao2-shared")
    colnames(out) <- c("Estimator", "Est_s.e.", paste(conf * 
        100, "\% Lower"), paste(conf * 100, "\% Upper"))
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
