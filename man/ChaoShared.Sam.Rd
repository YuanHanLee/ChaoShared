\name{ChaoShared.Sam}
\alias{ChaoShared.Sam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ChaoShared.Sam(y1, y2, method = c("all", "Chao2-shared", "Chao2-shared-bc", "Lower-bound", "Lower-bound-bc"), conf = 0.95, se = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
%%     ~~Describe \code{y1} here~~
}
  \item{y2}{
%%     ~~Describe \code{y2} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
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
function (y1, y2, method = c("all", "Chao2-shared", "Chao2-shared-bc", 
    "Lower-bound", "Lower-bound-bc"), conf = 0.95, se = TRUE) 
{
    method <- match.arg(method)
    if (method == "all") {
        a <- Chao2_sharedFun(y1, y2, conf)
        b <- Chao2_bcFun(y1, y2, conf)
        c <- PanFun.Sam(y1, y2, conf)
        d <- PanbcFun.Sam(y1, y2, conf)
        out <- rbind(a, b, c, d)
    }
    if (method == "Chao2-shared") 
        out <- Chao2_sharedFun(y1, y2, conf)
    if (method == "Chao2-shared-bc") 
        out <- Chao2_bcFun(y1, y2, conf)
    if (method == "Lower-bound") 
        out <- PanFun.Sam(y1, y2, conf)
    if (method == "Lower-bound-bc") 
        out <- PanbcFun.Sam(y1, y2, conf)
    if (se == FALSE) {
        out <- data.frame(Estimator = out[, 1], row.names = rownames(out))
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
