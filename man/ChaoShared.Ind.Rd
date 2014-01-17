\name{ChaoShared.Ind}
\alias{ChaoShared.Ind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ChaoShared.Ind(x1, x2, method = c("all", "Homogeneous", "Heterogeneous(ACE-shared)", "Chao1(shared)", "Chao1-shared-bc", "Lower-bound", "Lower-bound-bc"), B = 200, conf = 0.95, se = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
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
function (x1, x2, method = c("all", "Homogeneous", "Heterogeneous(ACE-shared)", 
    "Chao1(shared)", "Chao1-shared-bc", "Lower-bound", "Lower-bound-bc"), 
    B = 200, conf = 0.95, se = TRUE) 
{
    method <- match.arg(method)
    if (se == FALSE) 
        B <- 1
    if (method == "all") {
        a <- HomoFun(x1, x2, B, conf)
        b <- HeteroFun(x1, x2, B, conf)
        c <- Chao1_sharedFun(x1, x2, conf)
        d <- Chao1_bcFun(x1, x2, conf)
        e <- PanFun(x1, x2, conf)
        f <- PanbcFun(x1, x2, conf)
        out <- rbind(a, b, c, d, e, f)
    }
    if (method == "Homogeneous") 
        out <- HomoFun(x1, x2, B, conf)
    if (method == "Heterogeneous(ACE-shared)") 
        out <- HeteroFun(x1, x2, B, conf)
    if (method == "Chao1(shared)") 
        out <- Chao1_sharedFun(x1, x2, conf)
    if (method == "Chao1-shared-bc") 
        out <- Chao1_bcFun(x1, x2, conf)
    if (method == "Lower-bound)") 
        out <- PanFun(x1, x2, conf)
    if (method == "Lower-bound-bc") 
        out <- PanbcFun(x1, x2, conf)
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
