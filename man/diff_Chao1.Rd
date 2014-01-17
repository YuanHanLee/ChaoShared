\name{diff_Chao1}
\alias{diff_Chao1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diff_Chao1(ksi, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ksi}{
%%     ~~Describe \code{ksi} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
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
function (ksi, s) 
{
    if (sum(ksi == 0) != 0) {
        ksi[which(ksi == 0)] <- 1
    }
    D12 <- ksi[1]
    f11 <- ksi[2]
    fp1 <- ksi[3]
    f1p <- ksi[4]
    fp2 <- ksi[5]
    f2p <- ksi[6]
    if (s == 1) {
        d <- 1
    }
    else if (s == 2) {
        d <- f1p * fp1/(4 * f2p * fp2)
    }
    else if (s == 3) {
        d <- f11 * f1p/(4 * f2p * fp2) + fp1/fp2
    }
    else if (s == 4) {
        d <- f11 * fp1/(4 * f2p * fp2) + f1p/f2p
    }
    else if (s == 5) {
        d <- -f11 * f1p * fp1/(4 * f2p)/fp2^2 - (fp1/fp2)^2/2
    }
    else if (s == 6) {
        d <- -f11 * f1p * fp1/(4 * fp2)/f2p^2 - (f1p/f2p)^2/2
    }
    else {
        d <- 0
    }
    return(d)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
