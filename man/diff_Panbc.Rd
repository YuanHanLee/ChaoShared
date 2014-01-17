\name{diff_Panbc}
\alias{diff_Panbc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diff_Panbc(ksi, s)
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
    D12 <- ksi[1]
    f11 <- ksi[2]
    fp1 <- ksi[3]
    f1p <- ksi[4]
    fp2 <- ksi[5]
    f2p <- ksi[6]
    f22 <- ksi[7]
    K1 <- ksi[8]
    K2 <- ksi[9]
    if (s == 1) {
        d <- 1
    }
    else if (s == 2) {
        d <- K1 * K2 * (2 * f11 - 1)/(4 * (f22 + 1))
    }
    else if (s == 3) {
        d <- K2 * (2 * fp1 - 1)/(2 * (fp2 + 1))
    }
    else if (s == 4) {
        d <- K1 * (2 * f1p - 1)/(2 * (f2p + 1))
    }
    else if (s == 5) {
        d <- -K2 * fp1 * (fp1 - 1)/(fp2 + 1)^2/2
    }
    else if (s == 6) {
        d <- -K1 * f1p * (f1p - 1)/(f2p + 1)^2/2
    }
    else {
        d <- -K1 * K2 * f11 * (f11 - 1)/4/(f22 + 1)^2
    }
    return(d)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
