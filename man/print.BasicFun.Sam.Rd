\name{print.BasicFun.Sam}
\alias{print.BasicFun.Sam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.BasicFun.Sam(x)
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
    cat("(1)  BASIC DATA INFORMATION:", "\n\n")
    cat("                      (Number of samples from community 1)   t1  = ", 
        x$t1, "\n")
    cat("                      (Number of samples from community 2)   t2  = ", 
        x$t2, "\n")
    cat("               (Number of observed species in community 1)   D1  = ", 
        x$D1, "\n")
    cat("               (Number of observed species in community 2)   D2  = ", 
        x$D2, "\n")
    cat("     (Number of observed shared species in two communities)  D12 = ", 
        x$D12, "\n")
    cat("              (Bootstrap replications for s.e. estimate)    ", 
        x$B, "\n\n")
    cat("     Some Statistics:", "\n")
    cat("         --------------------------------------------------------------", 
        "\n")
    cat("         Q[11] =", x$Q11, "; ", "Q[1+]", x$Q1.plus, 
        ";", "Q[+1] =", x$Qplus.1, "; ", "Q[2+] =", x$Q2.plus, 
        "; ", "Q[+2] =", x$Qplus.2, "\n")
    cat("         --------------------------------------------------------------", 
        "\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
