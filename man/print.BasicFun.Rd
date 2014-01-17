\name{print.BasicFun}
\alias{print.BasicFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.BasicFun(x)
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
    cat("              (Number of observed individuals in community 1)   n1  = ", 
        x$n1, "\n")
    cat("              (Number of observed individuals in community 2)   n2  = ", 
        x$n2, "\n")
    cat("              (Number of observed species in community 1)       D1  = ", 
        x$D1, "\n")
    cat("              (Number of observed species in community 2)       D2  = ", 
        x$D2, "\n")
    cat("              (Number of observed shared species)               D12 = ", 
        x$D12, "\n")
    cat("              (Bootstrap replications for s.e. estimate)       ", 
        x$B, "\n\n")
    cat("     \"Rare\" Shared Species Group: (Both Frequencies can only up to 10)", 
        "\n")
    cat("         Some Statistics:", "\n")
    cat("         --------------------------------------------------------------", 
        "\n")
    cat("         f[11] =", x$f11, "; ", "f[1+]", x$f1.plus, 
        ";", "f[+1] =", x$fplus.1, "; ", "f[2+] =", x$f2.plus, 
        "; ", "f[+2] =", x$fplus.2, "\n")
    cat("         --------------------------------------------------------------", 
        "\n")
    cat("            (Number of observed individuals in community 1)     n1_rare   = ", 
        x$n1_rare, "\n")
    cat("            (Number of observed individuals in community 2)     n2_rare   = ", 
        x$n2_rare, "\n")
    cat("            (Number of observed shared species)                 D12_rare  = ", 
        x$D12_rare, "\n")
    cat("            (Estimated sample coverage)                         C12_rare  = ", 
        x$C12_rare, "\n")
    cat("            (Estimated CCVs)                                    CCV_1     = ", 
        x$CCV_1, "\n")
    cat("                                                                CCV_2     = ", 
        x$CCV_2, "\n")
    cat("                                                                CCV_12    = ", 
        x$CCV_12, "\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
