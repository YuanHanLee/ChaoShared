\name{BasicFun}
\alias{BasicFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BasicFun(x1, x2, B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
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
function (x1, x2, B) 
{
    n1 <- sum(x1)
    n2 <- sum(x2)
    D1 <- sum(x1 > 0)
    D2 <- sum(x2 > 0)
    D12 <- sum(x1 > 0 & x2 > 0)
    x1_share <- x1[which(x1 > 0 & x2 > 0)]
    x2_share <- x2[which(x1 > 0 & x2 > 0)]
    f11 <- sum(x1_share == 1 & x2_share == 1)
    f1.plus <- sum(x1_share == 1 & x2_share <= 10)
    fplus.1 <- sum(x2_share == 1 & x1_share <= 10)
    f2.plus <- sum(x1_share == 2 & x2_share <= 10)
    fplus.2 <- sum(x2_share == 2 & x1_share <= 10)
    D12_rare <- sum(x1_share <= 10 & x2_share <= 10)
    pos <- (x1 > 0 & x2 > 0) & (x1 > 10 | x2 > 10)
    n1_rare <- n1 - sum(x1[pos])
    n2_rare <- n2 - sum(x2[pos])
    pos_r <- (x1_share <= 10 & x2_share <= 10)
    pos1_r <- (x1_share == 1 & x2_share <= 10)
    pos2_r <- (x2_share == 1 & x1_share <= 10)
    tmp <- sum(x1_share[pos_r] * x2_share[pos_r])
    C12_rare <- 1 - (sum(x2_share[pos1_r]) + sum(x1_share[pos2_r]) - 
        f11)/tmp
    T10 <- sum(x1_share[x1_share <= 10 & x2_share <= 10])
    T01 <- sum(x2_share[x1_share <= 10 & x2_share <= 10])
    T11 <- tmp
    T21 <- sum(x1_share[pos_r] * (x1_share - 1)[pos_r] * x2_share[pos_r])
    T12 <- sum(x1_share[pos_r] * (x2_share - 1)[pos_r] * x2_share[pos_r])
    T22 <- sum(x1_share[pos_r] * x2_share[pos_r] * (x1_share - 
        1)[pos_r] * (x2_share - 1)[pos_r])
    S12_0 <- D12_rare/C12_rare
    CCV_1 <- S12_0 * n1_rare * T21/(n1_rare - 1)/T10/T11 - 1
    CCV_2 <- S12_0 * n2_rare * T12/(n2_rare - 1)/T01/T11 - 1
    CCV_12 <- n1_rare * n2_rare * S12_0^2 * T22/((n1_rare - 1) * 
        (n2_rare - 1) * T10 * T01 * T11) - S12_0 * T11/T10/T01 - 
        CCV_1 - CCV_2
    out <- list(n1 = n1, n2 = n2, D1 = D1, D2 = D2, D12 = D12, 
        B = B, f11 = f11, f1.plus = f1.plus, f1.plus = f1.plus, 
        fplus.1 = fplus.1, f2.plus = f2.plus, fplus.2 = fplus.2, 
        n1_rare = n1_rare, n2_rare = n2_rare, D12_rare = D12_rare, 
        C12_rare = C12_rare, CCV_1 = CCV_1, CCV_2 = CCV_2, CCV_12 = CCV_12)
    class(out) <- "BasicFun"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
