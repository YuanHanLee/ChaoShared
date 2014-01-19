\name{ChaoShared}
\alias{ChaoShared}
\title{
Estimation of the number of shared species in two communities
}
\description{
\kbd{ChaoShared} is a function to provide a shared species estimator in tow communities based on the following two sampling schemes. \cr
  \enumerate{
    \item{In each community, a random sample of individuals is taken and species frequencies or abundances are recored;}
    \item{Each communtiy is sample several times or the whole area is divided into several quadrats and species presence/absence data for multiple samples/quadrats are recorded.}
  }
}

\usage{
ChaoShared(data, datatype = c("abundance", "incidence"), se = TRUE, nboot = 200, 
           conf = 0.95)
}

\arguments{
  \item{data}{
  a numerical matrix or a data frame with two columns. Each column represent the species abundances or incidence frequencies of each community.
  If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequencies in each column (See examples).
}
  \item{datatype}{
  the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or presence/absence sample-base incidence data (\code{datatype = "incidence"}).

}
  \item{se}{
calculate bootstrap standard error and show confidence interval; default is \code{TRUE}.
}
  \item{nboot}{
the number of bootstrap resampling times, default is \code{200}.
}
  \item{conf}{
  a positive number \eqn{\le} 1. "\code{conf}" specifies the confidence level for confidence interval. The default is \code{0.95}.
}
}
\value{
\kbd{ChaoShared} returns a list of two. First is the basic data informaiton; the other is a table of various estimators, their standard error and \code{100*conf} \% confidence interval.
}

\author{
Y.H. Lee ; K.S. Tseng ; Anne Chao\cr
Maintainer : Anne Chao \email{chao@stat.nthu.edu.tw}
}

\references{
Chao, A., Hwang, W.-H., Chen, Y.-C. and Kuo. C.-Y. (2000). Estimating the number of shared species in two communities. \emph{Statistica Sinica} \bold{10}, 227-246. \cr\cr
Chao, A., Shen, T.-J. and Hwang, W. H. (2006). Application of Laplace's boundary-mode approximations to estimate species and shared species richness. \emph{Australian and New Zealand Journal of Statistics} \bold{48}, 117-128. \cr\cr
Pan, H. Y., Chao, A. and Foissner, W. (2009). A non-parametric lower bound for the number of species shared by multiple communities. Journal of Agricultural, \emph{Biological and Environmental Statistics}. \bold{14}, 452-468.
}

\examples{
# load the individual-base (abundance) data
data(Birds2Estuaries)

# Estimation of shared species
ChaoShared(Birds2Estuaries, datatype="abundance", se=TRUE, nboot=200, conf=0.95)

# load the presence/absence sample-base (incidence) data
data(BBR)

# Estimation of shared species
ChaoShared(BBR, datatype="incidence", se=TRUE, nboot=200, conf=0.95)
}
\keyword{shared species}
