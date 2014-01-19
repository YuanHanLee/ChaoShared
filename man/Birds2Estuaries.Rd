\name{Birds2Estuaries}
\alias{Birds2Estuaries}
\docType{data}
\title{
Birds Abundance Data in Two Estuaries
}
\description{
This data was analyzed in details by Chao et al. (2000). A total of 85867 and 59646 observations have been made from two estuaries (which are hereafter referred to as Community I and II respectively). \cr\cr
In these two communities, there were, respectively, 155 and 140 species observed, with 111 of these recorded for both areas (shared species).
}
\usage{data(Birds2Estuaries)}
\format{
  A numerical matrix or a data frame with two columns. The first column represent the species abundances of community 1, second column is the species abundance of community 2.
  \tabular{rlllll}{
      \tab \tab [Com 1]  \tab [Com 2] \cr
      \tab [Sp 1 ] \tab 39  \tab 42 \cr
      \tab [Sp 2 ] \tab 0  \tab 70 \cr
      \tab [Sp 3 ] \tab 0  \tab 1 \cr
      \tab [Sp 4 ] \tab 842  \tab 616 \cr
      \tab [.....] \tab ...  \tab ... \cr
      \tab [Sp200] \tab 63  \tab 8 \cr
      \tab [Sp201] \tab 1  \tab 1 \cr
    }

}
\source{
Chao, A., Hwang, W.-H., Chen, Y.-C. and Kuo. C.-Y. (2000). Estimating the number of shared species in two communities. \emph{Statistica Sinica} \bold{10}, 227-246.
}
\examples{
data(Birds2Estuaries)
}
\keyword{datasets}
