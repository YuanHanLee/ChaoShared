\name{BBR}
\alias{BBR}
\docType{data}
\title{
The Hong Kong Big Bird Race (BBR)
}
\description{
In 1999, a total of 217 species was observed by 19 teams competed in Hong Kong BBR. In 2000, a total of 220 species was observed. Merging the two-year data by species names, we found that there were 115 observed shared species. \cr\cr
In 1999, the winning team recorded 152 species. This means that the winning team missed 65 species that were observed by at least one of the other teams. In 2000, the winning team recorded 154 species; thus, 66 species that were observed by the other teams were missed by the winning team.
}
\usage{data(BBR)}
\format{
  A numerical matrix or a data frame with two columns. \cr
  In each column "the first entry is the total number of sampling units (t)", and followed by the incidences frequencies. The format is: \cr
  \tabular{rlllll}{
      \tab           \tab [Com 1]  \tab [Com 2] \cr
      \tab [t ]      \tab 19  \tab 20 \cr
      \tab [Sp1]     \tab 0   \tab 0 \cr
      \tab [Sp2]     \tab 0   \tab 0 \cr
      \tab [...]     \tab ..  \tab .. \cr
      \tab [Sp157]   \tab 8   \tab 1 \cr
      \tab [...]     \tab ..  \tab .. \cr
      \tab [Sp358]   \tab 19  \tab 20 \cr
    }
}
\source{
Chao, A., Shen, T.-J. and Hwang, W. H. (2006). Application of Laplace's boundary-mode approximations to estimate species and shared species richness. \emph{Australian and New Zealand Journal of Statistics} \bold{48}, 117-128.
}
\examples{
data(BBR)
}
\keyword{datasets}
