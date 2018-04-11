\name{Efficiency_calculation}
\alias{Efficiency_calculation}
\docType{data}
\title{
Raw data for primer amplification efficiency calculation.
}
\description{
This is a dataset containg the raw data (Cq values) of a qPCR run for the primer
 amplification efficiency calculation.
}
\usage{data(Efficiency_calculation)}
\format{
  A data frame with 15 observations on the following 3 variables.
  \describe{
    \item{\code{Samples}}{5 cDNA dilutions (1/1, 1/10, 1/100, 1/1000, 1/10000)}
    \item{\code{Gene.1}}{The first gene}
    \item{\code{Gene.2}}{The second gene}
  }
}
\details{
This data.frame is composed by 2 genes (Gene.1 and Gene.2). There are 5 cDNA
 dilutions (1/1, 1/10, 1/100, 1/1000, 1/10000).
}
\source{
S. Le Pape, IRTOMIT-INSERM U1082.
}
\references{
S. Le Pape, IRTOMIT-INSERM U1082.
}
\examples{
data(Efficiency_calculation)
}
\keyword{datasets}
