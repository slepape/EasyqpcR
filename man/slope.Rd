\name{slope}
\alias{slope}

\title{
Function to calculate the amplification efficiency
}
\description{
This function calculates the amplification efficiency from classical qPCR
 dilution experiment using the Cq values. The Cq values are plotted against the
 logarithmized concentration (or dilution) values, a linear regression line is
 fit and the efficiency calculated by E = 10^{-1/slope}. 
}
\usage{
slope(data, q, r, na.rm = F)
}

\arguments{
  \item{data}{
data.frame containing row datas (genes in columns, samples in rows, Cq values).
}
  \item{q}{
numeric, cDNA dilution values.
}
  \item{r}{
numeric, number of qPCR replicates.
}
  \item{na.rm}{
logical, indicating whether NA values should be stripped before the computation
 proceeds.
}
}

\value{
 \item{Efficiency}{Primer amplification efficiency.}
 \item{Slope}{Slope of the dilution curve.}
 \item{Intercept}{Intercept of the dilution curve}
}

\author{
Sylvain Le pape <sylvain.le.pape@univ-poitiers.fr>
}

\seealso{
You can also see the \code{qpcR} package with the \code{calib} and \code{calib2}
 functions.
}
\examples{
data(Efficiency_calculation)

slope(data = Efficiency_calculation, q=c(1000, 100 ,10, 1, 0.1), r=3, na.rm=TRUE)
}

\keyword{Efficiency}
