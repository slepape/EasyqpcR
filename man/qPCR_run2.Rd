\name{qPCR_run2}
\alias{qPCR_run2}
\docType{data}
\title{
Raw data from the second qPCR run.
}
\description{
This is a dataset containg the raw data (Cq values) of the second qPCR run (the
 second biological replicate).
}
\usage{data(qPCR_run2)}
\format{
  A data frame with 15 observations on the following 5 variables.
  \describe{
    \item{\code{Samples}}{a factor with levels \code{Calibrator} \code{Control
 1} \code{Control 2} \code{Treatment 1} \code{Treatment 2}}
    \item{\code{RG1}}{a numeric vector}
    \item{\code{RG2}}{a numeric vector}
    \item{\code{TG}}{a numeric vector}
    \item{\code{TGb}}{a numeric vector}
  }
}
\details{
This data.frame is composed by 2 reference genes (RG1 and RG2), and two target
 genes (TG and TGb). There are 5 samples, starting by the controls (1 and 2),
 the treatment groups (1 and 2), and a calibrator (calibrator). It is the second
 biological replicate and there is a qPCR triplicate by sample.
}
\source{
S. Le Pape, IRTOMIT-INSERM U1082.
}
\references{
S. Le Pape, IRTOMIT-INSERM U1082.
}
\examples{
data(qPCR_run2)
}
\keyword{datasets}
