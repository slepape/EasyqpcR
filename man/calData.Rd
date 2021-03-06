\name{calData}
\alias{calData}

\title{
Calculation of calibration factors
}
\description{
This function determines the calibration factors (CF) using the method described
 in Hellemans et al. (2007).
}
\usage{
calData(data)
}

\arguments{
  \item{data}{
data.frame containing the NRQs of your calibrator(s) for each gene obtained by
 the nrmData function of this package.
}
}
\details{
This function is necessary for comparing different quantitative real-time PCR
 runs to reduce the inter-run variability (Hellemans et al. (2007)). Then, the
 results obtained have to be included in an R object and then be inputed in the
 nrmData function (see the vignette for more informations).
}
\value{
This function returns the calibration factor associated to each gene for the
 whole runs.
}
\references{
Jan Hellemans, Geert Mortier, Anne De Paepe, Frank Speleman and Jo Vandesompele.
 qBase relative quantification framework and software for management and
 automated analysis of real-time quantitative PCR data. Genome Biology 2007,
 8:R19 (doi:10.1186/gb-2007-8-2-r19). <url:http://genomebiology.com/2007/8/2/R19>


}
\author{
Sylvain Le Pape (IRTOMIT-INSERM U1082) <sylvain.le.pape@univ-poitiers.fr>
}

\examples{
data(qPCR_run1,qPCR_run2,qPCR_run3)

nrmData(data = qPCR_run1 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)

nrmData(data = qPCR_run2 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)

nrmData(data = qPCR_run3 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)


## Isolate the calibrator NRQ values of the first biological replicate

a <- nrmData(data = qPCR_run1 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[[3]] 

## Isolate the calibrator NRQ values of the first biological replicate

b <- nrmData(data = qPCR_run2 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[[3]]
 
## Isolate the calibrator NRQ values of the first biological replicate

c <- nrmData(data = qPCR_run3 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[[3]]

## Regrouping the calibrator NRQ values of all the biological replicates

d <- rbind(a, b, c) 

## Calibration factor calculation

e <- calData(d) 
}

\keyword{Inter-run calibration}
\keyword{Calibration factors}
