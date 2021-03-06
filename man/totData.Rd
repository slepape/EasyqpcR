\name{totData}
\alias{totData}

\title{
Aggregation of qPCR biological replicates and data transformation
}
\description{
This function aggregates qPCR biological replicates and calculates the main
 parameters : mean (arithmetic or geometric), the standard deviation and the
 standard error from your biological replicates of your experience. This
 function has an algorithm published by Willems et al. (2008) which performs a
 standardization procedure that can be applied to data sets that display high
 variation between biological replicates. This enables proper statistical
 analysis and drawing relevant conclusions. The procedure is not new, it has
 been used in microarray data analysis and is based on a series of sequential
 corrections, including log transformation, mean centering, and autoscaling.
}
\usage{
totData(data, r, geo = TRUE, logarithm = TRUE, base, transformation = TRUE,
	nSpl, linear = TRUE, na.rm = na.rm)
}

\arguments{
  \item{data}{
data.frame containing row datas (genes in columns, samples in rows, Cq values).
}
  \item{r}{
numeric, number of qPCR replicates.
}
  \item{geo}{
logical, the function will use the geometrical mean of your biological
 replicates if TRUE or the arithmetic mean if FALSE.
}
  \item{logarithm}{
logical, the NRQs will be log-transformed.
}
  \item{base}{
numeric, the logarithmic base (2 or 10).
}
  \item{transformation}{
logical, if TRUE, the transformation procedure for highly variable biological
 replicates (but with the same tendency) will be done.
}
  \item{nSpl}{
numeric, the number of samples.
}
  \item{linear}{
logical, after the transformation procedure done, your raw data will be
 normalized (anti-log-transformed).
}
  \item{na.rm}{
logical, indicating whether NA values should be stripped before the computation
 proceeds.
}
}
\details{
The standardization procedure used in this function (if TRUE for the
 transformation argument) is based on the article of Willems et al. (2008).
 This function perform successively thEerror operations : log-transformation of
 your raw data, mean of your log-transformed data for each biological replicate,
 mean centering for each biological replicate, standard deviation of each
 mean-centered biological replicate, autoscaling of your data, i.e., your
 mean-centered data for each biological replicate will be divided by the
 standard deviation of the mean-centered biological replicate and then
 multiplicated by the mean of the standard deviation of all the biological
 replicates.

For more information for the way to use this function, please see the vignette.
}
\value{

 \item{Mean of your qPCR runs}{The geometric (if TRUE for geo) or arithmetic
 mean of your biological replicates.}
 \item{Standard deviations of your qPCR runs}{The standard deviation of your
 biological replicates.}
 \item{Standard errors of your qPCR runs}{The standard error of your biological
 replicates.}
 \item{Transformed data}{If TRUE for transformation, your raw data will be
 transformed by the algorithm of Willems et al. (2008).}
 \item{Reordered transformed data}{The transformed data reordered by rowname.}


}
\references{
Erik Willems Luc Leyns, Jo Vandesompele. Standardization of real-time PCR gene
 expression data from independent biological replicates. Analytical Biochemistry
 379 (2008) 127-129 (doi:10.1016/j.ab.2008.04.036).
 <url:http://www.sciencedirect.com/science/article/pii/S0003269708002649>
}
\author{
Sylvain Le pape <sylvain.le.pape@univ-poitiers.fr>
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


## Attenuation of inter-run variation thanks to the calibration factor for the
## first biological replicate

nrmData(data = qPCR_run1 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)

## Attenuation of inter-run variation thanks to the calibration factor for the
## second biological replicate

nrmData(data = qPCR_run2 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE) 

## Attenuation of inter-run variation thanks to the calibration factor for the
## third biological replicate

nrmData(data = qPCR_run3 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE) 

## Isolate the NRQs scaled to control of the first biological replicate

a1 <- nrmData(data = qPCR_run1 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[1] 

## Isolate the NRQs scaled to control of the second biological replicate 

b1 <- nrmData(data = qPCR_run2 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[1] 

## Isolate the NRQs scaled to control of the third biological replicate

c1 <- nrmData(data = qPCR_run3 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[1]

## Data frame transformation

a2 <- as.data.frame(a1)
b2 <- as.data.frame(b1)
c2 <- as.data.frame(c1)

## Aggregation of the three biological replicates

d2 <- rbind(a2, b2, c2)

totData(data=d2, r=3, geo=TRUE, logarithm=TRUE, base=2, 
	transformation=TRUE, nSpl=5, linear=TRUE,
	na.rm=TRUE)

}

\keyword{Biological replicates}
\keyword{Standardization procedure}
