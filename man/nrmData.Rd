\name{nrmData}
\alias{nrmData}

\title{
Determination of the NF, RQ, NRQ, NRQ scaled to control and their SE and SD.
}
\description{
This function determines the values of the normalization factors, the relative
 quantitues, the normalized relative quantities, the normalized relative
 quantities scaled to control and their respectives standard errors and
 standard deviations by the method described by Hellemans et al. (2007).}
\usage{
nrmData(data, r, E, Eerror, nSpl, nbRef, Refposcol, nCTL, CF, CalPos, 
	trace = FALSE, geo = FALSE, na.rm = na.rm)
}

\arguments{
  \item{data}{
data.frame containing row datas (genes in columns, samples in rows, Cq values).
}
  \item{r}{
numeric, number of qPCR replicates.
}
  \item{E}{
numeric, amplification efficiency values for each gene (follow the same order of
 the genes).
}
  \item{Eerror}{
numeric, standard errors of amplification efficiencies for each gene (follow the
 same order of the genes).
}
  \item{nSpl}{
numeric, number of samples to analyzed.
}
  \item{nbRef}{
numeric, number of reference genes used.
}
  \item{Refposcol}{
column position of your reference gene(s).
}
  \item{nCTL}{
numeric, number of samples forming your control group.
}
  \item{CF}{
numeric (or object if you have used the calData function from this package),
 values of the calibration factors for each gene (follow the same order of the
 genes).
}
  \item{CalPos}{
numeric, sample number of your calibrator(s).
}
  \item{trace}{
logical, print additional information.
}
  \item{geo}{	
logical, to scale to your control group, the function will use the geometrical
 mean if TRUE or the arithmetic mean if FALSE.
}
  \item{na.rm}{
a logical value indicating whether NA values should be stripped before the
 computation proceeds.
}
}
\details{

The algorithm used in this function is based on the article of Hellemans et al.
 (2007). This function calculates the expression value scaled to your control
 group and normalized to the calibration factor and the normalization factor.
 The limiting step is that you need to put the control samples on the top of the
 data frame otherwise, the algorithm will not work correctly. For more
 information for the way to use this function, please see the vignette.

}
\value{

  \item{NRQs normalized to control}{Gives the normlized relative quantities
 scaled to your control group.}
  \item{NRQs}{Gives the normlized relative quantities.}
  \item{NRQs of your calibrator for this run}{Gives the normlized relative
 quantities of your calibrator(s).}
}
\references{
Jan Hellemans, Geert Mortier, Anne De Paepe, Frank Speleman and Jo Vandesompele.
 qBase relative quantification framework and software for management and
 automated analysis of real-time quantitative PCR data. Genome Biology 2007,
 8:R19 (doi:10.1186/gb-2007-8-2-r19). <url:http://genomebiology.com/2007/8/2/R19>
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
}

\keyword{qPCR expression}
