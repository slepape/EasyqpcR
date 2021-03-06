\name{badCt}
\alias{badCt}

\title{
Evaluation of the qPCR technical replicates
}
\description{
This function allows you to evaluate your qPCR technical replicates, you only
 need to define the threshold (according to Hellemans et al. (2007), 0.5 is a
 good threshold value), the dataset and the number of technical replicates you
 have done. I recommand you to use the gWidgets package to easily exclude the
 failed replicates.
}
\usage{
badCt(data, r, threshold, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data.frame containing row datas (genes in columns, samples in rows, Cq values).
}
  \item{r}{
numeric, number of qPCR replicates.
}
  \item{threshold}{
numeric, the maximal variation between your qPCR replicates.
}
  \item{na.rm}{
a logical value indicating whether NA values should be stripped before the
 computation proceeds.
}
}
\details{
To facilitate the use of the function, I suggest you to use the gWidgets package
 as described in the vignette.
}
\value{
This function returns the position (sample position and column position) where
 the variation between qPCR replicates is superior to the threshold value.
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
data(qPCR_run1)

badCt(data=qPCR_run1, r=3, threshold=0.3, na.rm=TRUE)
}


\keyword{qPCR replicates}
