% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/antsify.R
\name{mincAddANTsAttrs}
\alias{mincAddANTsAttrs}
\title{Add ANTsR attributes to MINC data}
\usage{
mincAddANTsAttrs(x)
}
\arguments{
\item{x}{RMINC data}
}
\value{
RMINC data with ANTsR attributes
}
\description{
Add ANTsR attributes to MINC data
}
\examples{

# Here is an example converting MINC data to ANTsR image
statFilename=system.file('extdata/enrichment_stats.mnc',package="ABIgeneRMINC")
mincvec = mincGetVolume(statFilename)
antsvec = ANTsify(mincvec)

}
