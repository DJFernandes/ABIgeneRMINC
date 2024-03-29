% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allen_minc_vector_conversion.R
\name{allenVectorTOmincVector}
\alias{allenVectorTOmincVector}
\title{Rotates Allen data to MINC space}
\usage{
allenVectorTOmincVector(allen.data, xyzDimSize = NULL, orientation.check = T)
}
\arguments{
\item{allen.data}{data (either as array or 1D vector) in Allen Space Orientation}

\item{xyzDimSize}{dimensions of Allen Space (in x,y,z). Default is c(67,41,58), the dimensions of Allen Gene Expression data}

\item{orientation.check}{check orientation attribute.}
}
\value{
1-D vector in MINC space orientation
}
\description{
Allen Gene Expression is stored as a 1-D array going from X=Anterior-to-Posterior, Y=Superior-to-Inferior, and Z=Left-to-Right (dimensions written from fastest changing index to slowest)
MINC Data is stored as a 1-D array going from X=Left-to-Right, Y=Posterior-to-Anterior, Z=Inferior-to-Superior (dimensions written from fastest changing index to slowest)
To compare Allen Gene expression data and MINC data, you have to rotate the Allen Gene Expression by 90 degrees about Z-axis and then 90 degrees about Y-axis
This function rotates a 1-D allen vector to 1-D MINC space oriendation
}
\examples{
# Download gene expression files from the Allen Brain Institute
# In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
# I included it in this library as an example

# Read Gene expression file
filename=system.file('extdata/Pdyn_P56_coronal_71717084_200um.zip',package="ABIgeneRMINC")
gene.expression.allen=read.raw.gene(filename)

# Rotate it to MINC Orientation
gene.expression.minc=allenVectorTOmincVector(gene.expression.allen)

# Rotate it back to Allen Orientation
gene.expression=mincVectorTOallenVector(gene.expression.minc)

# Going to minc orientation and back should NOT change any voxel values. Check that...
all(gene.expression == gene.expression.allen)
}
