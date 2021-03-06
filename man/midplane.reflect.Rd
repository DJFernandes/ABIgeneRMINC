% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/midplane.reflect.R
\name{midplane.reflect}
\alias{midplane.reflect}
\title{Reflect Gene Expression Vector across midplane}
\usage{
midplane.reflect(
  gene.expression,
  reflect.dim = NULL,
  minc.array.sizes = NULL,
  data = "gene"
)
}
\arguments{
\item{gene.expression}{gene expression data to reflect (as either vector or array)}

\item{data}{If "gene" (default), then assume missing data has value -1, which is the syntax chosen by ABI for their gene expression data. If "atlas", then missing data is 0, which is the syntax chosen by ABI`s for their atlases.}

\item{reflect.dimension}{INTEGER. Which dimension to reflect. This is done automatically if you use vectors and arrays generated by this package. If you want to do it manually, put 1 for the slowest dimension, 3 for the fastest. To reflect Allen orientation vectors across the sagittal midplane, set this value to 1. To reflect MINC orientation vectors across the sagittal midplane, set this value to 3.}

\item{array.dimensions}{Vector of length 3 indicating the dimension sizes in order of the SLOWEST dimension to the FASTEST. Must be supplied if gene.expression.vector is missing 'sizes' attributes (execute attr(gene.expression.vector,'sizes')).}
}
\value{
Gene Expression as vector or array
}
\description{
Reflect Gene Expression Vector across midplane
}
\examples{
# Read Bdnf gene expression (Experiment 75695642) from Allen Brain Institute
# Look at a coronal slice
# Reflect across the sagittal midplane 
# Look at a coronal slice


## Read Bdnf expression (Experiment 75695642)
URLtoDownload='http://api.brain-map.org/grid_data/download/75695642?include=energy'
gene.expression.vector=read.raw.gene(URLtoDownload,url=TRUE)

## Coronal slice
c.slice=array(gene.expression.vector,rev(attr(gene.expression.vector,'sizes')))[30,,]
image(c.slice, ylab='Left-Right' ,xlab='Superior-Inferior')

## Reflect across the sagittal midplane 
reflected.vector=midplane.reflect(gene.expression.vector,1)

## Coronal slice
reflected.c.slice=array(reflected.vector,rev(attr(gene.expression.vector,'sizes')))[30,,]
image(reflected.c.slice, ylab='Left-Right' ,xlab='Superior-Inferior')

## Reflect MINC orientation vector across sagittal midplane
minc.vector=allenVectorTOmincVector(gene.expression.vector)
reflected.minc.vector=midplane.reflect(minc.vector,3)
}
