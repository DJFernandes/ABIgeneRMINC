#' @title Interpolate Gene Expression
#' @description
#' Missing gene expression data is interpolated with a marching average of nearest neighbours. Neighbours are averaged if they have gene expression. After several iterations it will cover the whole brain.
#' 
#' @param genevec 1D vector of gene expression data. Missing data has value of -1.
#' @param maskvec 1D logical vector indicating which elements of genevec to consider for analysis. When function finishes running and converges, all genevec elements within the mask will have positive values.
#' @param dimensions 1D integer vector of length 3 specifying array dimensions (fastest to slowest dimensions)
#' @param itermax Maximum number of iterations to run the marching average of nearest neighbours.
#'
#' @return 1D vector of gene expression data, with missing data interpolated from nearest neighbours.
#'
#' @examples
#' # In this example, we will download gene expression of Nrxn1 (Experiment 75988632).
#' # Missing data will be interpolated using marching average of nearest neighbours
#' # 
#' # Download and read Nrxn Data
#' genetest=read.raw.gene("http://api.brain-map.org/grid_data/download/75988632?include=energy",url=TRUE)
#' 
#' # Read label file and create a mask of elements in the brain
#' labelfile=system.file('extdata/gridAnnotation.raw',package="ABIgeneRMINC")
#' mask=read.raw.gene(labelfile,labels=TRUE)>0
#' 
#' # Make image of horizontal slice before interpolation
#' image(array(genetest,rev(attributes(genetest)$sizes))[,25,],ylab='Left-Right' ,xlab='Anterior-Posterior')
#' # Make image of sagittal slice before interpolation
#' image(array(genetest,rev(attributes(genetest)$sizes))[,,25], ylab='Superior-Inferior' ,xlab='Anterior-Posterior')
#' 
#' # Interpolate
#' interp.gene=interpolate.gene(genetest,mask)
#' 
#' # Make image of horizontal slice after interpolation
#' image(array(interp.gene,rev(attributes(interp.gene)$sizes))[,25,],ylab='Left-Right' ,xlab='Anterior-Posterior')
#' # Make image of sagittal slice before interpolation
#' image(array(interp.gene,rev(attributes(interp.gene)$sizes))[,,25], ylab='Superior-Inferior' ,xlab='Anterior-Posterior')
#' 
#'
#' @useDynLib ABIgeneRMINC
#' @export


interpolate.gene=function(genevec,maskvec,dimensions=rev(attr(genevec,'sizes')),itermax=1000) {

  if (is.null(dimensions)) {stop('MUST specifiy dimensions of genevec (integer vector of length 3), or have sizes attribute')}

  # Load Fortran Code
#  dyn.load("interpgene.so")
  l=0.9
  # Find number of missing voxels. Interpolate till this number reaches 0
  missing=sum(genevec[mask]==-1)
  # If iteration does not reduce missing voxels, broaden interpolation criteria.
  prevmissing=missing
  # Counter till you reach max iteration, then break
  counter=0
    while (missing!=0) {
    counter=counter+1 ; if (counter > itermax) {print("Failed Convergence")}
    ret=.Fortran("interpgene",
         gene=as.numeric(genevec),
         mask=as.logical(maskvec),
         nx=as.integer(dimensions[1]),
         ny=as.integer(dimensions[2]),
         nz=as.integer(dimensions[3]),
         l=as.numeric(l),
         genei=as.numeric(genevec))$genei
    genevec=ret
    missing=sum(genevec[mask]==-1)
    if (missing == prevmissing) {l=l*0.75}
    prevmissing=missing
    }
  attr(ret,'sizes')=rev(dimensions)
  return(ret)
}


