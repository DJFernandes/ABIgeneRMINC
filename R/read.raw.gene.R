#' Read Raw Allen Gene Expression File
#'
#' @param filename Filename of Raw Gene Expression File
#' @param labels logical. If TRUE, then file is read as a string of integers
#'
#' @return 1D vector of length 159326 containing Raw Gene Expression data
#'
#' @examples
#' # Download gene expression files from the Allen Brain Institute
#' # In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
#' # Then Unzipped it, and renamed the energy.raw file to Pdyn_energy_experiment71717084.raw 
#' # I included it in this library as an example
#'
#' filename=system.file('extdata/Pdyn_energy_experiment71717084.raw',package="ABIgeneRMINC")
#' gene.expression=read.raw.gene(filename)
#' 
#' ##Make coronal slice image
#' # First convert data from 1-D vector to 3-D array
#' gene.expression.array=array(gene.expression,c(67,41,58))
#' 
#' #Coronal slice
#' image(gene.expression.array[10,,], ylab='Left-Right' ,xlab='Superior-Inferior')
#' 
#' #Sagittal Slice
#' image(gene.expression.array[,,5], ylab='Superior-Inferior' ,xlab='Anterior-Posterior')
#' @export

read.raw.gene <- function(filename,labels=FALSE) {
	connec=file(filename, "rb")
	#dimensions of allen gene expression arrays are 67*41*58=159326
	if (labels) { typ="integer"} else {typ="numeric"}
	dat=readBin(connec,typ,size=4,n=159326)
	close(connec)
	return(dat)
}

