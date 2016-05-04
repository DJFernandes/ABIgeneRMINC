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
#' @export

read.raw.gene <- function(filename,labels=FALSE) {
	connec=file(filename, "rb")
	#dimensions of allen gene expression arrays are 67*41*58=159326
	if (labels) { typ="integer"} else {typ="numeric"}
	dat=readBin(connec,typ,size=4,n=159326)
	close(connec)
	return(dat)
}

