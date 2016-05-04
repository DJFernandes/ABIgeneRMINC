#' Write 1D vector data to ABI raw file
#'
#' @param gene.vector 1D vector of length 159326 containing Raw Gene Expression data
#' @param filename Name of output file
#'
#'
#' @examples
#' # In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
#' # Then Unzipped it, and renamed the energy.raw file to Pdyn_energy_experiment71717084.raw 
#' # I included it in this library as an example
#' 
#' # Read Gene Expression Data
#' filename=system.file('extdata/Pdyn_energy_experiment71717084.raw',package="ABIgeneRMINC")
#' gene.expression=read.raw.gene(filename)
#'
#' # Write Gene Expression Data
#' write.raw.gene(gene.expression,"outputfile.raw")
#'
#' # 'outputfile.raw' can be read again using "read.raw.gene" function
#' @export


write.raw.gene <- function(gene.vector,filename) {
	connec=file(filename, "wb")
	writeBin(gene.vector,con=connec,size=4)
	close(connec)
}
