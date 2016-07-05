#' @title Rotates 1-D Allen data to MINC space
#' @description
#' Allen Gene Expression is stored as a 1-D array going from X=Anterior-to-Posterior, Y=Superior-to-Inferior, and Z=Left-to-Right (dimensions written from fastest changing index to slowest)
#' MINC Data is stored as a 1-D array going from X=Left-to-Right, Y=Posterior-to-Anterior, Z=Inferior-to-Superior (dimensions written from fastest changing index to slowest)
#' To compare Allen Gene expression data and MINC data, you have to rotate the Allen Gene Expression by 90 degrees about Z-axis and then 90 degrees about Y-axis
#' This function rotates a 1-D allen vector to 1-D MINC space oriendation
#'
#' @param allen.vector 1-D vector in Allen Space Orientation
#' @param xyzDimSize dimensions of Allen Space (in x,y,z). Default is c(67,41,58), the dimensions of Allen Gene Expression data

#' @return 1-D vector in MINC space orientation
#'
#' @examples
#' # Download gene expression files from the Allen Brain Institute
#' # In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
#' # Then Unzipped it, and renamed the energy.raw file to Pdyn_energy_experiment71717084.raw 
#' # I included it in this library as an example
#'
#' # Read Gene expression file
#' filename=system.file('extdata/Pdyn_energy_experiment71717084.raw',package="ABIgeneRMINC")
#' gene.expression.allen=read.raw.gene(filename)

#' # Rotate it to MINC Orientation
#' gene.expression.minc=allenVectorTOmincVector(gene.expression.allen)

#' @export


allenVectorTOmincVector <- function(allen.vector,xyzDimSize=c(67,41,58)) {
	if (length(xyzDimSize)!=3) {stop("Please specify the three dimension lengths as a vector")}
	nx=xyzDimSize[1];ny=xyzDimSize[2];nz=xyzDimSize[3]
	if (nx*ny*nz!=length(allen.vector)) {
	stop("vector length and product of dimension lengths are not the same")}
	idx=unlist(lapply((ny-1):0,
		function(y) unlist(lapply((nx-1):0,
		function(x) (0:(nz-1))*(nx*ny)+y*nx+x))))
	idx=idx+1
        ret=allen.vector[idx]
        #set sizes attributes
        attr(ret,'sizes')=c(41,67,58)
	return(ret)
}

