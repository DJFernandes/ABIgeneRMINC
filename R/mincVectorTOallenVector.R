#' @title Rotates 1-D MINC data to Allen space
#' @description
#' MINC Data is stored as a 1-D array going from X=Left-to-Right, Y=Posterior-to-Anterior, Z=Inferior-to-Superior (dimensions written from fastest changing index to slowest)
#' Allen Gene Expression is stored as a 1-D array going from X=Anterior-to-Posterior, Y=Superior-to-Inferior, and Z=Left-to-Right (dimensions written from fastest changing index to slowest)
#' To compare Allen Gene expression data and MINC data, you have to rotate the Allen Gene Expression by 90 degrees about Z-axis and then 90 degrees about Y-axis
#' This function rotates a 1-D MINC vector to 1-D allen space oriendation
#'
#' @param minc.vector 1-D vector in MINC Space Orientation
#' @param xyzDimSize dimensions of MINC Space (in x,y,z). Default is c(58,67,41), the dimensions of Allen Gene Expression data
#'
#' @return 1-D vector in Allen space orientation
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
#'
#' # Rotate it to MINC Orientation
#' gene.expression.minc=allenVectorTOmincVector(gene.expression.allen)
#'
#' # Rotate it back to Allen Orientation
#' gene.expression=mincVectorTOallenVector(gene.expression.minc)
#' 
#' @export


mincVectorTOallenVector <- function(minc.vector,xyzDimSize=c(58,67,41)) {
	nx=xyzDimSize[1];ny=xyzDimSize[2];nz=xyzDimSize[3]
	x=0;z=nz
	idx=unlist(lapply(0:(nx-1),
		function(x) unlist(lapply((nz-1):0,
		function(z) z*ny*nx+((ny-1):0)*nx+x))))
	idx=idx+1
	return(minc.vector[idx])	
}
