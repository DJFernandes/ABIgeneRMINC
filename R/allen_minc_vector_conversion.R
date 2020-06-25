#' @title Rotates Allen data to MINC space
#' @description
#' Allen Gene Expression is stored as a 1-D array going from X=Anterior-to-Posterior, Y=Superior-to-Inferior, and Z=Left-to-Right (dimensions written from fastest changing index to slowest)
#' MINC Data is stored as a 1-D array going from X=Left-to-Right, Y=Posterior-to-Anterior, Z=Inferior-to-Superior (dimensions written from fastest changing index to slowest)
#' To compare Allen Gene expression data and MINC data, you have to rotate the Allen Gene Expression by 90 degrees about Z-axis and then 90 degrees about Y-axis
#' This function rotates a 1-D allen vector to 1-D MINC space oriendation
#'
#' @param allen.data data (either as array or 1D vector) in Allen Space Orientation
#' @param xyzDimSize dimensions of Allen Space (in x,y,z). Default is c(67,41,58), the dimensions of Allen Gene Expression data

#' @return 1-D vector in MINC space orientation
#'
#' @examples
#' # Download gene expression files from the Allen Brain Institute
#' # In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
#' # I included it in this library as an example
#'
#' # Read Gene expression file
#' filename=system.file('extdata/Pdyn_P56_coronal_71717084_200um.zip',package="ABIgeneRMINC")
#' gene.expression.allen=read.raw.gene(filename)
#'
#' # Rotate it to MINC Orientation
#' gene.expression.minc=allenVectorTOmincVector(gene.expression.allen)
#' 
#' # Rotate it back to Allen Orientation
#' gene.expression=mincVectorTOallenVector(gene.expression.minc)
#' 
#' # Going to minc orientation and back should NOT change any voxel values. Check that...
#' all(gene.expression == gene.expression.allen)
#  # TRUE
#' @export
allenVectorTOmincVector <- function(allen.data,xyzDimSize=NULL) {
        save_attributes = attributes(allen.data)
        save_attributes = save_attributes[!names(save_attributes) %in% 'dim']

        if (is.null(dim(allen.data))) { d1 = T } else { 
            d1 = F 
            allen.data = as.vector(allen.data)
            attributes(allen.data) = save_attributes
        }

        if ( attr(allen.data,'orientation') != 'ABI' ) {
           warning(
                    'This file does not appear to be in Allen orientation.',
                    'I will keep going but your results may be wrong.'
                  )
        }
        if ( is.null(attr(allen.data,'sizes')) & is.null(xyzDimSize) ) {
           stop("Either supply the sizes attribute in allen.data or give xyzDimSize argument")
        }
        if ( !is.null(attr(allen.data,'sizes')) & !is.null(xyzDimSize) ) {
           if (all(rev(attr(allen.data,'sizes')) == xyzDimSize)) { 
              cat(
               'You supplied both sizes attribute in allen.data and xyzDimSize argument.',
               'This is reduntant but who am I to stop you from being redundant.',
                '\n')
           } else {
              warning(
               'You supplied both sizes attribute in allen.data and xyzDimSize argument',
               'and they do not agree. I will be only using the xyzDimSize argument')
           }
        }
        if (is.null(xyzDimSize)) { xyzDimSize = rev(attr(allen.data,'sizes')) }

	nx=xyzDimSize[1];ny=xyzDimSize[2];nz=xyzDimSize[3]
	if (nx*ny*nz!=length(allen.data)) {
	stop("vector length and product of dimension lengths are not the same")}
	idx=unlist(lapply((ny-1):0,
		function(y) unlist(lapply((nx-1):0,
		function(x) (0:(nz-1))*(nx*ny)+y*nx+x))))
	idx=idx+1
        ret=allen.data[idx]

        if (!d1) { ret = array(ret, dim = xyzDimSize[c(3,1,2)]) }
        attributes(ret) = c(attributes(ret),save_attributes)
        attr(ret,'sizes') = xyzDimSize[c(2,1,3)]          
        attr(ret,'ANTsRparam')$imagesize = xyzDimSize[c(3,1,2)]
        attr(ret,'ANTsRparam')$spacing   = save_attributes$ANTsRparam$spacing[c(3,1,2)]
        attr(ret,'ANTsRparam')$origin    = save_attributes$ANTsRparam$origin[c(3,1,2)]
        attr(ret,'ANTsRparam')$direction = save_attributes$ANTsRparam$direction[c(3,1,2),c(3,1,2)]

        # orientation (minc orientation is 
        #          X=Left-to-Right
        #          Y=Posterior-to-Anterior
        #          Z=Inferior-to-Superior)
        attr(ret,'orientation') = 'MINC'       
	return(ret)
}

#' @title Rotates MINC data to Allen space
#' @description
#' MINC Data is stored as a 1-D array going from X=Left-to-Right, Y=Posterior-to-Anterior, Z=Inferior-to-Superior (dimensions written from fastest changing index to slowest)
#' Allen Gene Expression is stored as a 1-D array going from X=Anterior-to-Posterior, Y=Superior-to-Inferior, and Z=Left-to-Right (dimensions written from fastest changing index to slowest)
#' To compare Allen Gene expression data and MINC data, you have to rotate the Allen Gene Expression by 90 degrees about Z-axis and then 90 degrees about Y-axis
#' This function rotates a 1-D MINC vector to 1-D allen space oriendation
#'
#' @param minc.data data (either as array or 1D vector) in MINC Space Orientation
#' @param xyzDimSize dimensions of MINC Space (in x,y,z). Default is c(58,67,41), the dimensions of Allen Gene Expression data
#'
#' @return 1-D vector in Allen space orientation
#'
#' @examples
#' # Download gene expression files from the Allen Brain Institute
#' # In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
#' # I included it in this library as an example
#'
#' # Read Gene expression file
#' filename=system.file('extdata/Pdyn_P56_coronal_71717084_200um.zip',package="ABIgeneRMINC")
#' gene.expression.allen=read.raw.gene(filename)
#'
#' # Rotate it to MINC Orientation
#' gene.expression.minc=allenVectorTOmincVector(gene.expression.allen)
#' 
#' # Rotate it back to Allen Orientation
#' gene.expression=mincVectorTOallenVector(gene.expression.minc)
#' 
#' # Going to minc orientation and back should NOT change any voxel values. Check that...
#' all(gene.expression == gene.expression.allen)
#' # TRUE
#' @export
mincVectorTOallenVector <- function(minc.data,xyzDimSize=NULL) {
        save_attributes = attributes(minc.data)
        save_attributes = save_attributes[!names(save_attributes) %in% 'dim']

        if (is.null(dim(minc.data))) { d1 = T } else { 
            d1 = F 
            minc.data = as.vector(minc.data)
            attributes(minc.data) = save_attributes
        }

        if ( !any(grepl('^minc',attr(minc.data, 'class'))) ) {
           if ( attr(minc.data,'orientation') != 'MINC' ) {
              warning(
                       'This file does not appear to be in MINC orientation.',
                       'I will keep going but your results may be wrong.'
                     )
           }
        }
        if ( !is.null(attr(minc.data,'sizes')) & !is.null(xyzDimSize) ) {
           if (all(rev(attr(minc.data,'sizes')) == xyzDimSize)) { 
              cat(
               'You supplied both sizes attribute in minc.data and xyzDimSize argument.',
               'This is reduntant but who am I to stop you from being redundant.',
                '\n')
           } else {
              warning(
               'You supplied both sizes attribute in minc.data and xyzDimSize argument',
               'and they do not agree. I will be only using the xyzDimSize argument')
           }
        }
        if (is.null(xyzDimSize)) { xyzDimSize = rev(attr(minc.data,'sizes')) }

	nx=xyzDimSize[1];ny=xyzDimSize[2];nz=xyzDimSize[3]
	x=0;z=nz
	idx=unlist(lapply(0:(nx-1),
		function(x) unlist(lapply((nz-1):0,
		function(z) z*ny*nx+((ny-1):0)*nx+x))))
	idx=idx+1
        ret=minc.data[idx]

        if (!d1) { ret = array(ret, dim = xyzDimSize[c(2,3,1)]) }
        attributes(ret) = c(attributes(ret),save_attributes)
        attr(ret,'sizes') = xyzDimSize[c(1,3,2)]          
        attr(ret,'ANTsRparam')$imagesize = xyzDimSize[c(2,3,1)]
        attr(ret,'ANTsRparam')$spacing   = save_attributes$ANTsRparam$spacing[c(2,3,1)]
        attr(ret,'ANTsRparam')$origin    = save_attributes$ANTsRparam$origin[c(2,3,1)]
        attr(ret,'ANTsRparam')$direction = save_attributes$ANTsRparam$direction[c(2,3,1),c(2,3,1)]

        # orientation (Allen orientation is 
        #          X=Anterior-to-Posterior
        #          Y=Superior-to-Inferior
        #          Z=Left-to-Right)
        attr(ret,'orientation') = 'ABI'
	return(ret)
}

