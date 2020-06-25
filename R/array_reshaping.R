#' Pad array to fit a different shape
#' @param x Array to pad
#' @param new.dims new dimensions
#' @param side should padding be applied to the left, right, or both (ties are broken by the left)?
#' @param padval what value should be used for padding?
#' @return array padded so it has dimensions of new.dims
#' @examples
#' # read a gene expression file
#' filename=system.file('extdata/Pdyn_P56_coronal_71717084_200um.zip',package="ABIgeneRMINC")
#' gene.expression=read.raw.gene(filename)
#' 
#' # make it into an array 
#' gene.expression.array=array(gene.expression,rev(attr(gene.expression,'sizes')))
#' # synonymous with mincArray(gene.expression))
#' 
#' # dimensions of array
#' dim(gene.expression.array)
#' # 67 41 58
#' 
#' # pad the array
#' padded.gene.expression.array = gene_pad_arr(gene.expression.array,new.dims=c(100,99,98),side='both')
#' 
#' # dimensions of padded array
#' dim(padded.gene.expression.array)
#' # 100 99 98
#' @export

gene_pad_arr = function(x,new.dims,side=c('left','right','both'), padval=0) {
   save_attributes = attributes(x)
   if (length(side) != 1) {
      stop('side must be one of left, right, or both')
   } else if ( !side %in% c('left','right','both') ) {
      stop('side must be one of left, right, or both')
   }
   if (length(dim(x)) != length(new.dims)) {
      stop('The dimensions of the array must match the length of new.dims')
   }

   old.dims = dim(x)
   dim.change = new.dims - old.dims
   if (any(dim.change<0)) { 
      stop('Dimension',which(dim.change<0)[1],'is smaller in new.dims than x array')
   }
   for (i in 1:length(dim.change)) {
      if (dim.change[i] == 0) { next }
      pad.arr.dim = old.dims ; pad.arr.dim[i] = dim.change[i]
      pad.arr = array(padval,dim=pad.arr.dim)

      if (side == 'both') {
         boolvec = (floor(2*((1:dim.change[i])-1)/dim.change[i]) == 0)
         if (i==1) {
            pad.arr.left  = pad.arr[ boolvec,,]
            pad.arr.right = pad.arr[!boolvec,,]
         } else if (i==2) {
            pad.arr.left  = pad.arr[, boolvec,]
            pad.arr.right = pad.arr[,!boolvec,]
         } else if (i==3) {
            pad.arr.left  = pad.arr[,, boolvec]
            pad.arr.right = pad.arr[,,!boolvec]
         }
         x = abind::abind(pad.arr.left, x, pad.arr.right, along=i)
      } else if (side == 'right') {
         x = abind::abind(x, pad.arr, along=i)
      } else if (side == 'left') {
         x = abind::abind(pad.arr, x, along=i)
      }
      old.dims[i] = new.dims[i]
      dim.change[i] = 0
   }

   if ('sizes' %in% names(save_attributes)) {
      attr(x,'sizes') = rev(attr(x,'dim'))
   } 
   if ('ANTsRparam' %in% names(save_attributes)) {
      attr(x,'ANTsRparam') = save_attributes$ANTsRparam
      attributes(x)$ANTsRparam$imagesize = attr(x,'dim')
   } 
   if ('orientation' %in% names(save_attributes)) {
      attr(x,'orientation') = save_attributes$orientation
   } 
   return(x)
}


#' Resample gene expression to a different resolution
#' @param gene.expression Gene expression data read with read.raw.gene (can be 1-D vector or array). If 1-D vector, it must has the sizes attribute
#' @param spacing new resolution to sample either as a vector (one component for each dimension) or scalar. If scalar, new resolution is isotropic.
#' @param interpolation Default is 'cubic'. 'Linear' and 'nearest' are also currently supported.
#' @return resampled gene expression
#' @examples
#' # read a gene expression file
#' filename=system.file('extdata/Pdyn_P56_coronal_71717084_200um.zip',package="ABIgeneRMINC")
#' gene.expression=read.raw.gene(filename)
#' 
#' # convert data from 1-D vector to 3-D array
#' gene.expression.array=array(gene.expression,rev(attr(gene.expression,'sizes')))
#' # synonymous with mincArray(gene.expression))
#' 
#' #Sagittal Slice
#' image(gene.expression.array[,,29], ylab='Superior-Inferior' ,xlab='Anterior-Posterior')
#' title('Native resolution', font.main=2)
#' 
#' # The current resolution is isotropic 200 microns
#' attributes(gene.expression)$ANTsRparam$spacing
#' # 200 200 200
#' 
#' # resample to isotropic 50 microns
#' resampled.gene.expression = gene_resample_resolution(gene.expression, c(50,50,50))
#' 
#' # The new resolution is
#' attributes(resampled.gene.expression)$ANTsRparam$spacing
#' # 50 50 50
#' 
#' # make sagittal slice
#' resampled.gene.expression.array=array(resampled.gene.expression,rev(attr(resampled.gene.expression,'sizes')))
#' image(resampled.gene.expression.array[,,29*4], ylab='Superior-Inferior' ,xlab='Anterior-Posterior')
#' title('4x super-sampling', font.main=2)
#' @export

gene_resample_resolution = function(gene.expression, spacing = NULL, interpolation = 'cubic') {

   if (is.null(spacing)) { stop('must supply new resolution') }
   if ( is.null(dim(gene.expression)) & is.null(attr(gene.expression,'sizes')) ) {
         stop('gene.expression must either have the sizes attribute or be an array')
   }
   if (is.null(attr(gene.expression,'ANTsRparam')$spacing)) {
         stop('gene.expression must have spacing attribute')
   }
   if (is.null(dim(gene.expression))) {
     d1 = T
     gene.arr = array(gene.expression, dim = rev(attr(gene.expression,'sizes')))
   } else {
     d1 = F
     gene.arr = gene.expression
   }
   if (length(spacing) == 1) {
      spacing = rep(spacing, length(dim(gene.arr)))
   } else if ( length(spacing) != length(dim(gene.arr)) ) {
      stop('Spacing must either be scalar (for isotropic sampling), or its length equal to number of dimensions')
   }
   nlength = (dim(gene.arr)-1) * attr(gene.expression,'ANTsRparam')$spacing
   sampling_grid = lapply(1:length(dim(gene.arr)), function(i) {
      world_coord = seq(from=0, to=nlength[i], by=spacing[i])
      vxl_coord = ( world_coord/attr(gene.expression,'ANTsRparam')$spacing[i] ) + 1
   })

   if (interpolation == 'cubic') {
      retimg = mmand::resample(gene.arr, points=sampling_grid, kernel = mmand::mnKernel())
   } else if (interpolation == 'nearest') {
      retimg = mmand::resample(gene.arr, points=sampling_grid, kernel = mmand::boxKernel())
   } else if (interpolation == 'linear') {
      retimg = mmand::resample(gene.arr, points=sampling_grid, kernel = mmand::triangleKernel())
   } else {
      stop('Only interpolation parameters "cubic", "nearest", and "linear" are currently supported')
   }

   odims = dim(retimg)
   if (d1) {  retimg = as.vector(retimg) }
   
   save_attributes = attributes(gene.expression)
   save_attributes = save_attributes[names(save_attributes) != 'dim']

   attributes(retimg) = c(attributes(retimg), save_attributes)
   attr(retimg, 'sizes') = rev(odims)
   attr(retimg, 'ANTsRparam')$imagesize = odims
   attr(retimg, 'ANTsRparam')$spacing   = spacing

   return(retimg)

}

#' Rescale Allen data
#'
#' @param img Image you want to resample (either a 1D vector or array)
#' @param factor How much down/up-sampling. Less than 1 to downsample, greater than 1 to upsample.
#'
#' @return resampled image
#'
#' @examples
#' 
#' 
#' # download the Allen Brain Institute p56 nissl atlas
#' p56vec = read.raw.atlas('http://download.alleninstitute.org/informatics-archive/current-release/mouse_annotation/P56_atlasVolume.zip', url=T)
#' 
#' # convert to a volume
#' p56vol = array(p56vec, dim = attr(p56vec,'ANTsRparam')$imagesize)
#' 
#' # Native file is 25 microns isotropic. Downsample to 100 microns isotropic
#' p56vol.ds = rescale.abi(p56vol,0.25)
#' 
#' # plot
#' dev.new(width=9,height=9)
#' RMINC::mincImage(p56vol, dimension = 2, slice=200, low=0, high=0.5, underTransparent = FALSE, col = colorRampPalette(c('black','white'))(100) )
#' title('Native Resolution', font.main=2)
#' 
#' dev.new(width=9,height=9)
#' RMINC::mincImage(p56vol.ds, dimension = 2, slice=50, low=0, high=0.5, underTransparent = FALSE, col = colorRampPalette(c('black','white'))(100) )
#' title('4x Downsampled', font.main=2)
#' @export

rescale.abi <- function(img, factor, interpolation = 'cubic') {
   # is image a 1-D vector of higher dimensional data?
   if (is.null(dim(img))) {
      d1 = T
      if (length(img) != prod(attr(img,'ANTsRparam')$imagesize)) { 
         stop('Dimenstion mis-match found')
      }
   }  else {
      d1 = F
   }
   if (d1) {
      vol = array(img, dim = attr(p56vec,'ANTsRparam')$imagesize)
      attributes(vol) = attributes(img)
      img = vol
   }

   if (length(factor) != length(dim(img)) )  {
     if (length(factor) != 1) {
        stop('factor must either be of length 1 or image dimension (',length(dim(img)),')')
     } else {
        factor = rep(factor, length(dim(img)))
     } 
   }

   if (interpolation == 'cubic') {
      retimg = mmand::rescale(img, factor , mmand::mnKernel())
   } else if (interpolation == 'nearest') {
      retimg = mmand::rescale(img, factor , mmand::boxKernel())
   } else if (interpolation == 'linear') {
      retimg = mmand::rescale(img, factor , mmand::triangleKernel())
   } else {
      stop('Only interpolation parameters "cubic", "nearest", and "linear" are currently supported')
   }

   if (d1) { retimg = as.vector(img) }
   

   new_attributes = attributes(img)
   new_attributes = new_attributes[names(new_attributes) != 'dim']
   new_attributes$sizes = rev(dim(retimg))
   new_attributes$ANTsRparam$imagesize = dim(retimg)
   new_attributes$ANTsRparam$spacing = attributes(img)$ANTsRparam$spacing/factor

   
   attributes(retimg) = c(attributes(retimg), new_attributes)

   return(retimg)
}


