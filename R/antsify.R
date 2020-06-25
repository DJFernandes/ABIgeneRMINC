#' Convert RMINC and ABIgeneRMINC data to ANTsR images
#' @param x RMINC or ABIgeneRMINC data
#' @return ANTsR images
#' @examples
#' # Here is an example converting ABIgeneRMINC data to ANTsR image
#' nrp1Filename=system.file('extdata/Nrp1_P56_coronal_74272479_200um.zip',package="ABIgeneRMINC")
#' ABIvec = read.raw.gene(nrp1Filename)
#' antsvec1 = ANTsify(ABIvec)
#'
#' # Here is an example converting MINC data to ANTsR image
#' statFilename=system.file('extdata/enrichment_stats.mnc',package="ABIgeneRMINC")
#' mincvec = mincGetVolume(statFilename)
#' antsvec2 = ANTsify(mincvec)
#'
#' @export

ANTsify = function(x) {
   # add attributes if needed/possible
   if ( !'ANTsRparam' %in% names(attributes(x)) ) {
      if (is.minc(x)) {x = mincAddANTsAttrs(x)
      }  else {stop('No ANTsR attributes')}
   }

   # check if all attributes are there
   ANTsRparamNames = c('pixeltype','components','spacing','origin','direction')
   bool = ANTsRparamNames %in% names(attributes(x)$ANTsRparam)
   if (!all(bool)) {
      stop( ANTsRparamNames[!bool][1], 'not found in ANTsR attributes' )
   }

   # make 3D
   if (is.null(dim(x)))  {
      vol = array(x, dim = attr(x,'ANTsRparam')$imagesize)
   } else {
      vol = x
   }

  # strip classes if present
  if ('class' %in% names(attributes(x))) {
     attributes(vol)$class = NULL
  }
   
   as.antsImage(
              object     = vol                                                  ,
              pixeltype  = attr(x,'ANTsRparam')$pixeltype                       , 
              components = attr(x,'ANTsRparam')$components                      ,
              spacing    = attr(x,'ANTsRparam')$spacing                         ,
              origin     = attr(x,'ANTsRparam')$origin                          ,
              direction  = attr(x,'ANTsRparam')$direction                       )
}

#' Turn ANTs image into ABIgeneRMINC data
#' @param x ANTsR images
#' @return ABIgeneRMINC data
#' @examples
#' # Here is an example converting ABIgeneRMINC data to ANTsR image
#' nrp1Filename=system.file('extdata/Nrp1_P56_coronal_74272479_200um.zip',package="ABIgeneRMINC")
#' ABIvec = read.raw.gene(nrp1Filename)
#' antsvec = ANTsify(ABIvec)
#' y = deANTsify(antsvec)
#'
#' @export
deANTsify = function(x) {
   ret = as.array(x)
   attr(ret , 'sizes') = rev(dim(ret))
   attr(ret ,'ANTsRparam') = list(
                     imagesize  = dim(ret)                 ,
                     pixeltype  = attributes(x)$pixeltype  ,
                     components = attributes(x)$components ,
                     spacing    = spacing(x)               ,
                     origin     = origin(x)                ,
                     direction  = direction(x)             
                 )
   ret 
}

#' Add ANTsR attributes to MINC data
#' @param x RMINC data
#' @return RMINC data with ANTsR attributes
#' @examples
#'
#' # Here is an example converting MINC data to ANTsR image
#' statFilename=system.file('extdata/enrichment_stats.mnc',package="ABIgeneRMINC")
#' mincvec = mincGetVolume(statFilename)
#' antsvec = ANTsify(mincvec)
#'
#' @export

mincAddANTsAttrs = function(x) {
   # x = mincGetVolume('../mncfiles/mr_allen_space.mnc')
   ANTsRparam = list(
           imagesize = rev(minc.dimensions.sizes(likeVolume(x))),
           pixeltype = 'float',    # is it ever NOT a float?
           components = FALSE,     # I doubt anyone will ever NOT use a scalar
           spacing = rev( minc.separation.sizes(likeVolume(x)) * 1000 ),
           origin  = c(0,0,0),     # Usually when processing in R, user has already registered data
           direction = diag(3)    # we usually do things in orthogonal coordinate systems 
        )
   attr(x,'ANTsRparam')  = ANTsRparam
   if (is.null(attr(x,'orientation'))) {
      warning('Assuming data is in MINC orientation')
      attr(x,'orientation') = 'MINC'
   }
   x
}


