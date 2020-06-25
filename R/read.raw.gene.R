#' Read Raw Allen Gene Expression File
#'
#' @param filename Filename (or URL) of Raw Gene Expression File
#' @param labels logical. If TRUE, then file is read as a string of integers
#' @param url logical. If TRUE, then filename MUST be a valid URL pointing to the allen gene expression data. The zip file will be downloaded into tmp and the gene expression data extracted and read. It is recommended you download and unzip the file yourself rather than relying on this flag. Debugging becomes much harder when you use this flag.
#'
#' @return 1D vector containing Raw Gene Expression data
#'
#' @examples
#' # Download gene expression files from the Allen Brain Institute
#' # In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
#'
#' filename=system.file('extdata/Pdyn_P56_coronal_71717084_200um.zip',package="ABIgeneRMINC")
#' gene.expression=read.raw.gene(filename)
#' 
#' ##Make coronal slice image
#' # First convert data from 1-D vector to 3-D array
#' gene.expression.array=array(gene.expression,rev(attr(gene.expression,'sizes')))
#' 
#' #Coronal slice
#' image(gene.expression.array[10,,], ylab='Left-Right' ,xlab='Superior-Inferior')
#' 
#' #Sagittal Slice
#' image(gene.expression.array[,,5], ylab='Superior-Inferior' ,xlab='Anterior-Posterior')
#' @export

read.raw.gene <- function(filename,labels=FALSE,url=FALSE) {
        # check if filename exists (or you have internet connection)
        rrg.file.check(filename, url, raw.allowed=T)
        
        # open connection
       	if (url==TRUE) {
           filename = rrg.download.url(filename)
        }
        if (is.zip(filename)) {
           flist = unzip(filename,list=T)$Name
           connec = unz(filename, grep('.raw$',flist,value=T),'rb')

           # get info from header
           headerconnec = unz(filename, grep('.mhd$',flist,value=T))
           headerinfo = ABI.info.from.header(headerconnec)
           close(headerconnec)

           dimsize   = headerinfo[[1]]
           size      = headerinfo[[2]]
           datatype  = headerinfo[[3]]
           pixeltype = headerinfo[[4]]
           antsparam = headerinfo[[5]]

           what = ifelse(pixeltype %in% c('unsigned char','unsigned int'), 'integer', 'numeric')
           dat=readBin(connec,what,size=size,n=prod(dimsize))

           skip.ants=F
        } else {
           connec = file(filename, "rb")
           if (labels) { what="integer"} else {what="numeric"}

           # based on the File size, make assumptions about what age the data is from
           dimsize = rrg.get.dimsize.from.file(filename)
           dat=readBin(connec,what,size=4,n=prod(dimsize))
           
           skip.ants=T
        }
	
        #close connections
	close(connec)
        if (url==TRUE) { unlink(filename) }

        #set sizes attributes (reversed so it plays nice with RMINC)
        attr(dat,'sizes')=rev(dimsize)

        if (!skip.ants) {
           #parameters for ANTsR
           attr(dat,'ANTsRparam')=antsparam
        }

        # orientation (Allen orientation is 
        #          X=Anterior-to-Posterior
        #          Y=Superior-to-Inferior
        #          Z=Left-to-Right)
        attr(dat,'orientation') = 'ABI'

	return(dat)
}

#' Read Raw Allen Atlas File
#'
#' @param filename Filename (zip or URL) of Raw Atlas File
#' @param url logical. If TRUE, then filename MUST be a valid URL pointing to the allen gene expression data. The zip file will be downloaded into tmp and the gene expression data extracted and read. It is recommended you download and unzip the file yourself rather than relying on this flag. Debugging becomes much harder when you use this flag.
#'
#' @return 3D array containing atlas
#'
#' @examples
#' atlas_url = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_annotation/P56_atlasVolume.zip'
#' 
#' # Read the atlas
#' nissl_atlas=read.raw.atlas(atlas_url, url=T)
#' 
#' ##Make coronal slice image
#' nissl_atlas_volume=array(nissl_atlas,rev(attr(nissl_atlas,'sizes')))
#' 
#' #Coronal slice
#' colscl = colorRampPalette(c('black','white'))(100)
#' image(nissl_atlas_volume[300,,], ylab='Left-Right' ,xlab='Superior-Inferior', zlim=c(0,0.9), col = colscl)
#' 
#' #Sagittal Slice
#' image(nissl_atlas_volume[,,200], ylab='Superior-Inferior' ,xlab='Anterior-Posterior', zlim=c(0,0.5), col = colscl)
#' @export

read.raw.atlas <- function(filename,url=FALSE) {
        # check if filename exists (or you have internet connection)
        rrg.file.check(filename, url,raw.allowed=F)
        
        # open connection
       	if (url==TRUE) {
           filename = rrg.download.url(filename)
        }
        if (is.zip(filename)) {
           flist = unzip(filename,list=T)$Name
           connec = unz(filename, grep('.raw$',flist,value=T),'rb')

           # get info from header
           headerconnec = unz(filename, grep('.mhd$',flist,value=T))
           headerinfo = ABI.info.from.header(headerconnec)
           close(headerconnec)

           dimsize   = headerinfo[[1]]
           size      = headerinfo[[2]]
           datatype  = headerinfo[[3]]
           pixeltype = headerinfo[[4]]
           antsparam = headerinfo[[5]]
        }

	dat=readBin(connec,'integer',size=size,n=prod(dimsize))
        #close connections
	close(connec)
        if (url==TRUE) { unlink(filename) }

        if (datatype == 'atlas') {
           # the scaling for the nissl atlases are kinda annoying.
           # they are integers with some negative values
           # so I am going to scale them
           dat = abs(dat)
           dat = dat / max(abs(dat))
           antsparam$pixeltype='float'
        }

        #set sizes attributes (reversed so it plays nice with RMINC)
        attr(dat,'sizes')=rev(dimsize)

        #parameters for ANTsR
        attr(dat,'ANTsRparam')=antsparam

        # orientation (Allen orientation is 
        #          X=Anterior-to-Posterior
        #          Y=Superior-to-Inferior
        #          Z=Left-to-Right)
        attr(dat,'orientation') = 'ABI'

	return(dat)
}

# Stolen from https://stackoverflow.com/questions/55283744/how-can-i-check-if-a-file-is-zip-or-not-in-r
is.zip <- function(filepath){
  result <- tryCatch({
              unzip(filepath, list = TRUE)
              return(TRUE)
            }, error = function(e){
              return(FALSE)
            })
  return(result)
}

rrg.file.check = function(filename, url, raw.allowed=T) {
        # check if filename exists
       	if (url==TRUE) { 
          if (!RCurl::url.exists(filename)) {stop(paste(filename,"does not exist"))}
        } else {
          if (!file.exists(filename)) { stop(paste(filename,"does not exist"))}
          if (!raw.allowed) {
             if (!is.zip(filename)) {stop(paste(filename,'must be a zipfile or url'))}
          }
        }
}

rrg.download.url = function(url_link) {
           temp = tempfile()
           download.file(url_link,temp)
           return(temp)
}

ABI.info.from.header = function(header) {
           header = readLines(header)
           dimsize = as.integer(
                        strsplit(
                           gsub(
                              '^DimSize = ','',
                              grep('^DimSize = ', header,value=T)
                              ),split = ' ')[[1]])
           spacing   = as.numeric(
                        strsplit(
                           gsub(
                              '^ElementSpacing = ','',
                              grep('^ElementSpacing = ', header,value=T)
                              ),split = ' ')[[1]])
           origin    = as.numeric(
                        strsplit(
                           gsub(
                              '^CenterOfRotation = ','',
                              grep('^CenterOfRotation = ', header,value=T)
                              ),split = ' ')[[1]])
           # I don't know if allen like row-major or column major orientation
           #  hopefully they only have symmetric matrices as directions...
           direction = as.numeric(
                        strsplit(
                           gsub(
                              '^TransformMatrix = ','',
                              grep('^TransformMatrix = ', header,value=T)
                              ),split = ' ')[[1]])
           n = sqrt(length(direction))
           if (round(n) != n) {stop('Direction Matrix is not square')}
           direction = matrix(direction, nrow=n)

           elementtype = gsub( '^ElementType = ', '', 
                  grep('^ElementType = ', header,value=T))
           if        ( elementtype == 'MET_UCHAR' ) { 
                size = 1
                datatype  = 'atlas'
                pixeltype = 'unsigned char'
           } else if ( elementtype == 'MET_UINT'  ) {
                size = 4
                datatype  = 'annotations'
                pixeltype = 'unsigned int'
           } else  if ( elementtype == 'MET_FLOAT'  ) {
                size = 4
                datatype  = 'gene'
                pixeltype = 'float'
           } else { stop('Elements must be either integers or floats') }

           # parameters for ants
           antsparam = list(
                             imagesize  = dimsize   ,
                             pixeltype  = pixeltype , 
                             components = F         ,
                             spacing    = spacing   , 
                             origin     = origin    ,
                             direction  = direction
                        )

           return(list(
             dimsize = dimsize, 
             size = size, 
             datatype = datatype,
             pixeltype = pixeltype,
             antsparam = antsparam
            ))
}



rrg.get.dimsize.from.file = function(filename) {
           # based on the File size, make assumptions about what age the data is from
           flsize = file.info(filename)$size/4
           if        (flsize == 159326) { 
              dimsize = c(67,41,58)
              cat('Assuming gene expression data is from adult brains\n')
           } else if (flsize == 158629) { 
              dimsize = c(73,41,53)
              cat('Assuming gene expression data is from P28 brains\n')
           } else if (flsize == 136000) { 
              dimsize = c(68,40,50)
              cat('Assuming gene expression data is from P14 brains\n')
           } else if (flsize == 165550) { 
              dimsize = c(77,43,50)
              cat('Assuming gene expression data is from P4 brains\n')
           } else if (flsize == 115240) { 
              dimsize = c(67,43,40)
              cat('Assuming gene expression data is from E18.5 brains\n')
           } else if (flsize == 806520) { 
              dimsize = c(94, 132, 65)
              cat('Assuming gene expression data is from E15.5 brains\n')
           } else if (flsize == 669369) { 
              dimsize = c(89,109,69)
              cat('Assuming gene expression data is from E13.5 brains\n')
           } else if (flsize == 210000) { 
              dimsize = c(70,75,40)
              cat('Assuming gene expression data is from E11.5 brains\n')
           } else { 
              stop(
                 "I can't tell what age this data if from.",
                 "I suggest giving the zip file downloaded directly from",
                 "the Allen Brain Institute as the input (filename parameter)",
                 "Or giving the download URL (with the url=T flag)."
                  )
           }
           return(dimsize)
}



# TODO: delete this function
rra.info.from.header = function(header) {
           header    = readLines(header)
           dimsize   = as.integer(
                        strsplit(
                           gsub(
                              '^DimSize = ','',
                              grep('^DimSize = ', header,value=T)
                              ),split = ' ')[[1]])
           spacing   = as.numeric(
                        strsplit(
                           gsub(
                              '^ElementSpacing = ','',
                              grep('^ElementSpacing = ', header,value=T)
                              ),split = ' ')[[1]])
           origin    = as.numeric(
                        strsplit(
                           gsub(
                              '^CenterOfRotation = ','',
                              grep('^CenterOfRotation = ', header,value=T)
                              ),split = ' ')[[1]])
           # I don't know if allen like row-major or column major orientation
           #  hopefully they only have symmetric matrices as directions...
           direction = as.numeric(
                        strsplit(
                           gsub(
                              '^TransformMatrix = ','',
                              grep('^TransformMatrix = ', header,value=T)
                              ),split = ' ')[[1]])
           n = sqrt(length(direction))
           if (round(n) != n) {stop('Direction Matrix is not square')}
           direction = matrix(direction, nrow=n)

           elementtype = gsub( '^ElementType = ', '', 
                  grep('^ElementType = ', header,value=T))
           if        ( elementtype == 'MET_UCHAR' ) { 
                size = 1
                atlas_or_annotations = 'atlas'
                pixeltype = 'unsigned char'
           } else if ( elementtype == 'MET_UINT'  ) {
                size = 4
                atlas_or_annotations = 'annotations'
                pixeltype = 'unsigned int'
           } else { stop('Elements must be either integers or floats') }

           # parameters for ants
           antsparam = list(
                             imagesize  = dimsize   ,
                             pixeltype  = pixeltype , 
                             components = F         ,
                             spacing    = spacing   , 
                             origin     = origin    ,
                             direction  = direction
                        )

           return(list(
             dimsize = dimsize, 
             size = size, 
             atlas_or_annotations = atlas_or_annotations,
             antsparam = antsparam
            ))
}

# TODO: delete this function
rrg.info.from.header = function(header) {
           header = readLines(header)
           dimsize = as.integer(
                        strsplit(
                           gsub(
                              '^DimSize = ','',
                              grep('^DimSize = ', header,value=T)
                              ),split = ' ')[[1]])
           elementtype = gsub( '^ElementType = ', '', 
                  grep('^ElementType = ', header,value=T))
           if        ( elementtype == 'MET_FLOAT' ) { typ = 'numeric'
           } else if ( elementtype == 'MET_UINT'  ) { typ = 'integer'
           } else { stop('Elements must be either integers or floats') }
           return(list(dimsize = dimsize, typ = typ))
}

