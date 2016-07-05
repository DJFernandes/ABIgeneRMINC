#' Read Raw Allen Gene Expression File
#'
#' @param filename Filename (or URL) of Raw Gene Expression File
#' @param labels logical. If TRUE, then file is read as a string of integers
#' @param url logical. If TRUE, then filename MUST be a valid URL pointing to the allen gene expression data. The zip file will be downloaded into tmp and the gene expression data extracted and read. It is recommended you download and unzip the file yourself rather than relying on this flag. Debugging becomes much harder when you use this flag.
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

read.raw.gene <- function(filename,labels=FALSE,url=FALSE) {
        # check if filename exists
       	if (url==TRUE) { 
          library(RCurl)
          if (!url.exists(filename)) {stop(paste(filename,"does not exist"))}
        } else {
          if (!file.exists(filename)) { stop(paste(filename,"does not exist"))}
        }
        # open connection
       	if (url==TRUE) {
           temp=tempfile()
           download.file(filename,temp)
           connec=unz(temp, "energy.raw",'rb')
         } else {
           connec=file(filename, "rb")
         }
	#dimensions of allen gene expression arrays are 67*41*58=159326
	if (labels) { typ="integer"} else {typ="numeric"}
	dat=readBin(connec,typ,size=4,n=159326)
        #close connections
	close(connec)
        if (url==TRUE) { unlink(temp) }
        #set sizes attributes
        attr(dat,'sizes')=c(58,41,67)
	return(dat)
}

