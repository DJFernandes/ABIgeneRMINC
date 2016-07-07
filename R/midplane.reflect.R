#' Reflect Gene Expression Vector across midplane
#'
#' @param gene.expression.vector Vector of gene expression to reflect
#' @param reflect.dimension INTEGER. Which dimension to reflect. 1 for the slowest dimension, 3 for the fastest. To reflect Allen orientation vectors across the sagittal midplane, set this value to 1. To reflect MINC orientation vectors across the sagittal midplane, set this value to 3. 
#' @param array.dimensions Vector of length 3 indicating the dimension sizes in order of the SLOWEST dimension to the FASTEST. Must be supplied if gene.expression.vector is missing 'sizes' attributes (execute attr(gene.expression.vector,'sizes')).
#'
#' @return 1D vector of length 159326 containing Raw Gene Expression data
#'
#' @examples
#' # Read Bdnf gene expression (Experiment 75695642) from Allen Brain Institute
#' # Look at a coronal slice
#' # Reflect across the sagittal midplane 
#' # Look at a coronal slice
#'
#'
#' ## Read Bdnf expression (Experiment 75695642)
#' URLtoDownload='http://api.brain-map.org/grid_data/download/75695642?include=energy'
#' gene.expression.vector=read.raw.gene(URLtoDownload,url=TRUE)
#' 
#' ## Coronal slice
#' c.slice=array(gene.expression.vector,rev(attr(gene.expression.vector,'sizes')))[30,,]
#' image(c.slice, ylab='Left-Right' ,xlab='Superior-Inferior')
#' 
#' ## Reflect across the sagittal midplane 
#' reflected.vector=midplane.reflect(gene.expression.vector,1)
#' 
#' ## Coronal slice
#' reflected.c.slice=array(reflected.vector,rev(attr(gene.expression.vector,'sizes')))[30,,]
#' image(reflected.c.slice, ylab='Left-Right' ,xlab='Superior-Inferior')
#'
#' ## Reflect MINC orientation vector across sagittal midplane
#' minc.vector=allenVectorTOmincVector(gene.expression.vector)
#' reflected.minc.vector=midplane.reflect(minc.vector,3)
#' @export

midplane.reflect=function(gene.expression.vector,reflect.dim,array.dimensions=NULL){
     if (is.null(attr(gene.expression.vector,'sizes'))) {
     if (is.null(array.dimensions) | length(array.dimensions) != 3) { 
        stop(paste("no attribute sizes in gene.expression.vector. Must provide attribute OR",
        "array.dimensions vector of length 3"))
     } else {arrdim=rev(array.dimensions)}
     } else {arrdim=rev(attr(gene.expression.vector,'sizes'))}
     
     # put into an array. reverse sizes attribute for array dimensions
     gene.expression.array=array(gene.expression.vector,arrdim)
     # must reverse dimension to reflect
     # reflect.dim=1 --> rev.reflect.dim=3
     # reflect.dim=2 --> rev.reflect.dim=2
     # reflect.dim=3 --> rev.reflect.dim=1
     rev.reflect.dim=4-reflect.dim
     
     # find midplane
     ref.dim.halflength=floor(median(1:arrdim[rev.reflect.dim]))
     
     #Find indices to split array (depending on which dimension to reflect)
     #define index order to be closest to the mirror
     if (rev.reflect.dim==3) {
          idx1=paste0(",,",ref.dim.halflength,":1")
          idx2=paste0(",,(",ref.dim.halflength,"+1):",arrdim[rev.reflect.dim])
     
          } else if (rev.reflect.dim==2) {
          idx1=paste0(",",ref.dim.halflength,":1,")
          idx2=paste0(",(",ref.dim.halflength,"+1):",arrdim[rev.reflect.dim],",")
     
          } else if (rev.reflect.dim==1) {
          idx1=paste0(ref.dim.halflength,":1,,")
          idx2=paste0("(",ref.dim.halflength,"+1):",arrdim[rev.reflect.dim],",,")}
     
     
     genearr1=eval(parse(text=paste0("gene.expression.array[",idx1,"]")))
     genearr2=eval(parse(text=paste0("gene.expression.array[",idx2,"]")))
     
     ## Find where one half is missing data (missing data are elements with -1)
     ## Replace it with the same voxels in the other half
     bool=genearr1==-1
     genearr1[bool]=genearr2[bool]
     
     ## DO the same with the other half
     bool=genearr2==-1
     genearr2[bool]=genearr1[bool]
     
     
     ## Put reflected back into a return array
     ret.array=gene.expression.array*0
     eval(parse(text=paste0("ret.array[",idx1,"]=genearr1")))
     eval(parse(text=paste0("ret.array[",idx2,"]=genearr2")))
     
     ## Output as vector with the old attributes
     ret.vector=as.vector(ret.array)
     attr(ret.vector,'sizes')=rev(arrdim)
     
     return(ret.vector)
}
