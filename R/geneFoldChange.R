#' @title Preferential Expresion in target
#' @description
#' This function will mask a statistic vector into a target and contrast region (default, contrast is the whole vector). Then, the mean gene expression energy in both target and contrast region is calculated, and its quotient (preferential gene expression in target vs contrast) is returned. 
#' @param stats 1-D vector or MINC filename whose values will define the target and contrast regions. It is recommended that you provide a 1D vector instead of a filename, though the code will work regardless. If you are providing a filename, make sure it is in MINC orientation (X=Left-to-Right, Y=Posterior-to-Anterior, Z=Inferior-to-Superior). This is usually the case with files produced with MINC Tools.
#' @param gene 1-D vector (same length as that of statsvector) or filename of gene expression data (raw format). It is recommended that you provide a 1D vector instead of a filename, though the code will work regardless. If you are providing a filename, make sure it is in ABI orientation X=Anterior-to-Posterior, Y=Superior-to-Inferior, and Z=Left-to-Right. This is usually the case when you download Allen Data.
#' @param maskvector mask for both stats and gene, for which elements to analyze. If filenames are provided for stats and gene, mask must be in MINC orientation. Vector Elements are either TRUE or FALSE and vector length must be the same as length of statsvector. DEFAULT: all elements are TRUE, ie. all elements in statsvector and gene are analyzed.
#' @param tgt.thresh threshold for target regions (statsvector>tgt.thresh form the target region)
#' @param cntrst.thresh threshold for contrast regions (statsvector<=crst.thresh form the contrast region). DEFAULT: contrast includes all the elements in the mask
#' @param checkOrientation should I check orientation. You should only turn this off if you know for sure that stats, gene, and mask are in the same orientation (either all ABI or all MINC).
#' @return Fold-Change between mean expression in the target and contrast regions
#'
#' @examples
#' # I included several files in this package for this example: 
#' #    The raw gene expression energy files for Nrp1 (http://api.brain-map.org/grid_data/download/74272479)
#' nrp1Filename=system.file('extdata/Nrp1_P56_coronal_74272479_200um.zip',package="ABIgeneRMINC")
#' #    Allen Brain Institute annotations (http://download.alleninstitute.org/informatics-archive/current-release/mouse_annotation/P56_Mouse_gridAnnotation.zip)
#' annotFilename=system.file('extdata/P56_Mouse_gridAnnotation.zip',package="ABIgeneRMINC")
#' #    T-statistics mincfiles comparing neuroanatomical volumes of mice raised in enriched environments vs standard lab cages.
#' statFilename=system.file('extdata/enrichment_stats.mnc',package="ABIgeneRMINC")
#' # 
#' # The following example calculates the foldchange between Nrp1 expression 
#' # in significantly larger neuroanatomy (t-statistics > 2) verses the whole brain.
#' # Read gene file and rotate it to MINC orientation
#' nrp1expr=allenVectorTOmincVector(read.raw.gene(nrp1Filename))
#' 
#' # Read stats file
#' stats=RMINC::mincGetVolume(statFilename)
#' 
#' # Read annotations, rotate it to MINC orientation, and binarize to make a brain maskvector
#' annotations=allenVectorTOmincVector(read.raw.gene(annotFilename))
#' maskvector=annotations>0 ; attributes(maskvector) = attributes(annotations)
#' 
#' # Foldchange of Nrp1 expression in target (regions with stats>2) vs contrast (whole brain)
#' geneFoldChange(stats,nrp1expr,maskvector,tgt.thresh=2)
#' # 2.072129
#' 
#' # Foldchange of Nrp1 expression in target (regions with stats>2) vs contrast (regions with stats<=0.5 )
#' geneFoldChange(stats,nrp1expr,maskvector,tgt.thresh=2,cntrst.thresh=0.5)
#' # 2.415599
#' # 
#' # Foldchange of Nrp1 expression in target (regions with stats>2) vs contrast (whole brain)
#' # using strings with filenames instead of 1-D vectors
#' geneFoldChange(statFilename,nrp1Filename,mask = annotations,tgt.thresh=2)
#' # 2.072129
#' @export
geneFoldChange=function(stats,gene,maskvector=NULL,mask=NULL,tgt.thresh,cntrst.thresh=NULL, checkOrientation = TRUE) {

       ## if input is a file read it.
        file_reads = gfg.filereads(stats, gene, mask)
        statsvector = file_reads[[1]]
        gene.vector = file_reads[[2]]
        if (!is.null(mask)) { maskvector = file_reads[[3]] }

       ## If maskvector or mask is not provided, all voxels are valid
	if (length(maskvector)==0) {maskvector=rep(TRUE,length(gene.vector))}
       
       ## Check if stats, genes, mask have the same orientation
       if (checkOrientation) {
           gfg.checkOrientation(statsvector, gene.vector, maskvector)
       }

       ## Mask out negative values (Allen gene expression does not exists there)
	genemask=(!gene.vector<0)&maskvector

       ## Mask the statsvector and gene.vector
	stats.masked=statsvector[genemask];gene.masked=gene.vector[genemask]

       ## If crst.thresh is not given, contrast is all the voxels
	if (length(cntrst.thresh)==0) {
           tgt.mn = mean(gene.masked[stats.masked>tgt.thresh])
           full.mn = mean(gene.masked)
           fc=tgt.mn/full.mn
	} else {
           tgt.mn = mean(gene.masked[stats.masked>tgt.thresh])
           crst.mn = mean(gene.masked[stats.masked<=cntrst.thresh])
           fc=tgt.mn/crst.mn
        }
	return(fc)
}

is.mncfile <- function(filepath){
  if ( is.zip( filepath ) ) { return (FALSE)
  } else {
    # hacky way to determine is something is a readable minc file
    result <- tryCatch({
                RMINC::mincGetVolume(filepath)
                return(TRUE)
              }, error = function(e){
                return(FALSE)
              })
    return(result)
   }
  return(result)
}

gfg.filereads = function(stats, gene, mask) {

        ## if input is a file read it.
        if (is.character(gene)) {
             gene.vector=allenVectorTOmincVector(read.raw.gene(gene))
        } else {gene.vector=gene}
        if (is.character(stats)) {
             statsvector=RMINC::mincGetVolume(stats)
        } else {statsvector=stats}
        if (!is.null(mask)) { 
           if (is.character(mask)) {
             if (is.mncfile(mask)) {
                maskvector = RMINC::mincGetVolume(mask) > 0.5
                attr(maskvector,'orientation') = 'MINC'
             } else {
                maskvector = allenVectorTOmincVector(read.raw.gene(mask, labels=TRUE)) > 0.5
                attr(maskvector,'orientation') = 'MINC'
             }
           } else if (RMINC::is.minc(mask)) {
                maskvector = mask
                if (is.null(attr(maskvector,'orientation'))) { 
                   attr(maskvector,'orientation') = 'MINC'
                }
           } else { 
              stop(
                    'mask should be either a minc file or ABI file (zip and raw accepted). ',
                    'If you want to give an arbitrary vector, supply it as the maskvector argument.'
                  ) 
	   }
           return(list(statsvector,gene.vector,maskvector))
        } else {
           return(list(statsvector,gene.vector))
        }

}

gfg.checkOrientation = function(statsvector, gene.vector, maskvector) {
       # Check if stats, genes, mask have the same orientation
        ## orientation of minc file
        mincfile_orientation = NA
        if ( !is.null(attr(statsvector,'orientation')) ) {
           mincfile_orientation = attr(statsvector,'orientation')
        } else if (RMINC::is.minc(statsvector)) {
           mincfile_orientation = 'MINC'
        }

        ## orientation of allen file
        allen_orientation = attr(gene.vector,'orientation')

        ## mask orientation if provided
        mask_orientation = NA
        if ( !is.null(attr(maskvector,'orientation')) ) {
           mask_orientation = attr(maskvector,'orientation')
        } 

        if (is.na(mincfile_orientation)) {warning('Could not determine orientation of stats')}
        if (is.na(allen_orientation)) {warning('Could not determine orientation of gene')}
        if (is.na(mask_orientation)) {warning('Could not determine orientation of maskvector')}

        if ( !is.na(mincfile_orientation) & !is.na(allen_orientation) ) {
           if (mincfile_orientation != allen_orientation) {
              cat('stats is in',mincfile_orientation,'orientation','\n')
              cat('gene is in',allen_orientation,'orientation','\n')
              warning('stats and gene are not in the same orientation')
           }
        }
        if (!is.na(mincfile_orientation) & !is.na(mask_orientation)) {
           if (mincfile_orientation != mask_orientation) {
              cat('stats is in',mincfile_orientation,'orientation','\n')
              cat('mask is in',mask_orientation,'orientation','\n')
              warning('stats and mask are not in the same orientation')
           }
        }
        if (!is.na(allen_orientation) & !is.na(mask_orientation)) {
           if (allen_orientation != mask_orientation) {
              cat('gene is in',allen_orientation,'orientation','\n')
              cat('mask is in',mask_orientation,'orientation','\n')
              warning('gene and mask are not in the same orientation')
           }
        }
        return(invisible(NULL))
}

is.minc.silent = function(filename) {
   is.minc(filename)
}
