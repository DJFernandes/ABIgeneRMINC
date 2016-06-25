#' @title Preferential Expresion in target
#' @description
#' This function will mask a statistic vector into a target and contrast region (default, contrast is the whole vector). Then, the mean gene expression energy in both target and contrast region is calculated, and its quotient (preferential gene expression in target vs contrast) is returned. 
#' @param stats 1-D vector or MINC filename whose values will define the target and contrast regions. It is recommended that you provide a 1D vector instead of a filename, though the code will work regardless. If you are providing a filename, make sure it is in MINC orientation (X=Left-to-Right, Y=Posterior-to-Anterior, Z=Inferior-to-Superior). This is usually the case with files produced with MINC Tools.
#' @param gene 1-D vector (same length as that of statsvector) or filename of gene expression data (raw format). It is recommended that you provide a 1D vector instead of a filename, though the code will work regardless. If you are providing a filename, make sure it is in ABI orientation X=Anterior-to-Posterior, Y=Superior-to-Inferior, and Z=Left-to-Right. This is usually the case when you download Allen Data.
#' @param maskvector mask for both stats and gene, for which elements to analyze. If filenames are provided for stats and gene, mask must be in MINC orientation. Vector Elements are either TRUE or FALSE and vector length must be the same as length of statsvector. DEFAULT: all elements are TRUE, ie. all elements in statsvector and gene are analyzed.
#' @param tgt.thresh threshold for target regions (statsvector>tgt.thresh form the target region)
#' @param cntrst.thresh threshold for contrast regions (statsvector<=crst.thresh form the contrast region). DEFAULT: contrast includes all the elements in the mask
#' @return Fold-Change between mean expression in the target and contrast regions
#'
#' @examples
#' # I included several files in this package for this example: 
#' #    The raw gene expression energy files for Nrp1 (http://api.brain-map.org/grid_data/download/74272479)
#' nrp1Filename=system.file('extdata/Nrp1_energy_experiment74272479.raw',package="ABIgeneRMINC")
#' #    Allen Brain Institute annotations (http://download.alleninstitute.org/informatics-archive/current-release/mouse_annotation/P56_Mouse_gridAnnotation.zip)
#' annotFilename=system.file('extdata/gridAnnotation.raw',package="ABIgeneRMINC")
#' #    T-statistics mincfiles comparing neuroanatomical volumes of mice raised in enriched environments vs standard lab cages.
#' statFilename=system.file('extdata/enrichment_stats.mnc',package="ABIgeneRMINC")
#' # 
#' # The following example calculates the foldchange between Nrp1 expression 
#' # in significantly larger neuroanatomy (t-statistics > 2) verses the whole brain.
#' # Read gene file and rotate it to MINC orientation
#' nrp1expr=allenVectorTOmincVector(read.raw.gene(nrp1Filename))
#' 
#' # Read stats file
#' stats=mincGetVolume(statFilename)
#' 
#' # Read annotations, rotate it to MINC orientation, and binarize to make a brain maskvector
#' annotations=allenVectorTOmincVector(read.raw.gene(annotFilename,labels=TRUE))
#' maskvector=annotations>0
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
#' geneFoldChange(statFilename,nrp1Filename,maskvector,tgt.thresh=2)
#' # 2.072129
#' @export


geneFoldChange=function(stats,gene,maskvector=NULL,tgt.thresh,cntrst.thresh=NULL) {
        ## if input is a file read it.
        if (is.character(gene)) {gene.vector=allenVectorTOmincVector(read.raw.gene(gene))} else {gene.vector=gene}
        if (is.character(stats)) {statsvector=mincGetVolume(stats)} else {statsvector=stats}
        ## If maskvector is not provided, all voxels are valid
	if (length(maskvector)==0) {maskvector=rep(TRUE,length(gene.vector))}
       ## Mask out negative values (Allen gene expression does not exists there)
	mask=(!gene.vector<0)&maskvector
       ## Mask the statsvector and gene.vector
	stats.masked=statsvector[mask];gene.masked=gene.vector[mask]
       ## If crst.thresh is not given, contrast is all the voxels
	if (length(cntrst.thresh)==0) {fc=mean(gene.masked[stats.masked>tgt.thresh])/mean(gene.masked)
	} else {
        fc=mean(gene.masked[stats.masked>tgt.thresh])/mean(gene.masked[stats.masked<=cntrst.thresh])}
	return(fc)
}

