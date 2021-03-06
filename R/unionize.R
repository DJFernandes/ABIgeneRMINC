#' @title Sum Grid Values over Labels
#' @description
#' Using this function you can add the gene expression voxels within a particular label. This allows you to find gene expression of any gene in a particular structure.
#' 
#' @param grid.data 1D vector of voxels to sum over (fro example gene expression data)
#' @param labels.to.sum The struction labels you want to unionize over. This is parameter is vectorized so you can sum over multiple labels and get statistical data for each label seperately
#' @param labels.grid 1D vector of labels (must be registered to grid.data)
#' @param neg.rm Remove negative data for unionization (negative gene expression means there is no data there)
#'
#' @return data frame with labels as rows; and columns with sum, mean, and Standard Deviation of unionized grid.data
#'
#' @examples
#' # In this example, we will find the gene expression of Pdyn in the Paraventricular Nucleus of the Thalamus.
#' # 
#' # First we download the gene expression file from the Allen Brain Institute
#' # In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
#' # Then Unzipped it, and renamed the energy.raw file to Pdyn_energy_experiment71717084.raw 
#' # I included it in this library as an example

#' # Also, download Grid Labels from the Allen Brain Institute
#' # In this example, I downloaded Grid Labels from (http://download.alleninstitute.org/informatics-archive/current-release/mouse_annotation/P56_Mouse_gridAnnotation.zip), then unzipped it.
#' # I included it in this library as an example though it is a good idea to download it yourself to keep it updated
#'
#' # Read Gene expression file
#' genefilename=system.file('extdata/Pdyn_energy_experiment71717084.raw',package="ABIgeneRMINC")
#' gene.expression=read.raw.gene(genefilename)
#' 
#' # Read Labels 
#' labelfilename=system.file('extdata/gridAnnotation.raw',package="ABIgeneRMINC")
#' labels.grid=read.raw.gene(labelfilename,labels=TRUE)
#' 
#' # Find Label Number Corresponding to Paraventricular Nucleus of the Thalamus
#' data(GridStructureLabels)
#' label.num=subset(GridStructureLabels,name=='Paraventricular nucleus of the thalamus')$id
#' # Paraventricular nucleus of the thalamus corresponds to label 149 in labels.grid
#' 
#' # Unionize to find gene expression in the Paraventricular nucleus of the thalamus
#' unionize(gene.expression,label.num,labels.grid)
#'
#' @export


unionize <- function(grid.data,labels.to.sum,labels.grid,neg.rm=TRUE) {
	domain=labels.to.sum
	summand=grid.data
	labelsvec=labels.grid
	df=as.data.frame(matrix(NA,nrow=length(domain),ncol=4))
	colnames(df)=c("labels","sum","mean","stdev")
	df[,"labels"]=domain
	for (i in 1:length(df[,1])) {
		bool=(labelsvec==domain[i])
		sbool=sum(bool)
		if (sbool==0) { df[i,c("sum","mean","stdev")]=NA } else {
			vec=summand[bool]
#ignore negative values
	if (neg.rm) {	vec=vec[!vec<0] }
			df[i,"sum"]=sum(vec)
			df[i,"mean"]=df[i,"sum"]/sbool
			df[i,"stdev"]=sd(vec)}
	}
	return(df)
}

