#' @title Find Ontology Structure Descendants 
#' @description
#' Given a stucture label in the Allen Brain Atlas, this function will find all its descendants structures in the the grid annotations
#' 
#' @param parents Label number of Parent Structure. This parameter is vectorized, so you can specify multiply parents. The return will be a list of same length as the number of parents.
#'
#' @return The label number of all descendant structures. These structures can be found in the 200 um grid labels for gene expression analysis
#'
#' @examples
#' # In this example, we will find the gene expression of Pdyn in the Hippocampal formation.
#' # First we download the gene expression file from the Allen Brain Institute
#' # In this example, I downloaded Pdyn expression energy (http://api.brain-map.org/grid_data/download/71717084)
#' # Then Unzipped it, and renamed the energy.raw file to Pdyn_energy_experiment71717084.raw 
#' # I included it in this library as an example
#' 
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
#' # The hippocampal formation has many Ontology Structure Descendants. We need to find them.
#' # Find hippocampal formation structure number
#' data(AllStructureLabels)
#' parent.label=subset(AllStructureLabels,name=='Hippocampal formation')$id   #the sturcture number is 1089
#' 
#' # Find Descendant labels
#' child.labels=ont.children(parent.label)[[1]]
#' 
#' # Unionize over all the children
#' unionize(gene.expression,child.labels,labels.grid)
#' 
#' @export


ont.children <- function(parents) {
	np=length(parents)
	data(AllStructureLabels)
	parentdf=AllStructureLabels[,c("id","structure_id_path")]
	data(GridStructureLabels)
	contdf=GridStructureLabels[,c("id","structure_id_path")]
	children=vector("list", np)
	for (i in 1:np) {
		parent=parents[i]
		parentaddress=parentdf[parentdf$id==parent,"structure_id_path"]
		children[[i]]=contdf[grepl(parentaddress,contdf$structure_id_path),"id"]
	}
	return(children)
}
