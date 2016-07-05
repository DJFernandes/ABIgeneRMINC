#' List the Allen experiments assiciated with gene
#'
#' @param geneAcronym Acronym of gene to find experiments of
#' @return Data frame with Experiment number, slices (coronal or sagittal), and URLs to download data
#'
#' @examples
#' # Find all experiments associated with Itsn1
#' find.gene.experiment('Itsn1')
#' @export



find.gene.experiment=function(geneAcronym) {
  library(XML)

  URL=paste0("http://api.brain-map.org/api/v2/data/SectionDataSet/query.xml?criteria=products[id$eq1],genes[acronym$eq'",geneAcronym,"']")
  xmlfile=xmlRoot(xmlTreeParse(URL))
  
  ni=xmlAttrs(xmlfile)[['total_rows']]
  if (ni==0) {stop(paste('no experiments with gene acronym',geneAcronym))}
  sids=c() ; slices=c(); counter=0
  for (i in 1:ni) {
   if (xmlValue(xmlfile[['section-data-sets']][[i]][['failed']])=='true') {next} #If experiment failed, skip
   counter=counter+1
   sids[counter]=xmlValue(xmlfile[['section-data-sets']][[i]][['id']])
    if (xmlValue(xmlfile[['section-data-sets']][[i]][["plane-of-section-id"]])==1) {
      slices[counter]='coronal'
    } else {slices[counter]='sagittal'}
  }
  URLs=paste0("http://api.brain-map.org/grid_data/download/",sids,"?include=energy")
  
  df=data.frame(gene=rep(geneAcronym,length(sids)),slices,ExperimentID=sids,URLs)
  
  return(df)
}



