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
  sids=rep(NA,ni) ; slices=rep(NA,ni)
  for (i in 1:ni) {
   sids[i]=xmlValue(xmlfile[['section-data-sets']][[i]][['id']])
    if (xmlValue(xmlfile[['section-data-sets']][[i]][["plane-of-section-id"]])==1) {
      slices[i]='coronal'
    } else {slices[i]='sagittal'}
  }
  URLs=paste0("http://api.brain-map.org/grid_data/download/",sids,"?include=energy")
  
  df=data.frame(gene=rep(geneAcronym,ni),slices,ExperimentID=sids,URLs)
  
  return(df)
}



