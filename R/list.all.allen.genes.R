#' List all genes in the Allen Gene Expression Atlas
#'
#' @param outfile (OPTIONAL) Filename to output results
#' @param parallel Number of cores registered for parallel computation. 
#' Default is NULL, computations performed serially. 
#' If not NULL, requires the parallel package for computations.  
#'
#' @return Vector of all genes currently in the Allen Gene Expression Atlas
#'
#' @examples
#' # Find all genes in the Allen Gene Expression Atlas and
#' # Output results to file "genes.txt" and perform computations in parallel (4 cores).
#' \dontrun{
#' genes=list.all.allen.genes("genes.txt",parallel=4)
#' }
#' @export
list.all.allen.genes=function(outfile=NULL,parallel=NULL) {

  # use pbapply as progress bar if available
  if ("pbapply" %in% rownames(installed.packages())) {progb=T} else {probg=F}

  #find number of genes
  URL="http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Gene,rma::criteria,products[id$eq1],rma::options[start_row$eq0][num_rows$eq50]"
  xmlfile = xml2::read_xml(URL)
  total_rows = as.integer(xml2::xml_attrs(xmlfile)[['total_rows']])

  #API only displays max 50 rows. Repeatedly query to find all genes.
  start_row=seq(0,total_rows,50)
  URLS=paste0("http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Gene,rma::criteria,products[id$eq1],rma::options[start_row$eq",start_row,"][num_rows$eq50]")

  function_to_parallize = function(i) {
                URL=URLS[i]
                xmlfile = xml2::read_xml(URL)
                unlist(base::lapply(
                                xml2::as_list(
                                xml2::xml_children(
                                xml2::xml_children(
                                   xmlfile
                                ))), 
                                function(x) x[['acronym']][[1]]))
  }

  ni=length(URLS)
  if (!is.null(parallel)) {
    cl <- parallel::makeCluster(parallel)
    parallel::clusterExport(cl, c('URLS'), envir=environment())
    if (progb) {
       genes = unlist(pbapply::pblapply(cl=cl, 1:ni, function_to_parallize))
    } else {
       genes = unlist(parallel::parLapply(cl=cl, 1:ni, function_to_parallize))
    }
    parallel::stopCluster(cl)
  } else {
    if (progb) {
       genes = unlist(pbapply::pblapply(1:ni, function_to_parallize))
    } else {
       genes = unlist(base::lapply(1:ni, function_to_parallize))
    }
  }

  if (!is.null(outfile)) { cat(genes,file=outfile,sep='\n',append=FALSE) }
  return(genes)
}

#' List all genes in the Allen Developmental Gene Expression Atlas
#'
#' @param outfile (OPTIONAL) Filename to output results
#' @param parallel Number of cores registered for parallel computation. 
#' Default is NULL, computations performed serially. 
#' If not NULL, requires the parallel package for computations.  
#'
#' @return Vector of all genes currently in the Allen Developmental Gene Expression Atlas
#'
#' @examples
#' # Find all genes in the Allen Gene Expression Atlas and
#' # Output results to file "genes.txt" and perform computations in parallel (4 cores).
#' \dontrun{
#' genes=list.all.dev.allen.genes("genes.txt",parallel=4)
#' }
#' @export
list.all.dev.allen.genes = function(outfile=NULL,parallel=NULL) {
#URL="http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Gene,rma::criteria,products[id$eq1],rma::options[start_row$eq0][num_rows$eq50]"

  # use pbapply as progress bar if available
  if ("pbapply" %in% rownames(installed.packages())) {progb=T} else {probg=F}

  #find number of genes
   URL = "http://api.brain-map.org/api/v2/data/Gene/query.xml?criteria=products[id$eq3],rma::options[start_row$eq0][num_rows$eq50]"
   xmlfile = xml2::read_xml(URL)
   total_rows = as.integer(xml2::xml_attrs(xmlfile)[['total_rows']])
   
   #API only displays max 50 rows. Repeatedly query to find all genes.
   start_row=seq(0,total_rows,50)
   URLS=paste0("http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Gene,rma::criteria,products[id$eq3],rma::options[start_row$eq",start_row,"][num_rows$eq50]")

   function_to_parallize = function(i) {
                 URL=URLS[i]
                 xmlfile = xml2::read_xml(URL)
                 unlist(base::lapply(
                                 xml2::as_list(
                                 xml2::xml_children(
                                 xml2::xml_children(
                                    xmlfile
                                 ))), 
                                 function(x) x[['acronym']][[1]]))
   }

  ni=length(URLS)
  if (!is.null(parallel)) {
    cl <- parallel::makeCluster(parallel)
    parallel::clusterExport(cl, c('URLS'), envir=environment())
    if (progb) {
       genes = unlist(pbapply::pblapply(cl=cl, 1:ni, function_to_parallize))
    } else {
       genes = unlist(parallel::parLapply(cl=cl, 1:ni, function_to_parallize))
    }
    parallel::stopCluster(cl)
  } else {
    if (progb) {
       genes = unlist(pbapply::pblapply(1:ni, function_to_parallize))
    } else {
       genes = unlist(base::lapply(1:ni, function_to_parallize))
    }
  }

  if (!is.null(outfile)) { cat(genes,file=outfile,sep='\n',append=FALSE) }
  return(genes)
}

#' List the Allen adult mouse experiments associated with gene
#'
#' @param geneAcronym Acronym of gene to find experiments of
#' @return Data frame with Experiment number, slices (coronal or sagittal), and URLs to download data
#'
#' @examples
#' # Find all experiments associated with Itsn1
#' find.gene.experiment('Itsn1')
#' @export
find.gene.experiment=function(geneAcronym) {

  URL=paste0("http://api.brain-map.org/api/v2/data/SectionDataSet/query.xml?criteria=products[id$eq1],genes[acronym$eq'",geneAcronym,"']")
  xmlfile = xml2::read_xml(URL)
  
  ni   = as.integer(xml2::xml_attrs(xmlfile)[['num_rows']])
  ntot = as.integer(xml2::xml_attrs(xmlfile)[['total_rows']])
  if (ni != ntot) {
      warning('check this URL. It looks like it is not displaying as much data sets as it should:', URL)
      if (ntot>50) {
         cat('ABI API only displays the first 50 experiments. There are more you should consider')
      }
  }
  
  if (ni==0) {
     cat(paste('no experiments with gene acronym',geneAcronym,'\n'))
     return (NA)
  } else {
     sids=c() ; slices=c(); counter=0
     for (i in 1:ni) {
      if (xml2::as_list(xml2::xml_children(xml2::xml_children(xmlfile))[[i]])[['failed']][[1]]=='true') {next} #If experiment failed, skip
      counter=counter+1
      sids[counter]=as.integer(xml2::as_list(xml2::xml_children(xml2::xml_children(xmlfile))[[i]])[['id']][[1]])
       if (as.integer(xml2::as_list(xml2::xml_children(xml2::xml_children(xmlfile))[[i]])[['plane-of-section-id']][[1]])==1) {
         slices[counter]='coronal'
       } else {slices[counter]='sagittal'}
     }
     URLs=paste0("http://api.brain-map.org/grid_data/download/",sids,"?include=energy")
     
     df=data.frame(stringsAsFactors=F, gene=rep(geneAcronym,length(sids)),slices,ExperimentID=sids,URLs)
  
     return(df)
   }
}

#' List the Allen mouse development experiments associated with gene
#'
#' @param geneAcronym Acronym of gene to find experiments of
#' @return Data frame with Experiment number, slices (coronal or sagittal), and URLs to download data
#'
#' @examples
#' # Find all mouse development experiments associated with Ntng1
#' find.dev.gene.experiment('Ntng1')
#' @export
find.dev.gene.experiment=function(geneAcronym) {

  URL=paste0("http://api.brain-map.org/api/v2/data/SectionDataSet/query.xml?criteria=[failed$eqfalse],products[id$eq3],genes[acronym$eq'",geneAcronym,"']&include=specimen(donor(age))")
  xmlfile = xml2::read_xml(URL)
  
  ni   = as.integer(xml2::xml_attrs(xmlfile)[['num_rows']])
  ntot = as.integer(xml2::xml_attrs(xmlfile)[['total_rows']])
  if (ni != ntot) {
      warning('check this URL. It looks like it is not displaying as much data sets as it should:', URL)
      if (ntot>50) {
         cat('ABI API only displays the first 50 experiments. There are more you should consider')
      }
  }
  if (ni==0) {
     cat(paste('no development experiments with gene acronym',geneAcronym,'\n'))
     return (NA)
  } else {
    sids=c() ; slices=c(); age=c(); counter=0
    for (i in 1:ni) {
     if (xml2::as_list(xml2::xml_children(xml2::xml_children(xmlfile))[[i]])[['failed']][[1]]=='true') {next} #If experiment failed, skip
     counter=counter+1
     sids[counter]=as.integer(xml2::as_list(xml2::xml_children(xml2::xml_children(xmlfile))[[i]])[['id']][[1]])
      if (as.integer(xml2::as_list(xml2::xml_children(xml2::xml_children(xmlfile))[[i]])[['plane-of-section-id']][[1]])==1) {
        slices[counter]='coronal'
      } else {slices[counter]='sagittal'}
      age[counter] = xml2::as_list(xml2::xml_children(xml2::xml_children(xml2::xml_children(xml2::xml_children(xml2::xml_children(xml2::xml_children(xmlfile))[[i]])))))[[6]][[1]]
    }
    URLs=paste0("http://api.brain-map.org/grid_data/download/",sids,"?include=energy")
    
    df=data.frame(stringsAsFactors=F, gene=rep(geneAcronym,length(sids)),age=age,ExperimentID=sids,slices,URLs)
    fage = factor(df$age, levels = c("E11.5", "E13.5", "E15.5", "E18.5", "P4", "P14", "P28", "P56", "18M" , "24M"))
    df = df[order(fage),]
    return(df)
   }
}




# TODO: delete this
old.list.all.allen.genes=function(outfile=NULL,parallel=NULL) {
  library(XML)
  if (!is.null(parallel)) {
    library(foreach)
    library(doMC)
    registerDoMC(parallel)
  }
  #find number of genes
  URL="http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Gene,rma::criteria,products[id$eq1],rma::options[start_row$eq0][num_rows$eq50]"
  xmlfile=xmlRoot(xmlTreeParse(URL))
  total_rows=as.integer(xmlAttrs(xmlfile)[['total_rows']])

  #API only displays max 50 rows. Repeatedly query to find all genes.
  start_row=seq(0,total_rows,50)
  URLS=paste0("http://api.brain-map.org/api/v2/data/query.xml?criteria=model::Gene,rma::criteria,products[id$eq1],rma::options[start_row$eq",start_row,"][num_rows$eq50]")

  if (!is.null(parallel)) {
    ni=length(start_row);genes=c()
    out=foreach(i=1:ni,.errorhandling='stop') %dopar% {
      print(paste("still",total_rows-start_row[i],"Acronyms to get"))
      URL=URLS[i]
      xmlfile=xmlRoot(xmlTreeParse(URL))
      xmlobjs=xmlfile[[1]];nj=length(xmlobjs)
      gene=rep(NA,nj);for (j in 1:nj) { gene[j]=xmlValue(xmlobjs[[j]][["acronym"]]) }
      gene
     }
    genes=unlist(out)
  } else {
    ni=length(start_row);genes=c()
    for (i in 1:ni) {
     print(paste("still",total_rows-start_row[i],"Acronyms to get"))
     URL=URLS[i]
     xmlfile=xmlRoot(xmlTreeParse(URL))
     xmlobjs=xmlfile[[1]];nj=length(xmlobjs)
     gene=rep(NA,nj);for (j in 1:nj) { gene[j]=xmlValue(xmlobjs[[j]][["acronym"]]) }
     genes=c(genes,gene)
    }
  }
  if (!is.null(outfile)) { cat(genes,file=outfile,sep='\n',append=FALSE) }
  return(genes)
  }

