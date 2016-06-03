#' List all genes in the Allen Gene Expression Atlas
#'
#' @param outfile (OPTIONAL) Filename to output results
#' @param parallel Number of cores registered for parallel computation. 
#' Default is NULL, computations performed serially. 
#' If not NULL, requires the foreach & doMC pakages for parallel computations.  
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

