#' Run Gene Ontology enrichment analysis on GOrilla
#' @param target.genes genes of interest
#' @param background.genes pool of genes from which the target genes were selected
#' @param tmpdir path to directory where files will be stored
#' @return data frame of enriched process ontologies, significance from hypergeometric test (with FDR), and associated genes
#' @examples
#' # Target list of genes
#' target.genes = readLines(system.file('extdata/GOtest/target.txt',package="ABIgeneRMINC"))
#' # Background genes
#' background.genes = readLines(system.file('extdata/GOtest/background.txt',package="ABIgeneRMINC"))
#' 
#' # Run GOrilla
#' go_results = GOrilla_run(target.genes, background.genes)
#' @export

GOrilla_run = function(target.genes, background.genes, tmpdir = tempdir()) {
   # check if perl modules are installed
   check_perl_modules('WWW::Mechanize')
   check_perl_modules('Getopt::Long')

   # if directory does not exist, make it
   if (!file.exists(tmpdir)) {
      dir.create(tmpdir)
   }

   olddir = getwd()
   setwd(tmpdir)

   # copy necessary files
   file.copy(
     from = system.file('extdata/GOrilla.pl',package='ABIgeneRMINC'),
     to   = 'GOrilla.pl'
   )
   cat(target.genes, file = 'target.txt', sep='\n')
   cat(background.genes, file = 'background.txt', sep='\n')

   # run GOrilla
   results_url = system('perl GOrilla.pl', intern=T)[2]
   if (is.na(results_url)) {
      stop('GOrilla did not run')
   } else {
      cat('Results can be found here:','\n')
      cat(results_url,'\n')
   }

   # wait for results and download when available
   maxwait = 120     # 1 minute timeout
   i = 0 ; while (i<120) {
      i = i + 1
      keeptrying = tryCatch({
               xml2::download_html(
                       results_url,
                       file = 'GOrilla_web.html',
                       quiet = TRUE,
                       mode = "wb",
                       handle = curl::new_handle()
                )
             }, error = function(e) {
                      #cat('Waiting for results on attempt',i,'of', maxwait, '\n')
                      Sys.sleep(0.5)
                      return(TRUE)
             })
      if (keeptrying != TRUE) { break }
   }
   if (i == 120) {stop('Failed to get GOrilla results from the web')}

   # read results
   gorilla_res = xml2::read_html('GOrilla_web.html')
   GO_res = rvest::html_table(
             rvest::html_children(
              rvest::html_children(gorilla_res)[[2]])[[7]])

   # Remove First row and make it the header
   cnames = GO_res[1,]
   GO_res = GO_res[-1,]
   colnames(GO_res) = cnames

   # clean up gene list
   gene = pull(GO_res,'Genes')
   gene = lapply(gene, function(x) {
         strsplit(x,split='\r\n')[[1]][-1]
   })
   GO_res$Genes = gene
   
   # put into tibble if packase is installed
   if("tibble" %in% rownames(installed.packages())) {
      GO_res = tibble::as_tibble(GO_res)
   }

   # return to old directory
   setwd(olddir)

   return(GO_res)
}

check_perl_modules = function(module_name) {
   test = system(paste0('perl -e "use ',module_name,'"'), ignore.stdout = T, ignore.stderr = T)
   if (test == 2) {
      cat(paste0("Perl module '",module_name,"' is required"),"\n")
      cat("The easiest way to install it is to first install and run CPAN:","\n")
      cat("> ","\n")
      cat("> perl -MCPAN -e shell","\n")
      cat("> ","\n")
      cat("Then, install the missing module:","\n")
      cat("> ","\n")
      cat(paste0("> install module ",module_name),"\n")
      cat("> ","\n")
      stop(paste0(module_name," module is missing"))
   }
}


