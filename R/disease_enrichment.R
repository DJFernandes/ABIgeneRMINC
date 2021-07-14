#' gene to NCBI
#'
#' Converts gene acronym to NCBI gene IDs
#' @param genes gene acronyms
#  @param organism Which organism do there genes belong to? Only 'mouse' and 'human' is supported. 
#  @param updated_info Should data be re-downloaded? If FALSE, data stored in system file (downloaded June 14th, 2021) will be used. Only applies to mouse data.
#' @param parallel how many processors to run on (default=single processor). Specified as a two element vector, with the first element corresponding to the type of parallelization, and the second to the number of processors to use. For local running set the first element to "local" or "snowfall" for back-compatibility, anything else will be run with batchtools see pMincApply Leaving this argument NULL runs sequentially.
#' @param conf_file  A batchtools configuration file defaulting to getOption("RMINC_BATCH_CONF")
#' @return NCBI gene IDs
#' @importFrom geneSynonym mouseSyno
#' @importFrom parallel stopCluster makeCluster clusterExport
#' @export
gene_to_NCBI = function(
              genes, 
              organism = 'mouse',
              update_info = FALSE,
              parallel = NULL,
              conf_file = getOption("RMINC_BATCH_CONF")
    ) {

    if (length(organism) != 1) { 
       stop('Organism must be provided.',
            'Only "mouse" and "human" are allowed values') }
    if (! organism %in% c('mouse','human')) {
       cat('Only "mouse" and "human" are allowed values for organism')
       stop(organism, 'is not supported.' )
    }

    if (organism == 'mouse') {
       synoData = geneSynonym:::teval("geneSynonym::syno10090")
     } else {
       synoData = geneSynonym:::teval("geneSynonym::syno9606")
     }

    geneSearcher = function(x) {
        synonyms = strsplit(synoData[grepl(paste0("(^|[|])", 
            "\\Q", tolower(x), "\\E", "($|[|])"), tolower(synoData)) | 
            (names(synoData) %in% x)], split = "[|]")
        if (length(synonyms) != 0) {
           ret = tibble(gene=x, NCBI = as.integer(names(synonyms)))
        } else {
           ret = tibble(gene=x, NCBI = NA)
        }
        return(ret)
    }
    
    parfunc = function(group) bind_rows(base::lapply(group, geneSearcher))

    if (organism != 'mouse' & !update_info) {
       cat('update_info can only be set to FALSE if organism is "mouse"\n')
       warning('changing update_info to true\n')
       update_info = TRUE
    }
    
    if (!update_info) {
       col_types = cols(
          acronym = col_character(),
          cutting_plane = col_character(),
          Section_Data_ID = col_double(),
          URL = col_character(),
          entrez.id = col_double(),
          homologous_human_ID = col_character()
         )

       ncbi_ids = tibble(acronym=genes) %>% 
         inner_join(
             read_csv(
                  system.file('extdata/ABI_genes.csv',package="ABIgeneRMINC"),
                  col_types = col_types),
             by='acronym') %>% 
         rename(gene='acronym', NCBI='entrez.id') %>% 
         select(gene, NCBI)
    } else {
       if (is.null(parallel)) {
           if ("pbapply" %in% rownames(installed.packages())) {
                lapply = pbapply::pblapply
            }
           groups <- genes %>% as.list
           ncbi_ids = bind_rows(lapply( groups , parfunc ))
           lapply = base::lapply
       } 
       else {
           par_bool_local = parallel[1] %in% c("local", "snowfall")
           par_bool_pb = "pbapply" %in% rownames(installed.packages())
           n_groups <- as.numeric(parallel[2])
           if ( par_bool_local & par_bool_pb ) {
              cl <- makeCluster(n_groups)
              on.exit(stopCluster(cl))
              clusterExport(cl, c("bind_rows", "tibble"))
              ncbi_ids <- pbapply::pblapply(genes, geneSearcher, cl = cl) %>% 
                              bind_rows
           } 
           else if (par_bool_local) {
              if (length(genes) > 1) {
                  groups  <- split(
                     genes, 
                     RMINC:::groupingVector(length(genes), n_groups))
               } 
              else {
                  groups  <- list(genes)
               }
              ncbi_ids <- RMINC:::failing_mclapply(
                                     groups, 
                                     parfunc, 
                                     mc.cores = n_groups) %>% 
                              bind_rows
           } 
           else {
              if (length(genes) > 1) {
                  groups  <- split(
                     genes, 
                     RMINC:::groupingVector(length(genes), n_groups))
               } 
              else {
                  groups  <- list(genes)
               }
              regdir_name = RMINC:::new_file("geneSynonym_registry")
   
              cat('Setting up parallel jobs on cluster with the following configurations\n')
              cat(paste0('    Configuration file: ',conf_file,'\n'))
              cat(paste0('    Registry directory: ',getwd(),'/',regdir_name,'\n'))
              suppressMessages({
                 reg = makeRegistry(
                             regdir_name, 
                             packages = c('tidyverse','geneSynonym'), 
                             conf.file = conf_file
                           )
                 batchExport(reg=reg , list(
                          synoData = synoData
                    ))
               })
              on.exit({ 
                 cat('Cleaning up parallel jobs. Deleting jobs registry in 5 seconds...')
                 suppressMessages({
                    RMINC:::tenacious_remove_registry(reg) 
                 })
                 cat('Done\n')
               })
              suppressWarnings({suppressMessages({
                   ids <- batchMap(
                                reg = reg, 
                                fun = parfunc,
                                group = groups)
               })})
   #           suppressMessages({
                 submitJobs(ids, reg = reg)
                 waitForJobs(reg = reg)
   #            })
              not_done_ids = findNotDone(reg=reg)$job.id
              if (length(not_done_ids)>0) {
                 cat(length(not_done_ids),'have not been completed','\n')
                 cat('Here is the message in job',not_done_ids[1],':\n')
                 cat('------',tail(getLog(not_done_ids[1]),1),'------\n')
                 stop(paste0('Some ',parallel[1],' jobs have failed'))
              }
              suppressMessages({
                 ncbi_ids <- reduceResultsList(reg = reg) %>% bind_rows
               })
            }
       }
    }
    return(ncbi_ids)
}

#' Mouse to human genes
#'
#' Finds homologous human genes from mouse genes
#  @param mouse_genes Vector of mouse genes. If integers (or coercable to integers), assumed to be NCBI gene ids. If characters, assumed to be gene acronyms. 
#' @param parallel how many processors to run on (default=single processor). Only used if mouse_genes are characters. Specified as a two element vector, with the first element corresponding to the type of parallelization, and the second to the number of processors to use. For local running set the first element to "local" or "snowfall" for back-compatibility, anything else will be run with batchtools see pMincApply Leaving this argument NULL runs sequentially.
#' @param conf_file  A batchtools configuration file defaulting to getOption("RMINC_BATCH_CONF"). Only used if mouse_genes are characters. 
#' @return tibble with homologous genes
#' @import homologene
#' @importFrom parallel stopCluster makeCluster clusterExport
#' @export
mouse_genes_to_human_genes = function(
              mouse_genes, 
              parallel = NULL,
              conf_file = getOption("RMINC_BATCH_CONF")
    ) {
         # convert mouse_genes to integer of NCBI values
         convert_to_NCBI_ids_flag = tryCatch({
            x = as.integer(mouse_genes)
            FALSE
          },
          error = function(e) {
            stop(e)
          }, 
          warning = function(w) {
            TRUE
          }
         )
         if (convert_to_NCBI_ids_flag) {
            NCBI_df = gene_to_NCBI(
                           mouse_genes,
                           organism = 'mouse',
                           update_info=FALSE,
                           parallel=parallel, 
                           conf_file=conf_file)
            NCBI_ids = unique(na.omit(NCBI_df$NCBI))
          } else {
            NCBI_ids = as.integer(mouse_genes)
          }

         all_mouse_genes = homologene::homologeneData2 %>% 
                              as_tibble %>% filter(Taxonomy == '10090') %>% 
                              rename(mouse = 'Gene.Symbol', mouse_ID = 'Gene.ID') %>% 
                              select(HID, mouse, mouse_ID)
         all_human_genes = homologene::homologeneData2 %>% 
                              as_tibble %>% filter(Taxonomy == '9606') %>% 
                              rename(human = 'Gene.Symbol', human_ID = 'Gene.ID') %>% 
                              select(HID, human, human_ID)

         na_bool = is.na(NCBI_ids)
         ret = all_mouse_genes %>% 
                              filter( mouse_ID %in% NCBI_ids[!na_bool] ) %>% 
                              right_join( 
                                 tibble(mouse_ID = NCBI_ids[!na_bool]) , 
                                 by='mouse_ID' ) %>% 
                              left_join(all_human_genes, by='HID') %>% 
                              select(-HID) %>% 
                              right_join(
                                  tibble(
                                     idx = (1:length(NCBI_ids))[!na_bool], 
                                     mouse_ID = NCBI_ids[!na_bool]),
                                  by = 'mouse_ID'
                              ) %>% 
                              right_join(tibble(idx = (1:length(NCBI_ids))),by = 'idx') %>% 
                              arrange(idx) %>% select(-idx)

         if (convert_to_NCBI_ids_flag) {
            ret = ret %>% 
                     rename(NCBI = 'mouse_ID') %>% select(-mouse) %>% 
                     left_join(NCBI_df,., by = 'NCBI') %>% 
                     rename(mouse = 'gene', mouse_ID = 'NCBI')

         }

         return(ret)
}

#' Disease Ontology
#'
#' Finds DO of a gene set
#  @param genes Vector of genes that are of interest. If integers (or coercable to integers), assumed to be NCBI gene ids. If characters, assumed to be gene acronyms. 
#  @param organism Which organism do there genes belong to? If organism='human', disease enrichment is run directly on the genes provided. If organism='mouse', mouse genes will be converted to human genes before running disease enrichment analysis. Only 'mouse' and 'human' is supported. 
#' @param parallel how many processors to run on (default=single processor). Specified as a two element vector, with the first element corresponding to the type of parallelization, and the second to the number of processors to use. For local running set the first element to "local" or "snowfall" for back-compatibility, anything else will be run with batchtools see pMincApply Leaving this argument NULL runs sequentially.
#' @param conf_file  A batchtools configuration file defaulting to getOption("RMINC_BATCH_CONF").
#' @return tibble with homologous genes
#' @import disgenet2r
#' @export
disease_ontology = function(
              genes,
              organism = 'mouse',
              update_info = FALSE,
              api_key = NULL,
              parallel = NULL,
              conf_file = getOption("RMINC_BATCH_CONF")
    ) {
    
    if (length(organism) != 1) { 
       stop('Organism must be provided.',
            'Only "mouse" and "human" are allowed values') }
    if (! organism %in% c('mouse','human')) {
       cat('Only "mouse" and "human" are allowed values for organism')
       stop(organism, 'is not supported.' )
    }

    is_NCBI_id = function(x) { 
       tryCatch({y = as.integer(x) ; TRUE},
        error = function(e) {stop(e)}, 
        warning = function(w) {FALSE})
     }

    is_NCBI_id_flag = is_NCBI_id(genes)
        
    if (organism == 'mouse') {
       mouse_human_map = mouse_genes_to_human_genes(
                                 genes, 
                                 parallel = parallel, 
                                 conf_file = conf_file) %>% 
                           filter(!is.na(human_ID))
       genes = mouse_human_map %>% 
                      pull(human_ID) %>% unique %>% na.omit
       is_NCBI_id_flag = TRUE
    }

    if (!update_info) {
       col_types = cols(
                          gene = col_character(),
                          gene_id = col_double(),
                          disease = col_character(),
                          disease_id = col_character()
                     )
       gene_disease_df = read_csv(
                          system.file(
                             'extdata/disease_ontology.csv',
                              package="ABIgeneRMINC"),
                           col_types = col_types)

       if (!is_NCBI_id_flag) {
            genes = gene_disease_df %>% filter(gene %in% genes) %>% pull(gene_id) %>% unique
       }
       gene_disease_df = gene_disease_df %>% filter(gene_id %in% genes)
       gene_disease_df = gene_disease_df %>% 
                        mutate(gene_id = factor(
                                                 as.character(gene_id),
                                                 levels=unique(as.character(genes)))) %>% 
                        arrange(gene_id) %>% 
                        mutate(gene_id = as.integer(as.character(gene_id)))       
    } else {
       if (is.null(api_key)) {stop('must supply API key to be used by DisGeNET')}

       vocabulary = ifelse(is_NCBI_id_flag,'ENTREZ','HGNC')
   
       gene_to_disease = function(x, vocabulary) {
           res = gene2disease(
              gene = x, 
              vocabulary = vocabulary, 
              database = "ALL", 
              score = c(0.1,1), 
              api_key=api_key)
           if (is.character(res)) {
              if (res == 'no results for the query') { 
                 ret = tibble(gene = NA, gene_id = NA, disease = NA , disease_id = NA)
              } else {
                 stop("unexpected output")
              }
           } else {
                ret = disgenet2r::extract(res) %>% 
                  select(gene_symbol, geneid, disease_name, diseaseid) %>% 
                  rename(
                     gene = 'gene_symbol',
                     gene_id = 'geneid',
                     disease = 'disease_name',
                     disease_id = 'diseaseid'
                     ) %>% 
                  mutate(
                     gene       = as.character(gene),
                     disease    = as.character(disease),
                     disease_id = as.character(disease_id),
                     ) %>% 
                  select(gene,gene_id,disease,disease_id) %>% as_tibble
           }
           return(ret)
        }
   
       parfunc = function(group, vocabulary) {
          bind_rows(lapply(group,gene_to_disease,vocabulary = vocabulary))
       }
   
        if (length(genes) > 1) {
          chunks = split(
                     genes, 
                     ( seq_along(genes) %/% 100) )
        } else {
          chunks  = list(genes)
        }
   
       if (is.null(parallel)) {
           if ("pbapply" %in% rownames(installed.packages())) {
                lapply = pbapply::pblapply
            }
           gene_disease_df = bind_rows(lapply( 
                                   chunks , 
                                   gene_to_disease ,
                                   vocabulary = vocabulary))
           lapply = base::lapply
       } 
       else {
           par_bool_local = parallel[1] %in% c("local", "snowfall")
           par_bool_pb = "pbapply" %in% rownames(installed.packages())
           n_cores <- as.numeric(parallel[2])
           if ( par_bool_local & par_bool_pb ) {
              cl <- makeCluster(n_cores)
              on.exit(stopCluster(cl))
              clusterExport(cl, c(
                    "select", "tibble","gene2disease",
                    'rename','%>%','mutate','as_tibble'))
              gene_disease_df = pbapply::pblapply(
                                       chunks, 
                                       gene_to_disease, 
                                       vocabulary = vocabulary,
                                       cl = cl) %>% 
                              bind_rows
           } 
           else if (par_bool_local) {
              if (length(chunks) > 1) {
                  groups  <- split(
                     chunks, 
                     RMINC:::groupingVector(length(chunks), n_cores))
               } 
              else {
                  groups  <- list(chunks)
               }
              gene_disease_df <- RMINC:::failing_mclapply(
                                     groups, 
                                     parfunc, 
                                     vocabulary = vocabulary, 
                                     mc.cores = n_cores) %>% 
                              bind_rows
           } 
           else {
              if (length(chunks) > 1) {
                  groups  <- split(
                     chunks, 
                     RMINC:::groupingVector(length(chunks), n_cores))
               } 
              else {
                  groups  <- list(chunks)
               }
              regdir_name = RMINC:::new_file("diseaseOntology_registry")
   
              cat('Setting up parallel jobs on cluster with the following configurations\n')
              cat(paste0('    Configuration file: ',conf_file,'\n'))
              cat(paste0('    Registry directory: ',getwd(),'/',regdir_name,'\n'))
              suppressMessages({
                 reg = makeRegistry(
                             regdir_name, 
                             packages = c('tidyverse','disgenet2r'), 
                             conf.file = conf_file
                           )
                 batchExport(reg=reg , list(
                          gene_to_disease = gene_to_disease
                    ))
               })
              on.exit({ 
                 cat('Cleaning up parallel jobs. Deleting jobs registry in 5 seconds...')
                 suppressMessages({
                    RMINC:::tenacious_remove_registry(reg) 
                 })
                 cat('Done\n')
               })
              suppressWarnings({suppressMessages({
                   ids <- batchMap(
                                reg = reg, 
                                fun = parfunc,
                                group = groups,
                                more.args = list(
                                       vocabulary = vocabulary
                                ))
               })})
   #           suppressMessages({
                 submitJobs(ids, reg = reg)
                 waitForJobs(reg = reg)
                 i = 0 ; while (i<4) {
                    ndj = findNotDone(reg=reg)
                    if (nrow(ndj) == 0) break
                    if (nrow(getErrorMessages(reg=reg))!=0) break
                    killJobs(ids=ndj, reg=reg)
                    resetJobs(ids=ndj, reg=reg)
                    submitJobs(ids=ndj, reg=reg)
                    waitForJobs(reg = reg)
                    i = i + 1
                 }
   #            })
              not_done_ids = findNotDone(reg=reg)$job.id
              if (length(not_done_ids)>0) {
                 cat(length(not_done_ids),'jobs have not been completed','\n')
                 cat('Here is the message in job',not_done_ids[1],':\n')
                 cat('------',tail(getLog(not_done_ids[1]),1),'------\n')
                 stop(paste0('Some ',parallel[1],' jobs have failed'))
              }
              suppressMessages({
                 gene_disease_df = reduceResultsList(reg = reg) %>% bind_rows
               })
            }
       }  
    }
    return(gene_disease_df)
}


#' @export
disease_ontology_enrichment_analysis = function(
              genes_target, genes_background,
              organism = 'mouse',
              parallel = NULL,
              update_info = FALSE,
              api_key = NULL,
              conf_file = getOption("RMINC_BATCH_CONF")
    ) {

    if (length(organism) != 1) { 
       stop('Organism must be provided.',
            'Only "mouse" and "human" are allowed values') }
    if (! organism %in% c('mouse','human')) {
       cat('Only "mouse" and "human" are allowed values for organism')
       stop(organism, 'is not supported.' )
    }

    if (update_info) { if (is.null(api_key)) {
       stop('api_key for DisGeNET MUST be supplied if set update_info to TRUE')
    }}

    bool = genes_target %in% genes_background
    if ( !all(bool) ) {
       cat(sum(!bool), 'genes in target not found in background.\n')
       if (sum(!bool) > 5) { 
         cat('Here are the first 5:\n') 
         cat(paste0('    ',genes_target[!bool][1:5]),sep='\n')
        } else {
         cat('Here are the missing genes:\n') 
         cat(paste0('    ',genes_target[!bool]),sep='\n')
        }
       warning('missing genes will be excluded from target')
       genes_target = genes_target[bool]
    }

    is_NCBI_id = function(x) { 
       tryCatch({y = as.integer(x) ; TRUE},
        error = function(e) {stop(e)}, 
        warning = function(w) {FALSE})
     }

    if (!is_NCBI_id(genes_target))   { 
      genes_target     = gene_to_NCBI(
                              genes_target,
                              organism = organism,
                              update_info = update_info
                            )$NCBI
    }
    if (!is_NCBI_id(genes_background))  { 
       genes_background = gene_to_NCBI(
                              genes_background,
                              organism = organism,
                              update_info = update_info
                            )$NCBI
    }

    
    if (organism == 'mouse') {
       mouse_human_map = mouse_genes_to_human_genes(
                                 genes_background, 
                                 parallel = parallel, 
                                 conf_file = conf_file) %>% 
                           filter(!is.na(human_ID))
       genes_target = mouse_human_map %>% 
                           filter(mouse_ID %in% genes_target) %>% 
                           pull(human_ID)
       genes_background = mouse_human_map %>% 
                           pull(human_ID)
    }
    
   
   disease_ontology_df = disease_ontology(
              genes_background,
              organism = 'human',
              update_info = update_info,
              api_key = api_key,
              parallel = parallel,
              conf_file = conf_file
    ) 
   
   disease_ontology_df = disease_ontology_df %>% mutate( drawn = gene_id %in% genes_target )

   disease_gene_list = disease_ontology_df %>% filter(drawn) %>% 
        (function(df) split(df$gene,df$disease))

   # determine global hyper-geometric test parameters 
   # ( m_plus_n = total # of genes )
   # ( k = # of drawn genes )
    tmp_hypergeo_df = disease_ontology_df %>% 
        select(gene,drawn) %>% unique
    m_plus_n = nrow(tmp_hypergeo_df)
    k = sum(tmp_hypergeo_df$drawn)

    # determine hyper-geometric test parameters for each disease term
    # (q = number of drawn genes that are white)
    # (m = number of genes that are white)
    q_list = disease_gene_list %>% 
        lapply(function(x) length(x))
    m_list = disease_ontology_df %>%
        (function(df) split(df$gene,df$disease)) %>% 
        lapply(function(x) length(x)) %>% 
        .[names(q_list)]

    hyper_res = Map(list, q_list, m_list, m_plus_n, k, disease_gene_list) %>% 
      lapply(function(x) {
            q = x[[1]]
            m = x[[2]]
            n = x[[3]] - x[[2]]
            k = x[[4]]
            ret = tibble( 
                     p = phyper(q - 1,m,n,k,lower.tail=FALSE),
                     `Enrichment (N, B, n, b)` = paste0(
                                                   round((q/k)/(m/(m+n)),2),
                                                  ' (',m+n,',',m,',',k,',',q,')'),
                     Genes = x[5]
                    )
      }) %>% (function(x) {
          nx = names(x)
          tibble(disease = nx) %>% bind_cols(bind_rows(x))
      })
   hyper_res$q = p.adjust(hyper_res$p,method='fdr') 

   retdf = hyper_res %>% 
         inner_join(
            disease_ontology_df %>% select(disease,disease_id) %>% unique,
            by = 'disease'
           ) %>% 
         select(disease_id, disease, p, q, `Enrichment (N, B, n, b)`, Genes) %>% 
         arrange(q,p)

   return(retdf)
}








