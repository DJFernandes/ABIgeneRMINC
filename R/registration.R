#' Register two minc volumes using ANTs
#'
#' Uses system ANTs installation to register minc volumes
#' @param source_volume volume to move. Can be an RMINC object or filename.
#' @param target_volume reference volume. Can be an RMINC object or filename.
#' @param source_mask_volume volume for masking source_volume
#' @param target_mask_volume volume for masking target_volume
#' @param xfm_file location to save file. If NULL, file is saved in temporary directory.
#' @param verbose If true, ANTs messages are passed to console
#' @return location of xfm transform file
#' @importFrom RMINC mincWriteVolume
#' @import processx pbapply
#' @export
minc_ants_register = function( 
          source_volume , 
          target_volume , 
          source_mask_volume , 
          target_mask_volume , 
          xfm_file = NULL , 
          verbose = T ) {
   if (is.character(source_volume)) {
      RMINC:::mincFileCheck(source_volume)
   }
   if (is.character(target_volume)) {
      RMINC:::mincFileCheck(target_volume)
   }
   if (is.character(source_mask_volume)) {
      RMINC:::mincFileCheck(source_mask_volume)
   }
   if (is.character(target_mask_volume)) {
      RMINC:::mincFileCheck(target_mask_volume)
   }

   if (is.null(xfm_file)) {
      xfm_file = tempfile(pattern='transform_',tmpdir = tempdir(),fileext='.xfm')
      outdir = tempdir()
   } else {
      outdir = dirname(xfm_file)
      if ( ! outdir %in% c('.','') ) {
         dir.create(outdir, showWarnings = FALSE)
      }
   }
   
   tmpfl_source = tempfile(
                              pattern='source_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_target = tempfile(
                              pattern='target_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_target_padded = tempfile(
                              pattern='target_padded_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_target_mask = tempfile(
                              pattern='target_mask_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_target_mask_padded = tempfile(
                              pattern='target_mask_padded_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_source_mask = tempfile(
                              pattern='source_mask_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_source_affine = tempfile(
                              pattern='source_affine_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_source_affine_mask = tempfile(
                              pattern='source_affine_mask_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_xfm_affine = tempfile(
                              pattern='xfm_affine_',
                              tmpdir=outdir,
                              fileext='.xfm')
   tmpfl_xfm_affine_inv = gsub('.xfm$','_inverse.xfm',tmpfl_xfm_affine)
   tmpfl_xfm_warp = tempfile(
                              pattern='xfm_warp_',
                              tmpdir=outdir,
                              fileext='.xfm')
   tmpfl_xfm_warp_inv = gsub('.xfm$','_inverse.xfm',tmpfl_xfm_warp)
   ants_log_o = tempfile(
                              pattern='ants_',
                              tmpdir=outdir,
                              fileext='o.log')

   ants_log_e = tempfile(
                              pattern='ants_',
                              tmpdir=outdir,
                              fileext='e.log')

   rm_tmp_files = function() {   
      unlink(tmpfl_source)
      unlink(tmpfl_target)
      unlink(tmpfl_target_padded)
      unlink(tmpfl_target_mask)
      unlink(tmpfl_target_mask_padded)
      unlink(tmpfl_source_mask)
      unlink(tmpfl_source_affine)
      unlink(tmpfl_source_affine_mask)
      unlink(tmpfl_xfm_affine)
      unlink(tmpfl_xfm_affine_inv)
      unlink(tmpfl_xfm_warp)
      unlink(tmpfl_xfm_warp_inv)
      y = Sys.glob(gsub('.xfm$','_grid_*.mnc',tmpfl_xfm_warp))     ; if (length(y)>0) {unlink(y)}
      y = Sys.glob(gsub('.xfm$','_grid_*.mnc',tmpfl_xfm_warp_inv)) ; if (length(y)>0) {unlink(y)}
      unlink(ants_log_o)
      unlink(ants_log_e)
    }
    
   rm_tmp_files()
   on.exit( rm_tmp_files() )



   if (!is.character(source_volume)) {
      mincWriteVolume(source_volume, output.filename = tmpfl_source)
      source_volume = tmpfl_source
   }
   if (!is.character(target_volume)) {
      mincWriteVolume(target_volume, output.filename = tmpfl_target)
      target_volume = tmpfl_target
   }
   if (!is.character(source_mask_volume)) {
      mincWriteVolume(source_mask_volume, output.filename = tmpfl_source_mask)
      source_mask_volume = tmpfl_source_mask
   }
   if (!is.character(target_mask_volume)) {
      mincWriteVolume(target_mask_volume, output.filename = tmpfl_target_mask)
      target_mask_volume = tmpfl_target_mask
   }

  
   ants_affine_cmd = paste0(
          'antsRegistration -d 3 ',
          '-o ',gsub('.xfm$','',tmpfl_xfm_affine),' ',
          '-x [',source_mask_volume,',',target_mask_volume,'] ',
          '-a 1 -z 1 ',
          '-t Affine[0.25] ',
          '-m MI[',source_volume,',',target_volume,',1,32,Regular,0.25]',' ',
          '--convergence [1000,1e-6,10] ', 
           '--shrink-factors 1 ',
          '--smoothing-sigmas 0mm ',
          '--minc -v 1 --use-histogram-matching 1'
   )


   mincvol_pad_crop(
        infile = target_volume, 
        outfile = tmpfl_target_padded, 
        operation = 'pad')

   mincvol_pad_crop(
        infile = target_mask_volume, 
        outfile = tmpfl_target_mask_padded, 
        operation = 'pad')

   affine_resample_cmd = paste0(
          'mincresample ',
          '-quiet -like ',tmpfl_target_padded,' ',
          '-transform ',tmpfl_xfm_affine,' ',
          source_volume,' ',
          tmpfl_source_affine
   )

   affine_resample_mask_cmd = paste0(
          'mincresample ',
          '-quiet -nearest -labels -like ',tmpfl_target_mask_padded,' ',
          '-transform ',tmpfl_xfm_affine,' ',
          source_mask_volume,' ',
          tmpfl_source_affine_mask
   )

   ants_warp_cmd = paste0(
          'antsRegistration -d 3 ',
          '-o ',gsub('.xfm$','',tmpfl_xfm_warp),' ',
          '-x [',tmpfl_source_affine_mask,',',target_mask_volume,'] ',
          '-a 1 -z 1',' ', 
          '-t SyN[0.1,2,0] ',
          '-m MI[',tmpfl_source_affine,',',target_volume,',1,32,Regular,0.25] ',
          '--convergence [1000x1000x1000x1000,1e-6,10] ',
          '--shrink-factors 8x4x2x1 ',
          '--smoothing-sigmas 0.2x0.1x0.05x0mm ',
          '--minc -v 1 --use-histogram-matching 1'
   )

   concat_xfm_cmd =  paste0(
          'xfmconcat ',
          tmpfl_xfm_affine, ' ',
          tmpfl_xfm_warp, ' ',
          xfm_file
   )

   if (verbose) {
      x = system(ants_affine_cmd,intern=F)
      x = system(affine_resample_cmd,intern=F)
      x = system(affine_resample_mask_cmd,intern=F)
      x = system(ants_warp_cmd,intern=F)
      x = system(concat_xfm_cmd,intern=F)
   } else {
      stdout = ants_log_o ; stderr = ants_log_e
      extract_prog_from_log = function(stdout,epoch_maxes) {
         rl = readLines(stdout)
         
         bool = grepl('^  Elapsed time ',rl)
         if (any(bool)) { 
              val = 1 
          } else {
              grepstr = 'DIAGNOSTIC,Iteration,metricValue,convergenceValue,ITERATION_TIME_INDEX,SINCE_LAST'
              grepres = grep(grepstr,rl)
              if (length(grepres) != 0) {
                 str_idx =  min(grepres)
                 rl = rl[(str_idx+1):length(rl)]
                 rl = rl[grep(',',rl)]
   
                 val = -log(abs(as.numeric(strsplit(tail(rl,1),split=',')[[1]][4])),base=10)
                 if (val<0) {val = 0}
                 if (val>6) {val = 6}
                 val = val/6  # 6 is the -log10 of the convergence criteria
                 
                 if (length(epoch_maxes) > 1) {
                    i_epoch = length(grepres)
                    if (i_epoch == 1) {
                       val = val * epoch_maxes[i_epoch]
                    } else {
                       pre_epoch_val = epoch_maxes[i_epoch-1]
                       delta_epoch_val = epoch_maxes[i_epoch] - pre_epoch_val
                       val = val * ( delta_epoch_val ) + pre_epoch_val
                    }
                 }
               } else {
                  val = 0
               }
          }
          return(val)
       }
   
      system_cmd = function(
                    command,stdout,stderr,
                    proc_name,term_file,
                    startval,maxval,pb,
                    epoch_weights = 1,
                    downsample_and_retry = F,
                    ...) {
         unlink(stdout) ; unlink(stderr)
         xp = command %>% strsplit(split = ' ') %>% .[[1]]
         
         p = process$new(xp[1], xp[-1] ,stdout = stdout, stderr = stderr)
   
         while (TRUE) {
            if (file.exists(stderr)) {
               if (length(readLines(stderr)) != 0) {
                    stop('Error in ',proc_name)
               }
            }
            if (file.exists(term_file)) {
               val = maxval
               pbapply::setpb(pb,val)
               break
            } else {
               if (file.exists(stdout)) {
                  epoch_maxes = cumsum(epoch_weights)/sum(epoch_weights)
                  dval = (maxval - startval)
                  ival = extract_prog_from_log(stdout,epoch_maxes)
                  val = ival * dval + startval
                  pbapply::setpb(pb,val)
               }
            }
            if (!p$is_alive()) {
               message('\n',proc_name, ' died :(')
               if (downsample_and_retry) {
                  message(paste(
                     "The usual suspect is that the",
                     "files are too big for ANTs to register"))
                  n_iter = length(grep(paste0(
                          'XXDIAGNOSTIC,Iteration,metricValue,',
                          'convergenceValue,ITERATION_TIME_INDEX',
                          ',SINCE_LAST'),
                     readLines(stdout)))
                  message(paste(
                     "Looks like ANTs made it to",
                     "Level" , n_iter, 'before dying'))
                  message(paste(
                     "I will re-run the registration and stop at",
                     "Level" , n_iter-1))
                  
                  ln_iter = grep('--convergence',xp[-1])+1
                  iter_string = xp[-1][ln_iter] %>% 
                     gsub('^.|.$', '', .) %>% 
                     strsplit(split=',') %>% .[[1]]
                  iter_string_new = iter_string[1] %>% 
                     strsplit(split='x') %>% .[[1]] %>% 
                     ( function(x) { 
                        x[ (n_iter) : length(x) ] = 0
                        x
                     }) %>% 
                     paste(collapse='x') %>% 
                     c(iter_string[-1])
                  xp[-1][ln_iter] = paste0('[',paste(iter_string_new,collapse=','),']')

                  command = paste(xp,collapse=' ')
                  system_cmd(
                                command,stdout,stderr,
                                proc_name,term_file,
                                startval,maxval,pb,
                                epoch_weights = epoch_weights,
                                downsample_and_retry = F)
               } else {
                  stop(proc_name, ' died')
               }
            }
          Sys.sleep(5)
          }
      }


      pb = pbapply::startpb(0, 1)
      system_cmd(
            command = ants_affine_cmd,
            stdout = stdout,
            stderr = stderr,
            proc_name = 'affine registration',
            term_file = tmpfl_xfm_affine,
            startval = 0,
            maxval = 0.01,
            pb = pb,
            epoch_weights = 1
         )
      pbapply::setpb(pb,0.001)
   
      x = system(affine_resample_cmd,intern=T)
      pbapply::setpb(pb,0.0011)
      
      x = system(affine_resample_mask_cmd,intern=T)
      pbapply::setpb(pb,0.0012)
      
      system_cmd(
            command = ants_warp_cmd,
            stdout = stdout,
            stderr = stderr,
            proc_name = 'non-affine registration',
            term_file = tmpfl_xfm_warp,
            startval = 0.0012,
            maxval = 0.999,
            pb = pb,
            epoch_weights = c(7.5e-09, 3.8e-06, 0.001953125, 1),
            downsample_and_retry = T
         )
      pbapply::setpb(pb,0.999)
   
      x = system(concat_xfm_cmd,intern=T)
      pbapply::setpb(pb,1)
      pbapply::closepb(pb)
   }   


   return(xfm_file)

}


#' Resample and check registration quality
#'
#' Creates a PNG for manual assessment of registration quality.
#' @param source_file MINC file that is moved
#' @param target_file MINC file that is the target
#' @param source_transform XFM file mapping the source_file to the target_file
#' @param resampled_source_file MINC filename where the transformed sourcefile will be written. Temporary file is created if this is not supplied. 
#' @param reg_pik PNG filename where the image will be written. Temporary file is created if this is not supplied. 
#' @return PNG filename where the image will be written
#' @import RMINC MRIcrotome grid gridExtra
check_registration_quality_png = function(
                     source_file, 
                     target_file, 
                     source_transform, 
                     resampled_source_file = NULL, 
                     reg_pik = NULL) {
         if (is.null(resampled_source_file)) {
            resampled_source_file = tempfile(fileext='.mnc')
         }
         if (is.null(reg_pik)) {
            reg_pik = tempfile(fileext='.png')
         }
         resample_cmd = paste0(
               'mincresample -quiet ',
                  ' -like ',target_file,
                  ' -transform ', source_transform,' ',
                  source_file, ' ', resampled_source_file)
          
          system(resample_cmd)
          
          x = mincArray(mincGetVolume(resampled_source_file))
          y = mincArray(mincGetVolume(target_file))
          slcs = round(dim(y)/2)

          vxl_bool1 = array(F,dim(y))
             vxl_bool1[slcs[1],,] = T
             vxl_bool1[,slcs[2],] = T
             vxl_bool1[,,slcs[3]] = T
          vxl_bool2 = array(F,dim(y))
             vxl_bool2[y>8] = T
          vxl_bool = vxl_bool1 & vxl_bool2
          
          col_bounds = quantile(x[vxl_bool],c(0.01,0.99))
          countour_lvls = quantile(y[y>8],c(0.25,0.5))

          mrslc_grobs = lapply(1:3, function(i) {
             sliceSeries(nrow=1,slice=slcs[i],dimension=i) %>% 
                anatomy(y, low=0, high=max(y)) %>% 
                contours(y,levels=countour_lvls,lty=c(1,2,3)) %>% 
             sliceSeries(nrow=1,slice=slcs[i],dimension=i) %>% 
                anatomy(x, low=col_bounds[1], high=col_bounds[2]) %>% 
                contours(y,levels=countour_lvls,lty=c(1,2,3)) %>% 
                grobify()
           }) %>% 
           arrangeGrob(grobs=.,ncol=1)

          gp=gpar(col='white')
          row_ttl_grobs = c('Sagittal','Coronal','Transverse') %>% 
                lapply(function(txt) {
                     textGrob(txt,gp=gp,rot=90)
                }) %>% arrangeGrob(grobs=.,ncol=1)
          colm_ttl_grobs = c('ABI','resampled template') %>% 
                lapply(function(txt) {
                     textGrob(txt,gp=gp)
                }) %>% arrangeGrob(grobs=.,nrow=1)
          
          check_reg_grob = arrangeGrob(
                               mrslc_grobs, row_ttl_grobs, colm_ttl_grobs,
                               layout_matrix = matrix(c(
                                                   NA , 3 ,
                                                   2  , 1
                                               ),ncol=2, byrow=T),
                               widths=c(0.05,0.9), heights=c(0.05,0.9)
                            ) %>% 
                            grobTree(rectGrob(gp=gpar(col='black',fill='black')),.)
          png(file=reg_pik,width=7,height=7,res=150,units='in')
             grid.newpage()
             grid.draw(check_reg_grob)
          dev.off()
          cat('Check if registration was successful: ')
          cat(reg_pik,'\n')
          return(reg_pik)
 }

#' Download and register to ABI 50um template
#'
#' Downloads ABI 50um template and registers minc data to it
#' @param anatomy_template minc data to move
#' @param anatomy_mask mask used in the registration
#' @param allen_ccf3v_template Path to the ABI template. If NULL, ABI 50um template is downloaded and stored to a temporary file. If path does not point to an existing file, the template is downloaded and written here. If path points to an existing file, the file is used as the target for registration
#' @param xfm_files_dir Directory to store the XFM transformation files. If NULL, a temporary directory is used. 
#' @param keep_files If TRUE, pathnames to ABI template and XFM directories will be echoed. Its User's responsibility make sure these files are in a safe location. 
#' @return XFM transformation filepath
#' @importFrom nat read.nrrd
#' @import RMINC
#' @export
ABI_template_align = function(
      anatomy_template, 
      anatomy_mask, 
      allen_ccf3v_template = NULL, 
      allen_ccf3v_mask     = NULL, 
      xfm_files_dir = tempdir(),
      keep_files = F) {

         xfm_files_dir = ifelse(
                           grepl('/$',xfm_files_dir),
                           xfm_files_dir,
                           paste0(xfm_files_dir,'/'))

         if (is.character(anatomy_template)) {
            RMINC:::mincFileCheck(anatomy_template)
         } else {
            tmpfl_source = tempfile(pattern='anatomy_',fileext='.mnc')
            mincWriteVolume(anatomy_template, output.filename = tmpfl_source)
            anatomy_template = tmpfl_source
         }

         if (is.character(anatomy_mask)) {
            RMINC:::mincFileCheck(anatomy_mask)
         } else {
            tmpfl_source = tempfile(pattern='anatomy_',fileext='.mnc')
            mincWriteVolume(anatomy_mask, output.filename = tmpfl_source)
            anatomy_mask = tmpfl_source
         }
   

         if (length(allen_ccf3v_template) > 1) {stop('Invalid argument for allen_ccf3v_template. Must be either null or single-element.')}
         if (is.null(allen_ccf3v_template) | is.null(allen_ccf3v_mask)) {
            dld_flag = T
            allen_ccf3v_template = tempfile(pattern='allen_CCFV3',fileext='.mnc')
            allen_ccf3v_mask     = tempfile(pattern='allen_mask_CCFV3',fileext='.mnc')
          } else if (
                        (file.access(as.character(allen_ccf3v_template)) == -1) |
                        (file.access(as.character(allen_ccf3v_mask)) == -1) 
           ) {
            dld_flag = T
          } else {
            dld_flag = F
          }
          
          if (dld_flag) {
            cat('Downloading Allen CCFV3 template...')
            x = allen_s2p_template_download(allen_ccf3v_template)
            x = allen_s2p_template_download(allen_ccf3v_mask, mask=T) 
            cat(' Done\n')
            if (keep_files) {
               cat('   The Allen seriel two-photon template filepath is:')
               cat(allen_ccf3v_template,'\n')
               cat('   The Allen seriel two-photon template mask filepath is:')
               cat(allen_ccf3v_mask,'\n')
            }
          }

         cat('Registering anatomy to ABI...\n')
         cat(' ... This will probably take a while ... ')
         anatomy_transform = paste0(xfm_files_dir,'anatomy_to_abi.xfm')
         x = minc_ants_register(
                   anatomy_template, 
                   allen_ccf3v_template, 
                   anatomy_mask, 
                   allen_ccf3v_mask, 
                   xfm_file = anatomy_transform, 
                   verbose=F)
         cat('Finished registering anatomy to ABI\n')
         if (keep_files) {
            cat('   The transformation from anatomy to ABI is in this directory:\n')
            cat('       ',xfm_files_dir,'\n')
         }
         resampled_anatomy_template = paste0( xfm_files_dir, 'anatomy_resampled.mnc' )
         reg_pik = paste0( xfm_files_dir, 'check_registration_quality.png' )
         x = check_registration_quality_png(
                       anatomy_template, 
                       allen_ccf3v_template, 
                       anatomy_transform,
                       resampled_anatomy_template,
                       reg_pik)
         
         if ( ( dld_flag & keep_files ) | (!dld_flag) ) {
            attributes(anatomy_transform)$allen_ccf3v_template = allen_ccf3v_template
            attributes(anatomy_transform)$allen_ccf3v_mask = allen_ccf3v_mask
         } 
         
         return(anatomy_transform)
}

#' Gene expression spatial enrichment
#' 
#' Conducts a gene expression spatial enrichment analysis for arbitrary anatomy statistics. Statistics above the target threshold constitute the target region-of-interest (ROI). Statistics below a contrast threshold constitute the contrast region. Gene expression spatial enrichment is calculated for any number of genes in the ABI gene expression atlas. Enrichments is computed using a fold-change measure: expression in target ROI divided by expression in contrast region. If contrast threshold is not suppled, then the contrast regions is assumed to be the whole brain. 
#' @param anatomy_statistics MINC file or vector denoting statistics at each voxel
#' @param gene_expression_paths Filenames or URLs pointing to gene expression data. If NULL, then genome-wide gene expression analysis is conducted and gene expression data is downloaded to temporary files. 
#' @param target_threshold statistics greater than this value constitute the ROI
#' @param contrast_threshold statistics less than this value constitute the contract region. If NULL, then contrast region is assumed to be the whole brain. 
#' @param symmetric_statistics Should the absolute value of statistics at each voxel be considered instead of the signed value?
#' @param ABI_registration arguments that can be supplied if you want to perform a registration prior to analysis. If you provide an XFM filepath as 'anatomy_transform', then the statistics will be transformed and resampled prior to gene expression analysis. If you provide 'anatomy_template' (filepath or MINC vector denoting the anatomy where statistics were conducted), 'anatomy_mask' (filepath or MINC vector denoting brain mask for anatomy), then the anatomy is registered to the ABI template and the resulting transformation applied to the statistics prior to gene expression computation. In this case, ABI template and mask are downloaded unless 'allen_ccf3v_template' and 'allen_ccf3v_mask' point to filepaths containing the respective information. 
#' @param gene_expression_analyis_options options for gene expression analysis. 'interpolate_gene' flags whether nearest-neighbour imputation should be done to fill in missing gene expression voxels. If FALSE, then these voxels are ignored. 'reflect_gene' flags whether gene expression signal should be reflected across the sagittal midplane to fill in missing gene expression voxels in the opposite hemisphere. 'brain_mask' is a MINC file, RAW file, or vector identifying which voxels are in the brain. If NULL, the mask will be downloaded from the ABI. 
#' @param tmp_files_dir location to store temporary files from the registration
#' @param parallel how many processors to run on (default=single processor). Specified as a two element vector, with the first element corresponding to the type of parallelization, and the second to the number of processors to use. For local running set the first element to "local" or "snowfall" for back-compatibility, anything else will be run with batchtools. Leaving this argument NULL runs sequentially.
#' @param conf_file A batchtools configuration file defaulting to \code{getOption("RMINC_BATCH_CONF")}
#' @return gene expression enrichment expressed as fold-change 
#' @importFrom RCurl url.exists
#' @import RMINC batchtools
#' @export
adult_gene_expression_analysis = function(
       anatomy_statistics, 
       gene_expression_paths = NULL,
       target_threshold, 
       contrast_threshold = NULL,
       symmetric_statistics = F, 
       ABI_registration = list( 
            anatomy_template = NULL, 
            anatomy_mask = NULL, 
            allen_ccf3v_template = NULL,
            allen_ccf3v_mask = NULL,
            anatomy_transform = NULL
         ),
       gene_expression_analyis_options = list(
            interpolate_gene = F , 
            reflect_gene     = T ,
            brain_mask = NULL
        ),
       tmp_files_dir = NULL,
       parallel = NULL,
       conf_file = getOption("RMINC_BATCH_CONF")
 ) {

   # make tmp directory if it doesn't exist
   if (is.null(tmp_files_dir)) {
       tmp_files_dir = paste0( tempdir() , '/')
       unlink(tmp_files_dir,recursive = T)
       dir.create(tmp_files_dir)
       keep_files = F
    } else {
       tmp_files_dir = paste0( tmp_files_dir , '/')
       if (!dir.exists(tmp_files_dir)) { stop(tmp_files_dir, " already exists") } 
       keep_files = T
    }


   ABI_reg_flag = !all(sapply(ABI_registration, is.null))
   if (ABI_reg_flag) {
      cat('transforming statistics to ABI space\n')
      anatomy_transform    = ABI_registration[['anatomy_transform']]
      anatomy_mask         = ABI_registration[['anatomy_mask']]
      anatomy_template     = ABI_registration[['anatomy_template']]
      allen_ccf3v_template = ABI_registration[['allen_ccf3v_template']]
      allen_ccf3v_mask     = ABI_registration[['allen_ccf3v_mask']]
      
      ants_registration_flag = !( is.null(anatomy_template) | is.null(anatomy_mask) )
      resample_flag = !is.null(anatomy_transform)

      if ( ! ( ants_registration_flag | resample_flag ) ) {
         cat('If ABI_registration list has components that are not NULL\n')
         cat('   one of the following conditions must be met:\n')
         cat('    1) both anatomy_template and anatomy_mask must be supplied\n')
         cat('    2) anatomy_transform must be supplied\n')
         stop('Cannot perform registration')
      }

      if ( ants_registration_flag & resample_flag ) {
         cat('You provided both anatomy_template (with mask) and anatomy_transform. \n')
         cat("Since you gave me the transform, I am assuming you don't actually want me\n")
         cat("to use ANTS and find a transformation (which is a time-consuming process)\n")
         cat("Instead, I think you just want to resample using the transform provided.\n")
         warning('Ignoring anatomy_template and anatomy_mask')
         cat("Using anatomy_transform to resample.\n")
         ants_registration_flag = F
      }

      if (ants_registration_flag) {
         xfm_files_dir = paste0(tmp_files_dir,'anatomy_to_abi_transform/')
         if (!dir.exists(xfm_files_dir)) {
            dir.create(xfm_files_dir,recursive = T)
         } else {
            stop(xfm_files_dir, " already exists")
         }

         anatomy_transform = ABI_template_align(
              anatomy_template, 
              anatomy_mask, 
              allen_ccf3v_template, 
              allen_ccf3v_mask, 
              xfm_files_dir,
              keep_files)
         resample_flag = T
      }

      if (resample_flag) {
          anatomy_statistics = ABI_resample(
                source_volume = anatomy_statistics, 
                xfm_transform = anatomy_transform, 
                target_space_name = 'adult gene expression')
       }
   }

   if (is.character(anatomy_statistics)) {
      RMINC:::mincFileCheck(anatomy_statistics)
      anatomy_statistics = mincVectorTOallenVector(mincGetVolume(anatomy_statistics))
   }
   
   if (symmetric_statistics) {
      anatomy_statistics = abs(anatomy_statistics)
   }

   # process brain mask
   if ( is.null(gene_expression_analyis_options$brain_mask) ) {
      cat('Downloading ABI mask for gene expression...')
      brainmask = allen_grid_labels_download()
      cat('Done!\n')
    } else if (length(gene_interpolate_options$brain_mask) == 1) {
      if (grep(".mnc$",gene_interpolate_options$brain_mask)) { 
           brainmask = mincGetVolume(gene_interpolate_options$brain_mask) %>% 
               mincVectorTOallenVector
       } else if (grep(".raw$",gene_interpolate_options$brain_mask)) {
           brainmask = read.raw.gene(gene_interpolate_options$brain_mask,labels=TRUE,url=T)
       } else {
          stop(
            'Only *.mnc and *.raw can be valid gene expression mask file:',
            gene_interpolate_options$brain_mask)
       }
    }
   brainmask_attrs = attributes(brainmask)
   brainmask = brainmask > 0.5
   attributes(brainmask) = brainmask_attrs

   # if gene expression files are not provided, they will be queried and downloaded from ABI
   if (is.null(gene_expression_paths)) {
      download_gene_flag = T
      cat('Downloading and analysing all ABI genes...\n')
#      genes = readLines(system.file('extdata/ABI_genes.csv',package="ABIgeneRMINC")) #TODO
      genedf = read.csv('inst/extdata/ABI_genes.csv',stringsAsFactors=F)
      gene_expression_paths = genedf$URL
      isURL_vec=rep(TRUE,length(gene_expression_paths))
    } else {
      cat('Analysing ABI genes...\n')
      download_gene_flag = F
      isURL_vec = rep(FALSE,length(gene_expression_paths))
      isURL_vec[ ! file.exists(gene_expression_paths) ] = TRUE
      bool = !sapply( gene_expression_paths[ isURL_vec ] , url.exists )
      if ( any( bool ) ) {
         cat('The following path(s) cannot neither be found on disk nor are valid URLs:\n')
         cat( gene_expression_paths[ bool ] , sep = '\n' )
         cat('\n')
         stop('Invalid gene expression paths found')
      }
    }

   catcmd = function(download_gene_flag) {
      if (download_gene_flag) { cat('Downloading and analysing all ABI genes...')
      } else { cat('Analysing ABI genes...') }
   }

   if ( gene_expression_analyis_options$interpolate_gene ) {
   if ( !gene_expression_analyis_options$reflect_gene ) {
      cat('It does not really make sense to interpolate gene expression data\n')
      cat('  but not reflect gene expression data across the sagittal midplane.\n')
      cat('  So, analysis will proceed assuming gene expression data should be reflected\n')
      warning('Setting reflect_gene to TRUE')
      gene_expression_analyis_options$reflect_gene = TRUE
   }}

   compute_fold_change = function(
               anatomy_statistics, 
               gene_expression_file, 
               brainmask, 
               tgt.thresh, 
               cntrst.thresh,
               isURL,
               interpolate_gene_flag,
               reflect_gene_flag          
             ) {
          genedata=read.raw.gene(gene_expression_file,url=isURL)
          if (reflect_gene_flag) {
             genedata = midplane.reflect(genedata)
           }
          if (interpolate_gene_flag) {
             genedata = interpolate.gene(genedata, brainmask)
           }
          geneFoldChange(
                stats = anatomy_statistics,
                gene = genedata,
                maskvector = brainmask,
                tgt.thresh = tgt.thresh,
                cntrst.thresh = cntrst.thresh
             )
       }

    parfunc = function(group) {
      map2( group[[1]] , group[[2]] , list) %>% 
       sapply(function(x) {        
                      compute_fold_change(
                         anatomy_statistics,
                         x[[1]],
                         brainmask,
                         target_threshold,
                         contrast_threshold,
                         x[[2]],
                         gene_expression_analyis_options$interpolate_gene,
                         gene_expression_analyis_options$reflect_gene
                      )
       })}

    if (is.null(parallel)) {
        browser()
        if ("pbapply" %in% rownames(installed.packages())) {
             lapply = pbapply::pblapply
         }
        groups <- map2(list(gene_expression_paths), list(isURL_vec), list)
        catcmd(download_gene_flag)
        out = unlist(lapply( groups , parfunc ))
        cat('Done\n')
        lapply = base::lapply
    } else {
        n_groups <- as.numeric(parallel[2])
        if (length(gene_expression_paths) > 1) {
            gene_groups  <- split(
               gene_expression_paths, 
               RMINC:::groupingVector(length(gene_expression_paths), n_groups))
            isURL_groups <- split(
               isURL_vec, 
               RMINC:::groupingVector(length(isURL_vec), n_groups))
        }
        else {
            gene_groups  <- list(gene_expression_paths)
            isURL_groups <- list(isURL_vec)
        }
        groups <- map2(gene_groups, isURL_groups, list)
        if (parallel[1] %in% c("local", "snowfall")) {
            catcmd(download_gene_flag)
            if ("pbapply" %in% rownames(installed.packages())) {
                out = pbapply::pblapply(cl=n_groups, 
                       1:length(gene_expression_paths),
                       function(i) {
                           compute_fold_change(
                              anatomy_statistics,
                              gene_expression_paths[[i]],
                              brainmask,
                              target_threshold,
                              contrast_threshold,
                              isURL_vec[[i]],
                              gene_expression_analyis_options$interpolate_gene,
                              gene_expression_analyis_options$reflect_gene
                           )
                 }) %>% unname %>% unlist
             } else {
                 out = RMINC:::failing_mclapply(groups, function(group) {
                      parfunc(group)
                  }, mc.cores = n_groups) %>% unname %>% unlist
             }
             cat('Done\n')
        }
        else {
            regdir_name = RMINC:::new_file("gene_registry")

            cat('Setting up parallel jobs on cluster with the following configurations\n')
            cat(paste0('    Configuration file: ',conf_file,'\n'))
            cat(paste0('    Registry directory: ',getwd(),'/',regdir_name,'\n'))
            suppressMessages({
               reg = makeRegistry(
                           regdir_name, 
                           packages = c('tidyverse','ABIgeneRMINC'), 
                           conf.file = conf_file
                         )
               batchExport(reg=reg , list(
                        anatomy_statistics = anatomy_statistics,
                        brainmask = brainmask,
                        target_threshold = target_threshold,
                        contrast_threshold = contrast_threshold,
                        gene_expression_analyis_options = gene_expression_analyis_options
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
                              function(group) parfunc(group), 
                              group = groups)
             })})
            
            catcmd(download_gene_flag)
            suppressMessages({
               submitJobs(ids, reg = reg)
               waitForJobs(reg = reg)
               out <- reduceResultsList(reg = reg) %>% unname %>% unlist
             })
            cat('Done\n')
        }
    }

    if (download_gene_flag) {
       ret = genedf ; ret$fold_change = out
    } else {
       ret = out
    }
    return(ret)
}



