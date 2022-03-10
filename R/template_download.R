#' Download ABI 50um template
#'
#' Downloads and reads the ABI 50um template with an option to write to a mincfile
#' @param outfile Optional outfile to write
#' @param mask If true, then the mask is downloaded. Default is FALSE, so the template is downloaded. 
#' @param labels If true, then the labels are downloaded. Default is FALSE, so the template is downloaded. Warning: labels may have integers that can't be written to a MINC file.
#' @return template data as a 1D vector of class mincSingleDim
#' @importFrom nat read.nrrd
#' @import RMINC
allen_s2p_template_download = function(outfile=NULL, mask=F, labels = F) { 
   read.nrrd = nat::read.nrrd
   if (is.null(outfile)) {
      save_file=F
      outfile = tempfile(fileext='.mnc')
    } else {
      save_file=T
      if (file.exists(outfile)) stop(outfile,' is an existing file. Delete it before running.')
    }
   if (!mask & !labels) {
      url_link = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_50.nrrd'
    } else {
      url_link = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_50.nrrd'
    }
   tmpfl_nrrd = tempfile()
   download.file(url_link,tmpfl_nrrd,quiet=T)

   if (mask || labels) {
      # innocuous warning when reading allen labels 
      h = function(w) { 
           if ( w$message == "'signed = FALSE' is only valid for integers of sizes 1 and 2" ) {
              invokeRestart( "muffleWarning" )
           }}        
      vol = withCallingHandlers( read.nrrd(tmpfl_nrrd), warning = h )
    } else {
      vol = read.nrrd(tmpfl_nrrd)
    }

   attr_vol = attributes(vol)
   if (mask) { 
     vol = as.integer(vol>0.5) ; dim(vol) = attr_vol$dim 
   }
   imagesize = attr_vol$header$sizes
   spacing = diag(attr_vol$header$`space directions`)/1000
   origin_coord = c(-6.31 , -4.325 , -5.44)
   direction = diag(c(1,1,1))

   allen_data = vol ; attributes(allen_data)$sizes = rev(dim(vol))
   attr(allen_data,'ANTsRparam') = list(
          imagesize = imagesize,
          spacing = spacing,
          origin = origin_coord,
          direction = direction
    )
   allen_data = allenVectorTOmincVector(allen_data)

   template_file = system.file('extdata/Dorr_resampled_200um.mnc',package="ABIgeneRMINC")

   tmpfl1 = tempfile(fileext='.mnc')
   
   if (mask || labels) {
    tmpfl2 = tempfile(fileext='.mnc')
   }
   
   cmd1 = paste0(
            'mincresample -quiet ', 
            ifelse(mask,'-labels',''),' ',
            ' -nelements ',paste(attr(allen_data,'ANTsRparam')$imagesize,collapse=' '),
            ' -step ',paste(attr(allen_data,'ANTsRparam')$spacing,collapse=' '),
            ' -start ',paste(attr(allen_data,'ANTsRparam')$origin,collapse=' '),
            ' ', template_file, ' ', tmpfl1
   )
   system(cmd1)

   tmpvol = mincGetVolume(tmpfl1) ; attr_tmpvol = attributes(tmpvol)
   tmpvol = allen_data ; attributes(tmpvol) = attr_tmpvol

   if (mask) {
      silentmincwrite(tmpvol,output.filename=tmpfl2)
      system(paste0(
            'mincmorph ', 
            ' -successive CCCCCCCDDDDB[0.9:1.1:1:0] ',
            ' ', tmpfl2, ' ', outfile
      ))
      ret = mincGetVolume(outfile)
   } else {
      ret = tmpvol
   }
   
   if (labels & save_file) {
      silentmincwrite(tmpvol,output.filename=outfile)
      warning('labels may have integers that cannot be writen to mincFile')
   }

   if (!labels & !mask & save_file) {
      silentmincwrite(tmpvol,output.filename=outfile)
   }

   unlink(tmpfl_nrrd)
   unlink(tmpfl1)
   if (mask) {unlink(tmpfl2)}
   if (!save_file) {unlink(outfile)}

   return(ret)
}

#' Download ABI gene templates
#'
#' Downloads and reads the ABI gene template with an option to write to a mincfile
#' @param age_dataset Choose the age of the template to download. Options are 'adult','P56','P28','P14','P4','E18.5','E16.5','E15.5','E13.5','E11.5'. Note 'adult' and 'P56' are aliases. 
#' @param outfile Optional outfile to write
#' @param labels If true, then the labels are downloaded. Default is FALSE, so the template is downloaded. Warning: labels may have integers that can't be written to a MINC file.
#' @param grid_annot If true, then the grid annotations are downloaded. They are much lower resolution compared to template. Default is FALSE, so the template is downloaded. Warning: grid annotations may have integers that can't be written to a MINC file.
#' @param binarize If true, then values are binarized with threshold of 0.5
#' @return template data as a 1D vector of class mincSingleDim
#' @import RMINC
allen_gene_template_download = function(
                                 age_dataset = c(
                                                  'adult','P56','P28','P14','P4',
                                                  'E18.5','E16.5','E15.5','E13.5','E11.5'
                                                ),
                                 outfile=NULL, labels = F, grid_annot = F, binarize = F
                                ) { 

   if (length(age_dataset) > 1) {stop('Only 1 option must be supplied for age_dataset')}

   age_dataset_opts = c(
      'adult','P56','P28','P14','P4',
      'E18.5','E16.5','E15.5','E13.5','E11.5'
    )
   
   if (!(age_dataset %in% age_dataset_opts)) {stop(age_dataset,'is not a valid option')}
   
   if (is.null(outfile)) {
      save_file=F
      outfile = tempfile(fileext='.mnc')
    } else {
      save_file=T
      if (file.exists(outfile)) stop(outfile,' is an existing file. Delete it before running.')
    }

    url_link = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_annotation/'

    if (age_dataset == 'adult') { 
       url_link = paste0(url_link,'P56') 
    } else {
       url_link = paste0(url_link,gsub('\\.','pt',age_dataset)) 
    }
   
   
   if (!labels & !grid_annot) {
      url_link = paste0(url_link,'_atlasVolume.zip')
      readcmd = function(url_link) { read.raw.atlas(url_link, url=T) }
      dtype = 'float'
   } else {
      if ( age_dataset %in% c('adult','P56') ) {
         url_link = paste0(url_link,'_Mouse')
      } else {
         url_link = paste0(url_link,'_DevMouse2012')
      }
      if (grid_annot) {
         url_link = paste0(url_link,'_gridAnnotation.zip')
         readcmd = function(url_link) { read.raw.gene(url_link, labels=T, url=T) }
         dtype = 'int'
      } else {
         url_link = paste0(url_link,'_annotation.zip')
         readcmd = function(url_link) { read.raw.gene(url_link, labels=T, url=T) }
         dtype = 'int'
      }
   }

   vol = readcmd(url_link)
   attr_vol = attributes(vol)
   
   attr_vol$ANTsRparam$spacing = attr_vol$ANTsRparam$spacing / 1000
   attr_vol$ANTsRparam$origin = c(-6.31 , -4.325 , -5.44)
   attributes(vol) = attr_vol

   allen_data = allenVectorTOmincVector(vol)
   
   if (binarize) {
      attr_tmpvol = attributes(allen_data)
      allen_data = as.integer( allen_data > 0.5 )
      attributes(allen_data) = attr_tmpvol
   }
   
   
   if (save_file) {
      template_file = system.file('extdata/Dorr_resampled_200um.mnc',package="ABIgeneRMINC")
      tmpfl1 = tempfile(fileext='.mnc')
      
      cmd1 = paste0(
               'mincresample -quiet ', 
               ifelse(dtype == 'int','-labels',''),' ',
               ' -nelements ',paste(attr(allen_data,'ANTsRparam')$imagesize,collapse=' '),
               ' -step ',paste(attr(allen_data,'ANTsRparam')$spacing,collapse=' '),
               ' -start ',paste(attr(allen_data,'ANTsRparam')$origin,collapse=' '),
               ' ', template_file, ' ', tmpfl1
      )
      system(cmd1)
   
      tmpvol = mincGetVolume(tmpfl1) ; attr_tmpvol = attributes(tmpvol)
      tmpvol = allen_data ; attributes(tmpvol) = attr_tmpvol

      silentmincwrite(tmpvol,output.filename=outfile)   
      if ( ( labels | grid_annot ) & !binarize) {   
         warning('labels may have integers that cannot be writen to mincFile')
      }

      unlink(tmpfl1)
   } 
   
   return(allen_data)
   
}

#' Silently write minc file
#'
#' @importFrom RMINC mincWriteVolume
silentmincwrite = function(...) {
   sink('/dev/null') ; mincWriteVolume(...) ; sink()
}


#' Pad or crop volume
#'
#' Pad or crop 2mm to volume on all sides
#' @param infile input MINC file
#' @param outfile outfile MINC file
#' @param operation If 'pad', then 2mm is added to all sides. If crop, then 2mm is cropped from all sides.
#' @importFrom RMINC mincConvertVoxelToWorld minc.separation.sizes minc.dimensions.sizes
mincvol_pad_crop = function(infile, outfile, operation = c('pad','crop')) {

   if (length(operation) != 1) {
      stop('operation must be either pad OR crop')
   } else if ( ! operation %in% c('pad','crop') ) {
      stop('operation must be either pad OR crop, not', operation)
   }

   starts = RMINC::mincConvertVoxelToWorld(infile,0,0,0)
   steps = rev(RMINC::minc.separation.sizes(infile))
   nelements = rev(RMINC::minc.dimensions.sizes(infile))

   if (operation == 'pad') {
      new_starts = starts - 2
      new_nelements = nelements + ( 4 / steps )
    } else if (operation == 'crop') {
      new_starts = starts + 2
      new_nelements = nelements - ( 4 / steps )
    }
   
   cmd1 = paste0(
            'mincresample -quiet', 
            ' -nelements ', paste(new_nelements,collapse=' '),
            ' -step ', paste(steps,collapse=' '),
            ' -start ', paste(new_starts,collapse=' '),
            ' ', infile, ' ', outfile
   )
   system(cmd1)
   
   return(outfile)
}

allen_grid_labels_download = function() {
   url_link = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_annotation/P56_Mouse_gridAnnotation.zip'
   brainmask = read.raw.gene(url_link,labels=TRUE,url=T)
   return(brainmask)
}

#' Resample volume to ABI space
#'
#' Resample a volume to a target with optional transform using linear interpolation.
#' @param source_volume MINC file or vector you want to transform
#' @param target_volume MINC file or vector you want to resample to. If NULL, it is assumed you want to transform to ABI CCFv3 50 micron. 
#' @param target_space_name Instead of specifying target_volume, you can provide one of several strings representing target volumes. Options include "adult gene expression" (i.e. 200um isotropic voxels), and "CCFv3 50um" (i.e. 50um isotropic voxels). If NULL, it is assumed you want to transform to ABI CCFv3 50 micron.
#' @param xfm_transform xfm file denoting the transform. Default NULL means no transformation
#' @param source_volume_resampled_file filepath to save resampled source_volume. If NULL (default), file is not saved. 
#' @return source_volume resampled to target space
#' @import RMINC
#' @export
ABI_resample = function(
     source_volume, 
     target_volume = NULL,
     target_space_name = NULL,
     xfm_transform = NULL,
     source_volume_resampled_file = NULL
 ) {

   if (!is.null(source_volume_resampled_file)) {
      outdir = dirname(source_volume_resampled_file)
   } else {
      outdir = tempdir()
   }
   tmpfl_source_res = tempfile(
                              pattern='source_res_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_source = tempfile(
                              pattern='source_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_target = tempfile(
                              pattern='target_',
                              tmpdir=outdir,
                              fileext='.mnc')
   tmpfl_source_res = tempfile(
                              pattern='source_res_',
                              tmpdir=outdir,
                              fileext='.mnc')

   if (is.null(source_volume_resampled_file)) {
      source_volume_resampled_file = tmpfl_source_res
   }

   if (is.character(source_volume)) {
      RMINC:::mincFileCheck(source_volume)
   } else {
      mincWriteVolume(source_volume, output.filename = tmpfl_source)
      source_volume = tmpfl_source
   }

   if (!is.null(xfm_transform)) {
     if (!is.character(xfm_transform)) {
        stop('xfm_transform argument must be a character denoting a valid file-path')
     }
     if (!file.exists(xfm_transform)) {
        stop('xfm_transform argument be an existing file-path')
     }
     tfm_str = paste0('-transform ',xfm_transform)
   } else {
     tfm_str = ""
   }

   target_vol_check = F
   if (is.character(target_volume)) {
      RMINC:::mincFileCheck(target_volume)
      likeVol = target_volume
      target_vol_check = T
   }

   if (!target_vol_check & !is.null(attributes(target_volume)$likeVolume)) {      
      mincWriteVolume(target_volume, output.filename = tmpfl_target)
      target_volume = tmpfl_target
      likeVol = attributes(target_volume)$likeVolume
      target_vol_check = T
   }
   
   if (!target_vol_check & !is.null(attributes(xfm_transform)$allen_ccf3v_template) ) {
      bool = tryCatch( {
         RMINC:::mincFileCheck(attributes(xfm_transform)$allen_ccf3v_template) 
         return(TRUE)
        }, error=function(cond) {return(FALSE)} )
      if (bool) {
         target_volume = attributes(xfm_transform)$allen_ccf3v_template
         likeVol = target_volume
         target_vol_check = T
      }
   }
      
   if (!target_vol_check) {
      target_volume = tmpfl_target
      
      if (is.null(target_space_name)) {
         x = allen_s2p_template_download(target_volume)
      } else {
         target_space_name_options = c('CCFv3 50um',"adult gene expression")      
         if (!(target_space_name %in% target_space_name_options)) {
            stop(target_space_name, 'is not valid or currently supported')
         }
         if (target_space_name == 'CCFv3 50um') {
            x = allen_s2p_template_download(target_volume)
         }
         if (target_space_name == 'adult gene expression') {         
            x = allen_gene_template_download(
                    age_dataset='adult',
                    outfile=target_volume,
                    grid_annot=T , binarize = T)
         }
      }
   }
   on.exit( {
      unlink(tmpfl_target)
      unlink(tmpfl_source)
      unlink(tmpfl_source_res)
   } )

   resample_cmd = paste0(
          'mincresample ',
          '-quiet -like ',target_volume,' ',
          tfm_str,' ',
          source_volume,' ',
          source_volume_resampled_file
   )

   x = system(resample_cmd, intern=T)
   
   ret = mincGetVolume(source_volume_resampled_file)
   
   return(ret)
   
}


