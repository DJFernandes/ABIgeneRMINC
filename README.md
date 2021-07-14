# ABIgeneRMINC

ABIgeneRMINC is an R library for integrating Allen Brain Institute gene expression data with MINC-format neuroanatomy studies.



## Installation

Easiest way to install is using devtools

```R
devtools::install_github('DJFernandes/ABIgeneRMINC')
```



## Dependencies

ABIgeneRMINC uses several tidyverse packages. It also works best with MINC toolkit version 1.9.18 (https://bic-mni.github.io/) and the RMINC package. 

## Usage

Neuroanatomy statistics can be conducted using the RMINC package, prior to being analysed using this package for gene expression patterns. 

```R
library(tidyverse)       # for convenient syntax
library(RMINC)           # for conducting neuroanatomy statistics
library(MRIcrotome)      # for visualisation
library(ABIgeneRMINC)    # for gene expression analysis
```



The following is data from Scholz et al. (2015) looking at the neuroanatomical effect of environmental enrichment. 

```R
# Anatomy file and mask. This is typically the consensus average in the study
anatfile_study = "/projects/egerek/jscholz/enriched/Exp2/hr_reg/hr_masks/template.mnc"
maskfile_study = "/projects/egerek/jscholz/enriched/Exp2/hr_reg/hr_masks/mask_dil.mnc"

# Statistic file. This is in the same space as the anatomy file whose numerical value at each voxel represents a quantity of statistical interest (i.e. t-statistic)
statsfile = "/projects/egerek/jscholz/enriched/Exp2/hr_stats/standard_vs_maze_0.2_tvalue-conditionMaze.mnc"

# Visualise the statistics to illustrate the effect of environmental enrichment on the mouse brain.
dev.new(width=7.0, height=4.5)
sliceSeries(nrow = 5, ncol=5, begin=54, end=254) %>%
       anatomy(
           mincArray(mincGetVolume(anatfile_study)), 
           low=300, high=5600) %>%
       overlay(
           mincArray(mincGetVolume(statsfile)), 
           low=2.0, high=6.5, symmetric = T) %>%
       anatomySliceIndicator(
           mincArray(mincGetVolume(anatfile_study)), 
           low=300, high=5600) %>% 
       legend("t-statistics") %>%
       draw()
```

![enrichment_effect](https://wiki.mouseimaging.ca/download/attachments/10650346/enrichment_effect.png?version=2&modificationDate=1626296216408&api=v2)



Before running gene expression analysis, we need to register the consensus average (defining the study space) to the Allen Gene Expression Atlas template from the Allen Brain Institute (ABI). The transformation mapping the consensus average to the ABI template is stored as an XFM file.

```R
# Register anatomy to ABI template
anatomy_transform_path = ABI_template_align(
    anatfile_study, maskfile_study,
    xfm_files_dir = "study_to_abi_transforms", 
    keep_files = T
) 

anatomy_transform_path
# "study_to_abi_transforms/anatomy_to_abi.xfm"
```



There is also an image created to evaluate the quality of the registration. Contours of the ABI template are mapped onto the resampled consensus average. 

![check_registration_quality](https://wiki.mouseimaging.ca/download/attachments/10650346/check_registration_quality.png?version=1&modificationDate=1626296283407&api=v2)



The statistics can now be transformed to the gene expression atlas space, which is the final step before conducting gene expression analysis. 

```R
resampled_statistics = ABI_resample(
         anatomy = statsfile, 
         anatomy_transform = anatomy_transform_path, 
         ABI_statistics_file = 'resampled_statistics.mnc',
         keep_files = T
       ) 
```

Note that we can save the resampled statistics in a file (ex. 'resampled_statistics.mnc') so that, in the future, we can simply start the gene expression analysis from here. 

Let us run the gene expression analysis. We will set the target threshold to be 2 (i.e. regions with statistics above 2 are part of the ROI). We will also consider the statistics to be symmetric (i.e. -2 and below is also part of the ROI). We will use the whole brain as the contrast region. We will also run the jobs locally and in parallel on 4 cores. 

```R
gene_fold_change = adult_gene_expression_analysis(
   anatomy_statistics = 'resampled_statistics.mnc' ,
   gene_expression_paths = NULL                    ,
   target_threshold = 2.0                          ,
   contrast_threshold = NULL                       ,
   symmetric_statistics = T                        ,
   parallel = c('local',4)
)
```



To run Gene Ontology Enrichment Analysis, use the GOrilla_run function. In this example, we will run ontology on the top 7000 genes with the highest fold-change. 

```
background_genes = gene_fold_change %>% arrange(desc(fold_change)) %>% pull(acronym)
target_genes = background_genes[1:7000]

gene_ontology_enrichment_df = GOrilla_run(
     target.genes=target_genes,
     background.genes=background_genes
  )
```

![GO](https://wiki.mouseimaging.ca/download/attachments/10650346/GO.png?version=1&modificationDate=1626296493857&api=v2)

You can also find human diseases enriched in these mouse genes using the disease_ontology_enrichment_analysis function.

```R
disease_ontology_enrichment_df = disease_ontology_enrichment_analysis(
     target.genes=target_genes,
     background.genes=background_genes
  )

disease_ontology_enrichment_df
# # A tibble: 11,703 x 6
#    disease_id disease                     p        q `Enrichment (N, B,… Genes  
#    <chr>      <chr>                   <dbl>    <dbl> <chr>               <list> 
#  1 C0596887   mathematical abili…   9.36e-8 0.000583 1.29 (12174,679,39… <chr […
#  2 C0557874   Global development…   9.96e-8 0.000583 1.19 (12174,1434,3… <chr […
#  3 C1305855   Body mass index       2.09e-7 0.000817 1.25 (12174,843,39… <chr […
#  4 C0454644   Delayed speech and…   5.47e-7 0.00160  1.31 (12174,517,39… <chr […
#  5 C4048268   Cortical visual im…   7.13e-7 0.00167  1.68 (12174,115,39… <chr […
#  6 C3810365   Central visual imp…   1.29e-6 0.00252  1.72 (12174,98,397… <chr […
#  7 C1858120   Generalized hypoto…   3.97e-6 0.00664  1.21 (12174,899,39… <chr […
#  8 C0344482   Hypoplasia of corp…   4.97e-6 0.00708  1.35 (12174,335,39… <chr […
#  9 C0232466   Feeding difficulti…   5.94e-6 0.00708  1.3 (12174,445,397… <chr […
# 10 C0856975   Autistic behavior     6.05e-6 0.00708  1.44 (12174,222,39… <chr […
# # … with 11,693 more rows
```



## Contributing
Suggestions, concerns, and pull requests are welcome. 



## Citation

Fernandes, Darren J., et al. "Spatial gene expression analysis of neuroanatomical differences in mouse models." *Neuroimage* 163 (2017): 220-230.
