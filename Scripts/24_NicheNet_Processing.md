Processing of data for species integration
================
John Mariani
3/6/2023

## Load in Libraries

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(Seurat)
```

    ## Attaching SeuratObject

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(patchwork)
library(dplyr)
library(SeuratDisk)
```

    ## Registered S3 method overwritten by 'SeuratDisk':
    ##   method            from  
    ##   as.sparse.H5Group Seurat

``` r
library(biomaRt)
library(Matrix)


options(future.globals.maxSize = 16000 * 1024^2)

source("Scripts/HelperFunctions.R")
```

## Read in Human and Mouse Counts

``` r
mergedH <- readRDS("output/RDS/mergedH.rds")
mergedM <- readRDS("output/RDS/mergedM.rds")

# Remove ES Stage
mergedH <- subset(mergedH, subset = stage != "In Vitro - Pluripotent Stem Cell Stage")
table(mergedH$stage)
```

    ## 
    ## In Vitro - GPC Stage              In Vivo 
    ##                37805                 5337

``` r
table(mergedM$orig.ident)
```

    ## 
    ##        Sample_32        Sample_33        Sample_34   Sample_C27_Shi 
    ##             1493             1665             1852             1378 
    ## Sample_WA09_Shiv 
    ##             2010

## Read in gene names and ensembl IDs

``` r
humanFeatures <- read.delim("data_for_import/humanFeatures.txt")
mouseFeatures <- read.delim("data_for_import/mouseFeatures.txt")
```

## Subset out GPC4 so we can use all known human receptors/downstream genes in these cells

``` r
invitroInvivo <- readRDS("output/RDS/invitroInvivo.rds")
gpc4 <- subset(invitroInvivo, subset = cellType == "GPC4")

gpc4 <- subset(mergedH, cells = Cells(gpc4))
theRest <- subset(mergedH, subset = cellName %not in% Cells(gpc4))

dim(gpc4)
```

    ## [1] 39605  2728

## Make counts to edit. This seems safer than just merging and letting Seurat handle it

``` r
rawM <- mergedM@assays$RNA@counts
rawH <- theRest@assays$RNA@counts

dim(rawH)
```

    ## [1] 39605 40414

``` r
humanFeatures$seurat <- row.names(rawH)
mouseFeatures$seurat <- row.names(rawM)

#### Add If statements for these outputs

martH <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'http://apr2022.archive.ensembl.org/')
```

    ## Warning: Ensembl will soon enforce the use of https.
    ## Ensure the 'host' argument includes "https://"

``` r
ensemblGeneListHall <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type", "mmusculus_homolog_associated_gene_name"), filters = "ensembl_gene_id", values = humanFeatures$ensembl_gene_id, mart = martH)

#write.csv(ensemblGeneListHall, "data_for_import/NicheNet/ensemblGeneListHall.csv", quote = F, row.names = F)

ensemblGeneListH <- ensemblGeneListHall[ensemblGeneListHall$mmusculus_homolog_orthology_type == "ortholog_one2one",]
ensemblGeneListH <- merge(ensemblGeneListH, humanFeatures, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")

martM <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = 'http://apr2022.archive.ensembl.org/')
```

    ## Warning: Ensembl will soon enforce the use of https.
    ## Ensure the 'host' argument includes "https://"

``` r
ensemblGeneListMall <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_orthology_type", "hsapiens_homolog_associated_gene_name"), filters = "ensembl_gene_id", values = mouseFeatures$ensembl_gene_id, mart = martM)

#write.csv(ensemblGeneListMall, "data_for_import/NicheNet/ensemblGeneListMall", quote = F, row.names = F)
```

## H19 is a weird gene that is not a one to one ortholog but has the same name oddly…

``` r
sameNameH <- ensemblGeneListHall[ensemblGeneListHall$external_gene_name == ensemblGeneListHall$mmusculus_homolog_associated_gene_name,]
sameNameH <- sameNameH[sameNameH$external_gene_name != "",]

sameNameM <-  ensemblGeneListMall[ensemblGeneListMall$external_gene_name == ensemblGeneListMall$hsapiens_homolog_associated_gene_name,]
sameNameM <- sameNameM[sameNameM$external_gene_name != "",]

sameNameH2 <- humanFeatures[humanFeatures$external_gene_name %in% mouseFeatures$external_gene_name,]
```

## Make Dual Species object for all but GPC4

``` r
ensemblGeneListM <- ensemblGeneListMall[ensemblGeneListMall$hsapiens_homolog_orthology_type == "ortholog_one2one",]
ensemblGeneListM <- merge(ensemblGeneListM, mouseFeatures, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")

# These genes were filtered out because their biotypes differ between species and weren't any of the categories we kept
leftoutH <- ensemblGeneListH[ensemblGeneListH$ensembl_gene_id %not in% ensemblGeneListM$hsapiens_homolog_ensembl_gene,]
leftoutM <- ensemblGeneListM[ensemblGeneListM$ensembl_gene_id %not in% ensemblGeneListH$mmusculus_homolog_ensembl_gene,]


ensemblGeneListH <- ensemblGeneListH[ensemblGeneListH$ensembl_gene_id %in% ensemblGeneListM$hsapiens_homolog_ensembl_gene,]
ensemblGeneListM <- ensemblGeneListM[ensemblGeneListM$ensembl_gene_id %in% ensemblGeneListH$mmusculus_homolog_ensembl_gene,]


humanOne2One <- rawH[row.names(rawH) %in% ensemblGeneListH$seurat,]
mouseOne2One <- rawM[row.names(rawM) %in% ensemblGeneListM$seurat,]

rubric <- merge(ensemblGeneListH, ensemblGeneListM, by.x = "ensembl_gene_id", by.y = "hsapiens_homolog_ensembl_gene")

humanOne2One <- humanOne2One[match(rubric$seurat.x, row.names(humanOne2One)),]
mouseOne2One <- mouseOne2One[match(rubric$seurat.y, row.names(mouseOne2One)),]


head(row.names(humanOne2One))
```

    ## [1] "TSPAN6"   "TNMD"     "SCYL3"    "C1orf112" "FGR"      "FUCA2"

``` r
head(row.names(mouseOne2One))
```

    ## [1] "Tspan6"   "Tnmd"     "Scyl3"    "BC055324" "Fgr"      "Fuca2"

``` r
row.names(mouseOne2One) <- row.names(humanOne2One)

one2one <- cbind(humanOne2One, mouseOne2One)

humanOnly <- rawH[row.names(rawH) %not in% rubric$seurat.x,]
mouseOnly <- rawM[row.names(rawM) %not in% rubric$seurat.y,]

#Rename since they're apparently not orthologs
row.names(mouseOnly)[row.names(mouseOnly) == "H19"] <- "h19"

humanOnlyEmpty <- Matrix(nrow = nrow(humanOnly), ncol = ncol(mouseOnly), data = 0, sparse = TRUE)
colnames(humanOnlyEmpty) <- colnames(mouseOnly)

mouseOnlyEmpty <- Matrix(nrow = nrow(mouseOnly), ncol = ncol(humanOnly), data = 0, sparse = TRUE)
colnames(mouseOnlyEmpty) <- colnames(humanOnly)

humanOnly <- cbind(humanOnly, humanOnlyEmpty)
mouseOnly <- cbind(mouseOnlyEmpty, mouseOnly)

dim(one2one)
```

    ## [1] 16047 48812

``` r
dim(humanOnly)
```

    ## [1] 23558 48812

``` r
dim(mouseOnly)
```

    ## [1] 17992 48812

``` r
identical(colnames(one2one), colnames(humanOnly))
```

    ## [1] TRUE

``` r
identical(colnames(one2one), colnames(mouseOnly))
```

    ## [1] TRUE

``` r
dualSpecies <- rbind(one2one, humanOnly)
dualSpecies <- rbind(dualSpecies, mouseOnly)

dim(dualSpecies)
```

    ## [1] 57597 48812

``` r
dim(rawH)
```

    ## [1] 39605 40414

``` r
dim(rawM)
```

    ## [1] 34039  8398

## Reames cell names to species Human vs Mouse

``` r
colnames(dualSpecies)[1:ncol(rawH)] <- paste0(colnames(dualSpecies)[1:ncol(rawH)] , "H")

colnames(dualSpecies)[(ncol(rawH)+1):ncol(dualSpecies)] <- paste0(colnames(dualSpecies)[(ncol(rawH)+1):ncol(dualSpecies)], "M")

gpc4 <- gpc4@assays$RNA@counts

dim(gpc4)
```

    ## [1] 39605  2728

``` r
dim(dualSpecies)
```

    ## [1] 57597 48812

``` r
colnames(gpc4) <- paste0(colnames(gpc4), "H")
```

## Add GPC4 back with all human genes and add metadata

``` r
dualSpeciesSeurat <- CreateSeuratObject(dualSpecies)
dualSpeciesSeurat <- NormalizeData(dualSpeciesSeurat)
dim(dualSpeciesSeurat)
```

    ## [1] 57597 48812

``` r
gpc4Seurat <- CreateSeuratObject(gpc4)
gpc4Seurat <- NormalizeData(gpc4Seurat)

dualSpeciesOne2One <- subset(dualSpeciesSeurat, features =  row.names(one2one))
dim(dualSpeciesOne2One)
```

    ## [1] 16047 48812

``` r
finalNiche <- merge(dualSpeciesOne2One, gpc4Seurat)
dim(finalNiche)
```

    ## [1] 39605 51540

``` r
metaH <- mergedH@meta.data
metaM <- mergedM@meta.data

metaH$species <- "Human"
metaM$species <- "Mouse"

row.names(metaH) <- paste0(row.names(metaH), "H")
row.names(metaM) <- paste0(row.names(metaM), "M")

meta <- rbind(metaH, metaM)
meta <- meta[match(row.names(finalNiche@meta.data), row.names(meta)),]

mouse <- readRDS("output/RDS/mouse.rds")

invitroInvivoCellType <- data.frame(row.names = Cells(invitroInvivo), cellType = invitroInvivo$cellType)
mouseCellType <- data.frame(row.names = Cells(mouse), cellType = mouse$cellType)
mouseCellType$cellType <- paste0("Mouse ", mouseCellType$cellType)
row.names(invitroInvivoCellType) <- paste0(row.names(invitroInvivoCellType), "H")
row.names(mouseCellType) <- paste0(row.names(mouseCellType), "M")

cellTypes <- rbind(invitroInvivoCellType, mouseCellType)
cellTypes <- cellTypes[match(Cells(finalNiche), row.names(cellTypes)),, drop = F]

levels(cellTypes$cellType)
```

    ##  [1] "NPC"               "GPC1"              "GPC2"             
    ##  [4] "GPC3"              "GPC4"              "imOL"             
    ##  [7] "maOL"              "Astrocyte"         "Mouse maOL"       
    ## [10] "Mouse GPC"         "Mouse Astrocyte"   "Mouse imOL"       
    ## [13] "Mouse Microglia"   "Mouse Ependymal"   "Mouse Endothelial"
    ## [16] "Mouse Pericyte"    "Mouse NPC"         "Mouse Macrophage"

``` r
meta$cellType <- cellTypes$cellType

table(cellTypes$cellType)
```

    ## 
    ##               NPC              GPC1              GPC2              GPC3 
    ##             13784              8362             10826              3951 
    ##              GPC4              imOL              maOL         Astrocyte 
    ##              2728              1740              1004               747 
    ##        Mouse maOL         Mouse GPC   Mouse Astrocyte        Mouse imOL 
    ##              6763               152               151               278 
    ##   Mouse Microglia   Mouse Ependymal Mouse Endothelial    Mouse Pericyte 
    ##               321                50               483                39 
    ##         Mouse NPC  Mouse Macrophage 
    ##               148                13

``` r
finalNiche@meta.data <- meta
```

``` r
saveRDS(finalNiche, "output/RDS/finalNiche.rds")
```

``` r
sessionInfo()
```

    ## R version 4.2.3 (2023-03-15)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.2.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggplot2_3.4.4         Matrix_1.5-4          biomaRt_2.54.1       
    ##  [4] SeuratDisk_0.0.0.9020 patchwork_1.3.0.9000  devtools_2.4.5       
    ##  [7] usethis_2.1.6         SeuratObject_4.1.3    Seurat_4.3.0         
    ## [10] dplyr_1.1.1          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] BiocFileCache_2.6.1    plyr_1.8.8             igraph_2.0.3          
    ##   [4] lazyeval_0.2.2         sp_1.6-0               splines_4.2.3         
    ##   [7] listenv_0.9.0          scattermore_0.8        GenomeInfoDb_1.34.9   
    ##  [10] digest_0.6.31          htmltools_0.5.5        fansi_1.0.4           
    ##  [13] magrittr_2.0.3         memoise_2.0.1          tensor_1.5            
    ##  [16] cluster_2.1.4          ROCR_1.0-11            remotes_2.4.2         
    ##  [19] Biostrings_2.66.0      globals_0.16.2         matrixStats_0.63.0    
    ##  [22] spatstat.sparse_3.0-3  prettyunits_1.1.1      colorspace_2.1-0      
    ##  [25] rappdirs_0.3.3         blob_1.2.4             ggrepel_0.9.3         
    ##  [28] xfun_0.38              RCurl_1.98-1.12        callr_3.7.3           
    ##  [31] crayon_1.5.2           jsonlite_1.8.4         progressr_0.13.0      
    ##  [34] spatstat.data_3.0-4    survival_3.5-5         zoo_1.8-11            
    ##  [37] glue_1.6.2             polyclip_1.10-4        gtable_0.3.3          
    ##  [40] zlibbioc_1.44.0        XVector_0.38.0         leiden_0.4.3          
    ##  [43] pkgbuild_1.4.0         future.apply_1.10.0    BiocGenerics_0.44.0   
    ##  [46] abind_1.4-5            scales_1.3.0           DBI_1.1.3             
    ##  [49] spatstat.random_3.2-3  miniUI_0.1.1.1         Rcpp_1.0.10           
    ##  [52] progress_1.2.2         viridisLite_0.4.1      xtable_1.8-4          
    ##  [55] reticulate_1.34.0      bit_4.0.5              stats4_4.2.3          
    ##  [58] profvis_0.3.7          htmlwidgets_1.6.2      httr_1.4.5            
    ##  [61] RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3             
    ##  [64] urlchecker_1.0.1       pkgconfig_2.0.3        XML_3.99-0.14         
    ##  [67] farver_2.1.1           dbplyr_2.3.2           uwot_0.1.14           
    ##  [70] deldir_1.0-6           utf8_1.2.3             tidyselect_1.2.0      
    ##  [73] rlang_1.1.0            reshape2_1.4.4         later_1.3.0           
    ##  [76] AnnotationDbi_1.60.2   munsell_0.5.0          tools_4.2.3           
    ##  [79] cachem_1.0.7           cli_3.6.1              RSQLite_2.3.1         
    ##  [82] generics_0.1.3         ggridges_0.5.4         evaluate_0.20         
    ##  [85] stringr_1.5.0          fastmap_1.1.1          yaml_2.3.7            
    ##  [88] goftest_1.2-3          processx_3.8.0         knitr_1.42            
    ##  [91] bit64_4.0.5            fs_1.6.1               fitdistrplus_1.1-8    
    ##  [94] purrr_1.0.1            RANN_2.6.1             KEGGREST_1.38.0       
    ##  [97] pbapply_1.7-0          future_1.32.0          nlme_3.1-162          
    ## [100] mime_0.12              xml2_1.3.3             hdf5r_1.3.8           
    ## [103] compiler_4.2.3         rstudioapi_0.14        filelock_1.0.2        
    ## [106] curl_5.0.0             plotly_4.10.1          png_0.1-8             
    ## [109] spatstat.utils_3.1-0   tibble_3.2.1           stringi_1.7.12        
    ## [112] ps_1.7.4               lattice_0.21-8         vctrs_0.6.1           
    ## [115] pillar_1.9.0           lifecycle_1.0.3        spatstat.geom_3.2-9   
    ## [118] lmtest_0.9-40          RcppAnnoy_0.0.20       bitops_1.0-7          
    ## [121] data.table_1.14.8      cowplot_1.1.1          irlba_2.3.5.1         
    ## [124] httpuv_1.6.9           R6_2.5.1               promises_1.2.0.1      
    ## [127] KernSmooth_2.23-20     gridExtra_2.3          IRanges_2.32.0        
    ## [130] parallelly_1.35.0      sessioninfo_1.2.2      codetools_0.2-19      
    ## [133] MASS_7.3-58.3          pkgload_1.3.2          rprojroot_2.0.3       
    ## [136] withr_2.5.0            sctransform_0.3.5      GenomeInfoDbData_1.2.9
    ## [139] S4Vectors_0.36.2       hms_1.1.3              parallel_4.2.3        
    ## [142] grid_4.2.3             tidyr_1.3.0            rmarkdown_2.21        
    ## [145] Rtsne_0.16             spatstat.explore_3.2-7 Biobase_2.58.0        
    ## [148] shiny_1.7.4
