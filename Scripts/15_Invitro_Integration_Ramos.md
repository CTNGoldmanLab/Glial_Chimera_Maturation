Integration of in vitro hGPCs and Fetal snRNA-seq
================
John Mariani
12/6/2022

``` r
library(Seurat)
library(scPlottingTools)
library(ggplot2)
library(tidyr)
library(dplyr)
library(MAST)
library(plyr)
library(xlsx)
library(patchwork)
library(ggplot2)
library(scales)
library(ggVennDiagram)
library(data.table)
library(slingshot)
library(tradeSeq)
library(tidyr)
library(dplyr)
library(plyr)
library(magrittr)
library(viridis)
library(tidyr)
library(EnhancedVolcano)
library(ggalluvial)



axisTitleSize <- 20
axisTextSize <- 18
labelFont = 18
titleFont = 24
tagSize = 34

source("Scripts/HelperFunctions.R")


theme_manuscript <-  theme(axis.text = element_text(size = axisTextSize), 
        axis.title = element_text(size = axisTitleSize), 
        title = element_text(size = titleFont), 
        legend.title = element_text(size = axisTitleSize),
        legend.text = element_text(size = axisTitleSize),
        plot.tag = element_text(size = tagSize),
        plot.title = element_text(size = titleFont))

manuscriptPalette <- c("In Vivo" = "red2", 
                       "In Vitro - GPC Stage" = "#2E30FF",
                       "NPC" = "magenta",
                       "GPC1" = "forestgreen",
                       "GPC2" = "darkorange",
                       "GPC3" = "firebrick2",
                       "GPC4" = "turquoise",
                       "Astrocyte" = "dodgerblue2",
                       "imOL" = "gold",
                       "maOL" = "darkorchid4",
                       "GPC" = "turquoise",
                       "imAstrocyte" = "firebrick2",
                       "cGPC" = "darkorange",
                       "cAPC" = "forestgreen")
```

## Import and filter Fetal snRNA-seq data

``` r
sampleList <- list.files("data_for_import/External/Ramos/Matrices/")

sampleData <- read.csv("data_for_import/External/Ramos/Metadata/GlialMetadata.csv")

unique(sampleData$sample)
```

    ##  [1] "9G"  "11G" "24G" "26G" "23G" "30G" "34G" "56G" "63G" "3G"  "62G" "64G"
    ## [13] "6G"  "10G" "60G" "3C"  "62C" "64C" "6C"  "10C" "60C" "2C"  "20C" "31C"
    ## [25] "2G"  "20G" "31G"

``` r
sampleDF <- data.frame(sample = c("9C", "9G",
                      "24C", "24G",
                      "26C", "26G",
                      "11C", "11G",
                      "56C", "56G",
                      "34C", "34G",
                      "30C", "30G",
                      "63C", "63G",
                      "23C", "23G",
                      "62C", "62G",
                      "64C", "64G",
                      "3C", "3G",
                      "60C", "60G",
                      "10C", "10G",
                      "6C", "6G",
                      "31C", "31G",
                      "2C", "2G",
                      "20C", "20G"), GEO = unique(sort(sampleList)))

sampleDF[sampleDF$sample %not in% sampleData$sample,]
```

    ##    sample               GEO
    ## 1      9C Sample_GSM6720852
    ## 3     24C Sample_GSM6720854
    ## 5     26C Sample_GSM6720856
    ## 7     11C Sample_GSM6720858
    ## 9     56C Sample_GSM6720860
    ## 11    34C Sample_GSM6720862
    ## 13    30C Sample_GSM6720864
    ## 15    63C Sample_GSM6720866
    ## 17    23C Sample_GSM6720868

``` r
sampleData[sampleData$sample %not in% sampleDF$sample,]
```

    ##  [1] cellID          sample          nCount_RNA      nFeature_RNA   
    ##  [5] nCount_SCT      nFeature_SCT    percent.mt      seurat_clusters
    ##  [9] celltypes       Type           
    ## <0 rows> (or 0-length row.names)

``` r
rawH <- sapply(sampleList, function(x) {print(x) ; Read10X(paste0("data_for_import/External/Ramos/Matrices/",x, "/"))})
```

    ## [1] "Sample_GSM6720852"
    ## [1] "Sample_GSM6720853"
    ## [1] "Sample_GSM6720854"
    ## [1] "Sample_GSM6720855"
    ## [1] "Sample_GSM6720856"
    ## [1] "Sample_GSM6720857"
    ## [1] "Sample_GSM6720858"
    ## [1] "Sample_GSM6720859"
    ## [1] "Sample_GSM6720860"
    ## [1] "Sample_GSM6720861"
    ## [1] "Sample_GSM6720862"
    ## [1] "Sample_GSM6720863"
    ## [1] "Sample_GSM6720864"
    ## [1] "Sample_GSM6720865"
    ## [1] "Sample_GSM6720866"
    ## [1] "Sample_GSM6720867"
    ## [1] "Sample_GSM6720868"
    ## [1] "Sample_GSM6720869"
    ## [1] "Sample_GSM6720870"
    ## [1] "Sample_GSM6720871"
    ## [1] "Sample_GSM6720872"
    ## [1] "Sample_GSM6720873"
    ## [1] "Sample_GSM6720874"
    ## [1] "Sample_GSM6720875"
    ## [1] "Sample_GSM6720876"
    ## [1] "Sample_GSM6720877"
    ## [1] "Sample_GSM6720878"
    ## [1] "Sample_GSM6720879"
    ## [1] "Sample_GSM6720880"
    ## [1] "Sample_GSM6720881"
    ## [1] "Sample_GSM6720882"
    ## [1] "Sample_GSM6720883"
    ## [1] "Sample_GSM6720884"
    ## [1] "Sample_GSM6720885"
    ## [1] "Sample_GSM6720886"
    ## [1] "Sample_GSM6720887"

``` r
sum(rawH[[1]][,2])
```

    ## [1] 1

``` r
mean(colSums(rawH[[1]]))
```

    ## [1] 10.67638

``` r
for(i in 1:length(rawH)){
  rawH[[i]] <- rawH[[i]][,colSums(rawH[[i]]) > 1000]
}

for(i in 1:length(rawH)){ 
  colnames(rawH[[i]]) <- paste0(colnames(rawH[[i]]),"_",sampleList[i])
}


sampleData$cellID <- sub("_.*", "", sampleData$cellID)

sampleData <- merge(sampleData, sampleDF, by.x = "sample", by.y = "sample")

sampleData$cellID2 <- paste0(sampleData$cellID, "_", sampleData$GEO)

seurat.objectsH <- sapply(c(1:length(rawH)), function(x) CreateSeuratObject(rawH[[x]], project = sampleList[x]))

```

``` r
mergedH <- merge(seurat.objectsH[[1]], y = seurat.objectsH[2:length(seurat.objectsH)])

mergedH$cellID2 <- Cells(mergedH)

mergedHglia <- subset(mergedH, subset = cellID2 %in% sampleData$cellID2)

mergedHgliaMeta <- mergedHglia@meta.data
mergedHgliaMeta <- merge(mergedHgliaMeta, sampleData[,c("cellID2", "sample", "celltypes", "Type", "GEO", "percent.mt", "seurat_clusters")], by = "cellID2")


mergedHgliaMeta <- mergedHgliaMeta[match(Cells(mergedHglia), mergedHgliaMeta$cellID2),]

identical(mergedHgliaMeta$cellID2, Cells(mergedHglia))
```

    ## [1] TRUE

``` r
row.names(mergedHgliaMeta) <- mergedHgliaMeta$cellID2

mergedHglia@meta.data <- mergedHgliaMeta

Idents(mergedHglia) <- mergedHglia$celltypes
```

``` r
#Filter to fetal data
external <- subset(mergedHglia, subset = Type %in% c("CorticalPlate", "GerminalMatrix"))
external <- subset(external, subset = celltypes %not in% c("EPD", "CP"))

### clean up memory
rm(mergedHglia)
rm(rawH)
rm(seurat.objectsH)
rm(mergedH)
gc()
```

    ##             used   (Mb) gc trigger    (Mb) limit (Mb)   max used  (Mb)
    ## Ncells   8144141  435.0   24821835  1325.7         NA   26475662  1414
    ## Vcells 171069665 1305.2 3475714303 26517.6     102400 4344630994 33147

``` r
# Load and in vitro data
invitro <- readRDS("output/RDS/invitroInvivo.rds")
invitro <- subset(invitro, subset = stage ==  "In Vitro - GPC Stage")
invitro <- subset(invitro, subset = cellType %not in% c("Astrocyte", "maOL"))

DefaultAssay(external) <- "RNA"
external <- NormalizeData(external)

DefaultAssay(invitro) <- "RNA"
invitro <- NormalizeData(invitro)
```

## Subset to common genes and then scTransform datasets individually

``` r
table(external$orig.ident)
```

    ## 
    ## Sample_GSM6720853 Sample_GSM6720855 Sample_GSM6720857 Sample_GSM6720859 
    ##               628              1487               392              1061 
    ## Sample_GSM6720861 Sample_GSM6720863 Sample_GSM6720865 Sample_GSM6720867 
    ##              1011               341               892               602 
    ## Sample_GSM6720869 Sample_GSM6720870 Sample_GSM6720871 Sample_GSM6720872 
    ##               771               460               839               616 
    ## Sample_GSM6720873 Sample_GSM6720874 Sample_GSM6720875 Sample_GSM6720876 
    ##              2658               497               333               888 
    ## Sample_GSM6720877 Sample_GSM6720878 Sample_GSM6720879 Sample_GSM6720880 
    ##               935               523               758               971 
    ## Sample_GSM6720881 
    ##              5280

``` r
table(invitro$orig.ident)
```

    ## 
    ##     Sample_C27_CD140_3      Sample_C27_DAPI_3          Sample_CD-140 
    ##                   9723                   7084                   2389 
    ##    Sample_CD140_WA09_5        Sample_unsorted Sample_Unsorted_WA09_5 
    ##                   8202                   2008                   8397

``` r
genes_external <- rownames(external)
genes_invitro <- rownames(invitro)

common_genes <- intersect(genes_external, genes_invitro)

externalFilt<- subset(external, features = common_genes)

invitroFilt <- subset(invitro, features = common_genes)

mergedBoth <- merge(externalFilt, invitroFilt)
mergedBoth <- SplitObject(mergedBoth, split.by = "orig.ident")
mergedBoth <- lapply(mergedBoth, SCTransform, verbose = T)
```

``` r
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
```

   
``` r
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
```

    ## Computing nearest neighbor graph
    ## Computing SNN

``` r
integrated <- FindClusters(integrated, resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 59746
    ## Number of edges: 3078364
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9120
    ## Number of communities: 22
    ## Elapsed time: 15 seconds

    ## 1 singletons identified. 21 final clusters.

``` r
integrated$stage <- ifelse(!is.na(integrated$stage), yes = "In Vitro - GPC Stage", no = "Fetal")
```

## Predict labels from external data set and transfer to our hGPC data

``` r
invitroInt <- subset(integrated, subset = stage != "Fetal")
externalInt <- subset(integrated, subset = stage == "Fetal")


(DimPlotCustom(invitroInt, group.by = "cellType", label = T) + ggtitle("In vitro - GPC Stage")) | 
  DimPlotCustom(externalInt, group.by = "celltypes", label = T) + ggtitle("Fetal - Ramos et al")
```

![](Invitro_Integration_Ramos_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
anchors <- FindTransferAnchors(
  reference = externalInt,
  query = invitroInt,    
  reference.assay = "integrated",   
  query.assay = "integrated",
  reduction = "rpca",                
  dims = 1:30                        
)
```

    ## Performing PCA on the provided reference using 2000 features as input.

    ## Centering and scaling data matrix

    ## Performing PCA on the provided query using 2000 features as input.

    ## Centering and scaling data matrix

    ## Projecting new data onto SVD
    ## Projecting new data onto SVD

    ## Projecting cell embeddings

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 4764 anchors

    ## Filtering anchors

    ##  Retained 2541 anchors

``` r
predictions <- TransferData(
  anchorset = anchors,
  refdata = externalInt$celltypes,
  query = invitroInt,
  dims = 1:30,
  weight.reduction = "rpca.ref"
)
```

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Predicting cell labels

    ## Warning: Layer counts isn't present in the assay object; returning NULL

``` r
table(predictions$predicted.id, predictions$cellType)
```

    ##         
    ##          GPC1 GPC2 GPC3 GPC4 imOL  NPC
    ##   AC        6   18  669    3    0    6
    ##   AC-p    128   19  237    6    0  412
    ##   EN        0    0    0    1    0  366
    ##   gIPC     25  941  996   12    0   31
    ##   gIPC-A    0   42   23    0    0    0
    ##   gIPC-O    4   18  129   68    0    0
    ##   IN       96    7    3   20    0 2468
    ##   mIPC   3465 3323  443   55    0 4002
    ##   nIPC   4457 6145   61   86    0 5662
    ##   OPC      10    9    8  601    1    9
    ##   oRG       2  154  147    6    0    0
    ##   preOL     0    0    0   20   26   20
    ##   RG/AC   169  130  858    2    0  806
    ##   tRG       0    0  371    0    0    1

``` r
integrated <- AddMetaData(
  object = integrated,
  metadata = predictions$predicted.id,
  col.name = "predicted.id"
)

round((prop.table(table(integrated$predicted.id, integrated$cellType), margin = 2)*100),1)
```

    ##         
    ##          GPC1 GPC2 GPC3 GPC4 imOL  NPC
    ##   AC      0.1  0.2 17.0  0.3  0.0  0.0
    ##   AC-p    1.5  0.2  6.0  0.7  0.0  3.0
    ##   EN      0.0  0.0  0.0  0.1  0.0  2.7
    ##   gIPC    0.3  8.7 25.2  1.4  0.0  0.2
    ##   gIPC-A  0.0  0.4  0.6  0.0  0.0  0.0
    ##   gIPC-O  0.0  0.2  3.3  7.7  0.0  0.0
    ##   IN      1.1  0.1  0.1  2.3  0.0 17.9
    ##   mIPC   41.4 30.8 11.2  6.2  0.0 29.0
    ##   nIPC   53.3 56.9  1.5  9.8  0.0 41.1
    ##   OPC     0.1  0.1  0.2 68.3  3.7  0.1
    ##   oRG     0.0  1.4  3.7  0.7  0.0  0.0
    ##   preOL   0.0  0.0  0.0  2.3 96.3  0.1
    ##   RG/AC   2.0  1.2 21.7  0.2  0.0  5.8
    ##   tRG     0.0  0.0  9.4  0.0  0.0  0.0

``` r
integrated <- AddMetaData(
  object = integrated,
  metadata = predictions$predicted.id.score, 
  col.name = "predicted.score"
)

saveRDS(integrated, "output/RDS/integratedRamos.rds")
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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggalluvial_0.12.5           EnhancedVolcano_1.16.0     
    ##  [3] ggrepel_0.9.3               viridis_0.6.2              
    ##  [5] viridisLite_0.4.1           magrittr_2.0.3             
    ##  [7] tradeSeq_1.12.0             slingshot_2.6.0            
    ##  [9] TrajectoryUtils_1.6.0       princurve_2.1.6            
    ## [11] data.table_1.14.8           ggVennDiagram_1.2.2        
    ## [13] scales_1.3.0                patchwork_1.3.0.9000       
    ## [15] xlsx_0.6.5                  plyr_1.8.8                 
    ## [17] MAST_1.24.1                 SingleCellExperiment_1.20.1
    ## [19] SummarizedExperiment_1.28.0 Biobase_2.58.0             
    ## [21] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
    ## [23] IRanges_2.32.0              S4Vectors_0.36.2           
    ## [25] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
    ## [27] matrixStats_0.63.0          dplyr_1.1.1                
    ## [29] tidyr_1.3.0                 ggplot2_3.4.4              
    ## [31] scPlottingTools_0.0.0.9000  SeuratObject_5.0.2         
    ## [33] Seurat_4.3.0               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] spam_2.10-0            igraph_2.0.3           lazyeval_0.2.2        
    ##   [4] sp_1.6-0               splines_4.2.3          BiocParallel_1.32.6   
    ##   [7] listenv_0.9.0          scattermore_0.8        digest_0.6.31         
    ##  [10] htmltools_0.5.5        fansi_1.0.4            tensor_1.5            
    ##  [13] cluster_2.1.4          ROCR_1.0-11            limma_3.54.2          
    ##  [16] globals_0.16.2         R.utils_2.12.2         spatstat.sparse_3.0-3 
    ##  [19] RVenn_1.1.0            colorspace_2.1-0       xfun_0.38             
    ##  [22] RCurl_1.98-1.12        jsonlite_1.8.4         progressr_0.13.0      
    ##  [25] spatstat.data_3.0-4    survival_3.5-5         zoo_1.8-11            
    ##  [28] glue_1.6.2             polyclip_1.10-4        gtable_0.3.3          
    ##  [31] zlibbioc_1.44.0        XVector_0.38.0         leiden_0.4.3          
    ##  [34] DelayedArray_0.24.0    future.apply_1.10.0    abind_1.4-5           
    ##  [37] edgeR_3.40.2           DBI_1.1.3              spatstat.random_3.2-3 
    ##  [40] miniUI_0.1.1.1         Rcpp_1.0.10            xtable_1.8-4          
    ##  [43] reticulate_1.43.0      dotCall64_1.1-1        htmlwidgets_1.6.2     
    ##  [46] httr_1.4.5             RColorBrewer_1.1-3     ellipsis_0.3.2        
    ##  [49] ica_1.0-3              R.methodsS3_1.8.2      pkgconfig_2.0.3       
    ##  [52] rJava_1.0-6            farver_2.1.1           uwot_0.1.14           
    ##  [55] deldir_1.0-6           locfit_1.5-9.7         utf8_1.2.3            
    ##  [58] labeling_0.4.2         tidyselect_1.2.0       rlang_1.1.0           
    ##  [61] reshape2_1.4.4         later_1.3.0            munsell_0.5.0         
    ##  [64] tools_4.2.3            cli_3.6.1              generics_0.1.3        
    ##  [67] ggridges_0.5.4         evaluate_0.20          stringr_1.5.0         
    ##  [70] fastmap_1.1.1          yaml_2.3.7             goftest_1.2-3         
    ##  [73] knitr_1.42             fitdistrplus_1.1-8     purrr_1.0.1           
    ##  [76] RANN_2.6.1             pbapply_1.7-0          future_1.32.0         
    ##  [79] nlme_3.1-162           mime_0.12              R.oo_1.25.0           
    ##  [82] compiler_4.2.3         rstudioapi_0.14        plotly_4.10.1         
    ##  [85] png_0.1-8              spatstat.utils_3.1-0   tibble_3.2.1          
    ##  [88] stringi_1.7.12         highr_0.10             lattice_0.21-8        
    ##  [91] Matrix_1.6-4           vctrs_0.6.1            pillar_1.9.0          
    ##  [94] lifecycle_1.0.3        spatstat.geom_3.2-9    lmtest_0.9-40         
    ##  [97] RcppAnnoy_0.0.20       cowplot_1.1.1          bitops_1.0-7          
    ## [100] irlba_2.3.5.1          httpuv_1.6.9           R6_2.5.1              
    ## [103] promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3         
    ## [106] parallelly_1.35.0      codetools_0.2-19       MASS_7.3-58.3         
    ## [109] xlsxjars_0.6.1         rprojroot_2.0.3        withr_2.5.0           
    ## [112] sctransform_0.3.5      GenomeInfoDbData_1.2.9 mgcv_1.8-42           
    ## [115] parallel_4.2.3         grid_4.2.3             rmarkdown_2.21        
    ## [118] Rtsne_0.16             spatstat.explore_3.2-7 shiny_1.7.4
