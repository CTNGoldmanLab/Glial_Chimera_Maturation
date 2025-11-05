Integration of in vitro hGPCs and Fetal scRNA-seq
================
John Mariani
07/30/2025

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

## Import Fetal and in vitro scRNA-seq data

``` r
invitroInvivo <- readRDS("output/RDS/invitroInvivo.rds")

invitro <- subset(invitroInvivo, subset = stage == "In Vitro - GPC Stage")
invitro <- subset(invitro, subset = cellType %in% c("NPC", "GPC1", "GPC2", "GPC3", "GPC4", "imOL"))
invitro$cellType <- droplevels(invitro$cellType)


external <- readRDS("data_for_import/External/VanBruggen/scRNA_human_dev.rds")
external <- UpdateSeuratObject(external)
```

``` r
table(external$Clusters)
```

    ## 
    ##                                 Cortical Interneurons 
    ##                                                  1449 
    ##                                    Cortical Pyramidal 
    ##                                                  1241 
    ##                                           Endothelial 
    ##                                                    58 
    ##                             Excitatory neurons cortex 
    ##                                                  1741 
    ##                           Excitatory neurons midbrain 
    ##                                                  1531 
    ##                  Excitatory neurons possibly midbrain 
    ##                                                  1449 
    ##         Forebrain early neuroblast possibly GABAergic 
    ##                                                  1474 
    ##                       Forebrain inhibitory neuroblast 
    ##                                                  1230 
    ##                      Forebrain neural progenitor EMX1 
    ##                                                  1400 
    ##                                   GABAergic forebrain 
    ##                                                   847 
    ## GABAergic or interneuron neuroblast probably midbrain 
    ##                                                   351 
    ##                                             Glioblast 
    ##                                                  1343 
    ##                                     Glioblast/Pre-OPC 
    ##                                                   732 
    ##                                  Hindbrain neuroblast 
    ##                                                   428 
    ##       Inhibitory neurons Midbrain, possibly GABAergic 
    ##                                                  1736 
    ##                                             Microglia 
    ##                                                    43 
    ##                              Mid/Hindbrain neuroblast 
    ##                                                   812 
    ##                        Midbrain inhibitory neuroblast 
    ##                                                    89 
    ##              Midbrain/Hindbrain inhibitory neuroblast 
    ##                                                    15 
    ##                     Neuroblast motorneuron/GABAergic? 
    ##                                                   663 
    ##                                                  OPCs 
    ##                                                    98 
    ##                                           Radial Glia 
    ##                                                   844 
    ##                                   Radial Glia cycling 
    ##                                                   957 
    ##                       Radial Glia potential glioblast 
    ##                                                  1736 
    ##                              Radial Glia VLMC primed? 
    ##                                                   502 
    ##                                 Radial Glia/Glioblast 
    ##                                                  1336 
    ##       Radial glia/Glioblast/Forebrain progenitor EMX1 
    ##                                                   155 
    ##                             Striatum/Cortical neurons 
    ##                                                   300 
    ##                                                     U 
    ##                                                   504 
    ##                                                 VLMCs 
    ##                                                    97

``` r
external <- subset(external, subset = Clusters %not in% c("Endothelial", "Microglia", "U", "VLMCs", "Radial Glia VLMC primed?"))
external$Clusters <- droplevels(external$Clusters)

DefaultAssay(external) <- "RNA"
external <- NormalizeData(external)

DefaultAssay(invitro) <- "RNA"
invitro <- NormalizeData(invitro)
```

## Subset to common genes and then scTransform datasets individually

``` r
genes_external <- rownames(external)
genes_invitro <- rownames(invitro)

common_genes <- intersect(genes_external, genes_invitro)

externalFilt<- subset(external, features = common_genes)

invitroFilt <- subset(invitro, features = common_genes)

externalFilt$orig.ident <- externalFilt$Sample

mergedBoth <- merge(externalFilt, invitroFilt)
mergedBoth <- SplitObject(mergedBoth, split.by = "orig.ident")
mergedBoth <- lapply(mergedBoth, SCTransform, verbose = T)
```

## Integrate data across captures

``` r
features <- SelectIntegrationFeatures(object.list = mergedBoth)

mergedBoth <- PrepSCTIntegration(object.list = mergedBoth, anchor.features = features)

mergedBoth <- lapply(X = mergedBoth, FUN = RunPCA, features = features)
```

``` r
anchors <- FindIntegrationAnchors(object.list = mergedBoth, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
```

 
``` r
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
```


``` r
# Integrate the data
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

``` r
integrated$stage <- ifelse(!is.na(integrated$stage), yes = "In Vitro - GPC Stage", no = "Fetal")

#saveRDS(integrated, "output/RDS/integratedVanBrueggen.rds")
```

## Predict labels from external data set and transfer to our hGPC data

``` r
#integrated <- readRDS("output/RDS/integratedVanBrueggen.rds")
integrated$cellType <- factor(integrated$cellType, c("NPC", "GPC1", "GPC2", "GPC3", "GPC4", "imOL"))

groupDF <- data.frame(Level1 = unique(integrated$Level1)[1:7], group = NA)
groupDF$Level1
```

    ## [1] "NeuralProgenitor" "EarlyInhibitory"  "RadialGlia"       "EarlyExcitatory" 
    ## [5] "Glioblast"        "OPC"              "31"

``` r
groupDF$group <- c("NPC", "Early IN", "RG", "Early EN", "Glioblast", "OPC", "Early IN")

integrated$group <- plyr::mapvalues(integrated$Level1, from = groupDF$Level1, to = groupDF$group)

invitroInt <- subset(integrated, subset = stage != "Fetal")
externalInt <- subset(integrated, subset = stage == "Fetal")



((DimPlotCustom(externalInt, group.by = "group", label = T) + ggtitle("8-10wk GA Fetal - Van Brueggen et al 2022") + guides(fill = guide_legend(override.aes = list(size = 5)))) |
(DimPlotCustom(invitroInt, group.by = "cellType", label = T) + ggtitle("In vitro - GPC Stage") + guides(fill = guide_legend(override.aes = list(size = 5))))) & theme_manuscript & labs(fill = "Cell Type")
```

![](Invitro_Integration_van_Bruggen_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
DefaultAssay(externalInt)
```

    ## [1] "integrated"

``` r
DefaultAssay(invitroInt)
```

    ## [1] "integrated"

``` r
DefaultAssay(integrated) <- "integrated"



anchors <- FindTransferAnchors(
  reference = externalInt, # Use the full integrated object as the source for anchor finding
  query = invitroInt,     # Use the full integrated object as the target
  reference.assay = "integrated",    # Specify the assay used for integration
  query.assay = "integrated",
  reduction = "rpca",                 # Or "rpca", depending on your integration method
  dims = 1:30                        # Number of dimensions to use
)
```

``` r
##

predictions <- TransferData(
  anchorset = anchors,
  refdata = externalInt$group, # This should be the cell type column from your original dataset1
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
    ##              NPC GPC1 GPC2 GPC3 GPC4 imOL
    ##   Early EN   559    9   11   14    5    1
    ##   Early IN  9294 6816 2693   60   64    0
    ##   Glioblast  137   69 1473 3370   44    0
    ##   NPC       1867  370  121   39  128    0
    ##   OPC         25    2    7   57  592   26
    ##   RG        1901 1096 6501  405   47    0

``` r
round(prop.table(table(predictions$predicted.id, predictions$cellType), margin = 2)*100,1)
```

    ##            
    ##              NPC GPC1 GPC2 GPC3 GPC4 imOL
    ##   Early EN   4.1  0.1  0.1  0.4  0.6  3.7
    ##   Early IN  67.4 81.5 24.9  1.5  7.3  0.0
    ##   Glioblast  1.0  0.8 13.6 85.4  5.0  0.0
    ##   NPC       13.5  4.4  1.1  1.0 14.5  0.0
    ##   OPC        0.2  0.0  0.1  1.4 67.3 96.3
    ##   RG        13.8 13.1 60.2 10.3  5.3  0.0

``` r
integrated <- AddMetaData(
  object = integrated,
  metadata = predictions$predicted.id,
  col.name = "predicted.id"
)

round((prop.table(table(integrated$predicted.id, integrated$cellType), margin = 2)*100),1)
```

    ##            
    ##              NPC GPC1 GPC2 GPC3 GPC4 imOL
    ##   Early EN   4.1  0.1  0.1  0.4  0.6  3.7
    ##   Early IN  67.4 81.5 24.9  1.5  7.3  0.0
    ##   Glioblast  1.0  0.8 13.6 85.4  5.0  0.0
    ##   NPC       13.5  4.4  1.1  1.0 14.5  0.0
    ##   OPC        0.2  0.0  0.1  1.4 67.3 96.3
    ##   RG        13.8 13.1 60.2 10.3  5.3  0.0

``` r
integrated <- AddMetaData(
  object = integrated,
  metadata = predictions$predicted.id.score, 
  col.name = "predicted.score"
)


integrated$predicted.id <- factor(integrated$predicted.id, levels = c("Early IN", "Early EN", "NPC", "RG", "Glioblast", "OPC"))

saveRDS(integrated, "output/RDS/integratedVanBruggen.rds")
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
    ##  [16] globals_0.16.2         spatstat.sparse_3.0-3  RVenn_1.1.0           
    ##  [19] colorspace_2.1-0       xfun_0.38              RCurl_1.98-1.12       
    ##  [22] jsonlite_1.8.4         progressr_0.13.0       spatstat.data_3.0-4   
    ##  [25] survival_3.5-5         zoo_1.8-11             glue_1.6.2            
    ##  [28] polyclip_1.10-4        gtable_0.3.3           zlibbioc_1.44.0       
    ##  [31] XVector_0.38.0         leiden_0.4.3           DelayedArray_0.24.0   
    ##  [34] future.apply_1.10.0    abind_1.4-5            edgeR_3.40.2          
    ##  [37] DBI_1.1.3              spatstat.random_3.2-3  miniUI_0.1.1.1        
    ##  [40] Rcpp_1.0.10            xtable_1.8-4           reticulate_1.43.0     
    ##  [43] dotCall64_1.1-1        htmlwidgets_1.6.2      httr_1.4.5            
    ##  [46] RColorBrewer_1.1-3     ellipsis_0.3.2         ica_1.0-3             
    ##  [49] pkgconfig_2.0.3        rJava_1.0-6            farver_2.1.1          
    ##  [52] uwot_0.1.14            deldir_1.0-6           locfit_1.5-9.7        
    ##  [55] utf8_1.2.3             labeling_0.4.2         tidyselect_1.2.0      
    ##  [58] rlang_1.1.0            reshape2_1.4.4         later_1.3.0           
    ##  [61] munsell_0.5.0          tools_4.2.3            cli_3.6.1             
    ##  [64] generics_0.1.3         ggridges_0.5.4         evaluate_0.20         
    ##  [67] stringr_1.5.0          fastmap_1.1.1          yaml_2.3.7            
    ##  [70] goftest_1.2-3          knitr_1.42             fitdistrplus_1.1-8    
    ##  [73] purrr_1.0.1            RANN_2.6.1             pbapply_1.7-0         
    ##  [76] future_1.32.0          nlme_3.1-162           mime_0.12             
    ##  [79] compiler_4.2.3         rstudioapi_0.14        plotly_4.10.1         
    ##  [82] png_0.1-8              spatstat.utils_3.1-0   tibble_3.2.1          
    ##  [85] stringi_1.7.12         highr_0.10             lattice_0.21-8        
    ##  [88] Matrix_1.6-4           vctrs_0.6.1            pillar_1.9.0          
    ##  [91] lifecycle_1.0.3        spatstat.geom_3.2-9    lmtest_0.9-40         
    ##  [94] RcppAnnoy_0.0.20       cowplot_1.1.1          bitops_1.0-7          
    ##  [97] irlba_2.3.5.1          httpuv_1.6.9           R6_2.5.1              
    ## [100] promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3         
    ## [103] parallelly_1.35.0      codetools_0.2-19       MASS_7.3-58.3         
    ## [106] xlsxjars_0.6.1         rprojroot_2.0.3        withr_2.5.0           
    ## [109] sctransform_0.3.5      GenomeInfoDbData_1.2.9 mgcv_1.8-42           
    ## [112] parallel_4.2.3         grid_4.2.3             rmarkdown_2.21        
    ## [115] Rtsne_0.16             spatstat.explore_3.2-7 shiny_1.7.4
