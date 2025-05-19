Initial Analysis of PSC vs GPC scRNA
================
John Mariani
12/6/2023

``` r
library(Seurat)
library(tidyr)
library(MAST)
library(ggplot2)
library(scPlottingTools)
library(patchwork)
```

``` r
source("Scripts/HelperFunctions.R")
source("Scripts/StyleSettings.R")
```

# Load In Vivo Data for v3.1 samples

``` r
mergedH <- readRDS("output/RDS/mergedH.rds")

table(mergedH$stage)
```

    ## 
    ##                   In Vitro - GPC Stage In Vitro - Pluripotent Stem Cell Stage 
    ##                                  37805                                  17105 
    ##                                In Vivo 
    ##                                   5337

``` r
# Only use v3.1 for this analysis
GPC_PSC <- subset(mergedH, subset = stage != "In Vivo")
GPC_PSC <- subset(GPC_PSC, subset = chemistry == "v3.1")

table(GPC_PSC$stage)
```

    ## 
    ##                   In Vitro - GPC Stage In Vitro - Pluripotent Stem Cell Stage 
    ##                                  33408                                  17105

``` r
DefaultAssay(GPC_PSC) <- "RNA"
GPC_PSC <- NormalizeData(GPC_PSC)
```

## Determine Gene Expresion Fractions

``` r
expressionFractionsPSC <- calcExpressionFractions(GPC_PSC, "stage")


expressionFractionsPSCFilt <- expressionFractionsPSC[expressionFractionsPSC$pct.exp > 10,]
highFractionPSC <- unique(expressionFractionsPSCFilt$features.plot)
```

## Setup MAST ZLM for DE

``` r
# GPC_PSC.sca <- makeSCA(GPC_PSC, highFractionPSC)
# 
# 
# modelMAST <- as.formula(object = "~stage+line+ngeneson+(1|orig.ident)")
# 
# options(mc.cores=8)
# getOption("mc.cores")
# 
# ZLM.GPC_PSC <-MAST::zlm(formula = modelMAST, sca = GPC_PSC.sca, method='glmer',ebayes = F,
#                        strictConvergence = FALSE, parallel = T)
# 
# colnames(ZLM.GPC_PSC@coefC)
# 
# 
# 
# saveRDS(ZLM.GPC_PSC, "output/DE/ZLM.GPC_PSC.rds")
```

## Run DE

``` r
# ZLM.GPC_PSC <- readRDS("output/DE/ZLM.GPC_PSC.rds")
# 
# 
# runLR(ZLM.GPC_PSC, c(0,-1,0,0),
#       contrast0 = c(1,1,1/2,0),
#       contrast1 = c(1,0,1/2,0),
#       fileName = "GPC.vs.PSC")
```

# Add Infinite FC genes that LM assigns NAs

``` r
GPC.vs.PSC <- read.delim("output/DE/GPC.vs.PSC.txt")
#GPC.vs.PSC$logFC <- GPC.vs.PSC$logFC * -1
tempNA <- GPC.vs.PSC[!complete.cases(GPC.vs.PSC),]
tempNA <- tempNA[tempNA$FDR < 0.01,]
tempNA <- tempNA[!is.na(tempNA$FDR),]
expressionFractionsLine <- DotPlot(GPC_PSC, assay = "RNA", features = tempNA$gene, group.by = "stage", split.by = "line")$data
```

    ## Warning: Scaling data with a low number of groups may produce misleading
    ## results

``` r
expressionFractionsLineDF <- pivot_wider(data = expressionFractionsLine, values_from = pct.exp, names_from = id, id_cols = "features.plot")
expressionFractionsLineDF$C27 <- rowMeans(expressionFractionsLineDF[,c(2,3)])
expressionFractionsLineDF$WA09 <- rowMeans(expressionFractionsLineDF[,c(4,5)])

# Keep genes that are not line dependent 
expressionFractionsLineDFfilt <- expressionFractionsLineDF[expressionFractionsLineDF$C27 > .5 & expressionFractionsLineDF$WA09 > .5,]

# Set them to the absolute max log2 fc + .1 for visualization
expressionFractionsLineDFfilt$logFC <- ifelse(expressionFractionsLineDFfilt$`In Vitro - GPC Stage_C27` > expressionFractionsLineDFfilt$`In Vitro - Pluripotent Stem Cell Stage_C27`, max(abs(GPC.vs.PSC$logFC), na.rm = T) + .1, (max(abs(GPC.vs.PSC$logFC), na.rm = T) + .1)*-1)

tempNA <- tempNA[tempNA$gene %in% expressionFractionsLineDFfilt$features.plot,]

tempNA$logFC <- plyr::mapvalues(tempNA$gene, from = expressionFractionsLineDFfilt$features.plot, to = expressionFractionsLineDFfilt$logFC)
tempNA$logFC  <- as.numeric(tempNA$logFC )

GPC.vs.PSC.sig <- GPC.vs.PSC[complete.cases(GPC.vs.PSC) & GPC.vs.PSC$FDR < 0.01 & abs(GPC.vs.PSC$logFC) > 0.25,]
GPC.vs.PSC.sig <- rbind(GPC.vs.PSC.sig, tempNA)

write.table(GPC.vs.PSC.sig, paste0("output/DE/GPC.vs.PSC.sig.txt"), sep = "\t", row.names = F, quote = F)

####
GPC.vs.PSC.sig <- read.delim("output/DE/GPC.vs.PSC.sig.txt")

GPC.vs.PSC.update <- rbind(GPC.vs.PSC[GPC.vs.PSC$gene %not in% GPC.vs.PSC.sig$gene,], GPC.vs.PSC.sig)

write.table(GPC.vs.PSC.update, paste0("output/DE/GPC.vs.PSC.allFC.txt"), sep = "\t", row.names = F, quote = F)
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
    ##  [1] patchwork_1.3.0.9000        scPlottingTools_0.0.0.9000 
    ##  [3] ggplot2_3.4.4               MAST_1.24.1                
    ##  [5] SingleCellExperiment_1.20.1 SummarizedExperiment_1.28.0
    ##  [7] Biobase_2.58.0              GenomicRanges_1.50.2       
    ##  [9] GenomeInfoDb_1.34.9         IRanges_2.32.0             
    ## [11] S4Vectors_0.36.2            BiocGenerics_0.44.0        
    ## [13] MatrixGenerics_1.10.0       matrixStats_0.63.0         
    ## [15] tidyr_1.3.0                 SeuratObject_4.1.3         
    ## [17] Seurat_4.3.0               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-6          
    ##   [4] ellipsis_0.3.2         ggridges_0.5.4         rprojroot_2.0.3       
    ##   [7] XVector_0.38.0         rstudioapi_0.14        spatstat.data_3.0-4   
    ##  [10] leiden_0.4.3           listenv_0.9.0          farver_2.1.1          
    ##  [13] ggrepel_0.9.3          fansi_1.0.4            codetools_0.2-19      
    ##  [16] splines_4.2.3          knitr_1.42             polyclip_1.10-4       
    ##  [19] jsonlite_1.8.4         ica_1.0-3              cluster_2.1.4         
    ##  [22] png_0.1-8              uwot_0.1.14            shiny_1.7.4           
    ##  [25] sctransform_0.3.5      spatstat.sparse_3.0-3  compiler_4.2.3        
    ##  [28] httr_1.4.5             Matrix_1.5-4           fastmap_1.1.1         
    ##  [31] lazyeval_0.2.2         cli_3.6.1              later_1.3.0           
    ##  [34] htmltools_0.5.5        tools_4.2.3            igraph_2.0.3          
    ##  [37] GenomeInfoDbData_1.2.9 gtable_0.3.3           glue_1.6.2            
    ##  [40] RANN_2.6.1             reshape2_1.4.4         dplyr_1.1.1           
    ##  [43] Rcpp_1.0.10            scattermore_0.8        vctrs_0.6.1           
    ##  [46] spatstat.explore_3.2-7 nlme_3.1-162           progressr_0.13.0      
    ##  [49] lmtest_0.9-40          spatstat.random_3.2-3  xfun_0.38             
    ##  [52] stringr_1.5.0          globals_0.16.2         mime_0.12             
    ##  [55] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5.1         
    ##  [58] goftest_1.2-3          future_1.32.0          zlibbioc_1.44.0       
    ##  [61] MASS_7.3-58.3          zoo_1.8-11             scales_1.3.0          
    ##  [64] promises_1.2.0.1       spatstat.utils_3.1-0   parallel_4.2.3        
    ##  [67] RColorBrewer_1.1-3     yaml_2.3.7             reticulate_1.34.0     
    ##  [70] pbapply_1.7-0          gridExtra_2.3          stringi_1.7.12        
    ##  [73] bitops_1.0-7           rlang_1.1.0            pkgconfig_2.0.3       
    ##  [76] evaluate_0.20          lattice_0.21-8         ROCR_1.0-11           
    ##  [79] purrr_1.0.1            tensor_1.5             htmlwidgets_1.6.2     
    ##  [82] cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.35.0     
    ##  [85] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
    ##  [88] R6_2.5.1               generics_0.1.3         DelayedArray_0.24.0   
    ##  [91] DBI_1.1.3              withr_2.5.0            pillar_1.9.0          
    ##  [94] fitdistrplus_1.1-8     RCurl_1.98-1.12        survival_3.5-5        
    ##  [97] abind_1.4-5            sp_1.6-0               tibble_3.2.1          
    ## [100] future.apply_1.10.0    KernSmooth_2.23-20     utf8_1.2.3            
    ## [103] spatstat.geom_3.2-9    plotly_4.10.1          rmarkdown_2.21        
    ## [106] grid_4.2.3             data.table_1.14.8      digest_0.6.31         
    ## [109] xtable_1.8-4           httpuv_1.6.9           munsell_0.5.0         
    ## [112] viridisLite_0.4.1
