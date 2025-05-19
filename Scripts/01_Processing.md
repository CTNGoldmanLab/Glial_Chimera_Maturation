Processing of all scRNA-Seq Data
================
John Mariani
3/6/2023

## Load in Libraries

``` r
library(dplyr)
library(Seurat)
library(devtools)
library(patchwork)
library(SeuratDisk)

options(future.globals.maxSize = 16000 * 1024^2)
```

## Read in Human and Mouse Counts

``` r
sampleList <- list.files("data_for_import/Matrices/")

sampleData <- read.csv("data_for_import/sampleData.csv")

sampleList <- sampleList[sampleList %in% sampleData$sample]

rawH <- sapply(sampleList, function(x) {print(x) ; Read10X(paste0("data_for_import/Matrices/",x,"/human"))})
```

    ## [1] "GOL-2512A5"
    ## [1] "GOL2976A5"
    ## [1] "Sample_32"
    ## [1] "Sample_33"
    ## [1] "Sample_34"
    ## [1] "Sample_C27_CD140_3"
    ## [1] "Sample_C27_DAPI_3"
    ## [1] "Sample_C27_Shi"
    ## [1] "Sample_C27_undiff_1"
    ## [1] "Sample_CD-140"
    ## [1] "Sample_CD140_WA09_5"
    ## [1] "Sample_unsorted"
    ## [1] "Sample_Unsorted_WA09_5"
    ## [1] "Sample_WA09_Shiv"
    ## [1] "Sample_WA09_undiff_1"

``` r
rawM <- sapply(sampleList, function(x) {print(x) ; Read10X(paste0("data_for_import/Matrices/",x,"/mouse"))})
```

    ## [1] "GOL-2512A5"
    ## [1] "GOL2976A5"
    ## [1] "Sample_32"
    ## [1] "Sample_33"
    ## [1] "Sample_34"
    ## [1] "Sample_C27_CD140_3"
    ## [1] "Sample_C27_DAPI_3"
    ## [1] "Sample_C27_Shi"
    ## [1] "Sample_C27_undiff_1"
    ## [1] "Sample_CD-140"
    ## [1] "Sample_CD140_WA09_5"
    ## [1] "Sample_unsorted"
    ## [1] "Sample_Unsorted_WA09_5"
    ## [1] "Sample_WA09_Shiv"
    ## [1] "Sample_WA09_undiff_1"

``` r
# Remove EGFP from tagged C27s
tail(rownames(rawH[[1]]))
```

    ## [1] "ENSG00000275063" "ENSG00000277856" "ENSG00000271254" "ENSG00000268674"
    ## [5] "ENSG00000277475" "eGFP"

``` r
tail(rownames(rawM[[1]]))
```

    ## [1] "ENSMUSG00002074970" "ENSMUSG00002075729" "ENSMUSG00002074899"
    ## [4] "ENSMUSG00002076890" "ENSMUSG00000095742" "eGFP"

``` r
dim(rawH[[1]])
```

    ## [1] 39606   689

``` r
dim(rawM[[1]])
```

    ## [1] 34040   689

``` r
rawH[[1]] <- rawH[[1]][-39606,]
rawH[[2]] <- rawH[[2]][-39606,]

rawM[[1]] <- rawM[[1]][-34040,]
rawM[[2]] <- rawM[[2]][-34040,]


sets <- length(rawH)

#Update Cell names with sample name appended
for(i in 1:sets){ 
  colnames(rawH[[i]]) <- paste0(colnames(rawH[[i]]),"_",sampleList[i])
  colnames(rawM[[i]]) <- paste0(colnames(rawM[[i]]),"_",sampleList[i])
}

head(colnames(rawH[[1]]))
```

    ## [1] "AAACCCAAGGTAATCA_GOL-2512A5" "AAACCCAGTAGGCAAC_GOL-2512A5"
    ## [3] "AAACGAATCAGACCTA_GOL-2512A5" "AAACGCTTCAACCGAT_GOL-2512A5"
    ## [5] "AAACGCTTCAATGCAC_GOL-2512A5" "AAAGAACCATAACTCG_GOL-2512A5"

## Filter for quality and merge datasets

``` r
seurat.objectsH <- sapply(c(1:sets), function(x) CreateSeuratObject(rawH[[x]], project = sampleList[x]))
seurat.objectsH <- sapply(c(1:sets), function(x) PercentageFeatureSet(seurat.objectsH[[x]], pattern = "^MT-", col.name = "percent.mt"))

seurat.objectsM <- sapply(c(1:sets), function(x) CreateSeuratObject(rawM[[x]], project = sampleList[x]))
seurat.objectsM <- sapply(c(1:sets), function(x) PercentageFeatureSet(seurat.objectsM[[x]], pattern = "^mt-", col.name = "percent.mt"))


for (i in 1:sets) {
    seurat.objectsH[[i]] <- subset(x = seurat.objectsH[[i]], subset = nFeature_RNA > 500 & percent.mt < 15)
}


#Subset to only samples that actually have mouse cells
seurat.objectsM <- seurat.objectsM[c(3:5,8,14)]

for (i in 1:length(seurat.objectsM)) {
  seurat.objectsM[[i]] <- subset(x = seurat.objectsM[[i]], subset = nFeature_RNA > 500 & percent.mt < 15)
}


mergedH <- merge(seurat.objectsH[[1]], y = seurat.objectsH[2:length(seurat.objectsH)])
mergedM <- merge(seurat.objectsM[[1]], y = seurat.objectsM[2:length(seurat.objectsM)])
```

## Update metadata

``` r
metaMergedH <- mergedH@meta.data
metaMergedH$cellName <- row.names(metaMergedH)
identical(metaMergedH$cellName, Cells(mergedH))
```

    ## [1] TRUE

``` r
metaMergedH <- merge(metaMergedH, sampleData, by.x = "orig.ident", by.y = "sample")
row.names(metaMergedH) <- metaMergedH$cellName
metaMergedH <- metaMergedH[match(Cells(mergedH), metaMergedH$cellName),]
identical(metaMergedH$cellName, Cells(mergedH))
```

    ## [1] TRUE

``` r
mergedH@meta.data <- metaMergedH


#### Mouse

metaMergedM <- mergedM@meta.data
metaMergedM$cellName <- row.names(metaMergedM)
identical(metaMergedM$cellName, Cells(mergedM))
```

    ## [1] TRUE

``` r
metaMergedM <- merge(metaMergedM, sampleData, by.x = "orig.ident", by.y = "sample")
row.names(metaMergedM) <- metaMergedM$cellName
metaMergedM <- metaMergedM[match(Cells(mergedM), metaMergedM$cellName),]
identical(metaMergedM$cellName, Cells(mergedM))
```

    ## [1] TRUE

``` r
mergedM@meta.data <- metaMergedM

dim(mergedH)
```

    ## [1] 39605 60247

``` r
dim(mergedM)
```

    ## [1] 34039  8398

## Output Data for SCVI integration

``` r
saveRDS(mergedH, "output/RDS/mergedH.rds")
saveRDS(mergedM, "output/RDS/mergedM.rds")

library(SeuratDisk)
SaveH5Seurat(mergedH, filename = "output/H5AD/mergedH.h5Seurat", overwrite = T)
```

    ## Warning: Overwriting previous file output/H5AD/mergedH.h5Seurat

    ## Creating h5Seurat file for version 3.1.5.9900

    ## Adding counts for RNA

    ## Adding data for RNA

    ## No variable features found for RNA

    ## No feature-level metadata found for RNA

``` r
Convert("output/H5AD/mergedH.h5Seurat", dest = "h5ad", overwrite = T)
```

    ## Validating h5Seurat file

    ## Adding data from RNA as X

    ## Adding counts from RNA as raw

    ## Transfering meta.data to obs

``` r
SaveH5Seurat(mergedM, filename = "output/H5AD/mergedM.h5Seurat", overwrite = T)
```

    ## Warning: Overwriting previous file output/H5AD/mergedM.h5Seurat

    ## Creating h5Seurat file for version 3.1.5.9900

    ## Adding counts for RNA

    ## Adding data for RNA

    ## No variable features found for RNA

    ## No feature-level metadata found for RNA

``` r
Convert("output/H5AD/mergedM.h5Seurat", dest = "h5ad", overwrite = T)
```

    ## Validating h5Seurat file

    ## Adding data from RNA as X

    ## Adding counts from RNA as raw

    ## Transfering meta.data to obs

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
    ## [1] SeuratDisk_0.0.0.9020 patchwork_1.3.0.9000  devtools_2.4.5       
    ## [4] usethis_2.1.6         SeuratObject_4.1.3    Seurat_4.3.0         
    ## [7] dplyr_1.1.1          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-6          
    ##   [4] ellipsis_0.3.2         ggridges_0.5.4         rprojroot_2.0.3       
    ##   [7] fs_1.6.1               rstudioapi_0.14        spatstat.data_3.0-4   
    ##  [10] leiden_0.4.3           listenv_0.9.0          farver_2.1.1          
    ##  [13] remotes_2.4.2          bit64_4.0.5            ggrepel_0.9.3         
    ##  [16] fansi_1.0.4            R.methodsS3_1.8.2      codetools_0.2-19      
    ##  [19] splines_4.2.3          cachem_1.0.7           knitr_1.42            
    ##  [22] pkgload_1.3.2          polyclip_1.10-4        jsonlite_1.8.4        
    ##  [25] ica_1.0-3              cluster_2.1.4          R.oo_1.25.0           
    ##  [28] png_0.1-8              uwot_0.1.14            shiny_1.7.4           
    ##  [31] sctransform_0.3.5      spatstat.sparse_3.0-3  compiler_4.2.3        
    ##  [34] httr_1.4.5             Matrix_1.5-4           fastmap_1.1.1         
    ##  [37] lazyeval_0.2.2         cli_3.6.1              later_1.3.0           
    ##  [40] prettyunits_1.1.1      htmltools_0.5.5        tools_4.2.3           
    ##  [43] igraph_2.0.3           gtable_0.3.3           glue_1.6.2            
    ##  [46] RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.10           
    ##  [49] scattermore_0.8        vctrs_0.6.1            spatstat.explore_3.2-7
    ##  [52] nlme_3.1-162           progressr_0.13.0       lmtest_0.9-40         
    ##  [55] spatstat.random_3.2-3  xfun_0.38              stringr_1.5.0         
    ##  [58] ps_1.7.4               globals_0.16.2         mime_0.12             
    ##  [61] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5.1         
    ##  [64] goftest_1.2-3          future_1.32.0          MASS_7.3-58.3         
    ##  [67] zoo_1.8-11             scales_1.3.0           promises_1.2.0.1      
    ##  [70] spatstat.utils_3.1-0   parallel_4.2.3         RColorBrewer_1.1-3    
    ##  [73] yaml_2.3.7             memoise_2.0.1          reticulate_1.34.0     
    ##  [76] pbapply_1.7-0          gridExtra_2.3          ggplot2_3.4.4         
    ##  [79] stringi_1.7.12         pkgbuild_1.4.0         rlang_1.1.0           
    ##  [82] pkgconfig_2.0.3        matrixStats_0.63.0     evaluate_0.20         
    ##  [85] lattice_0.21-8         ROCR_1.0-11            purrr_1.0.1           
    ##  [88] tensor_1.5             htmlwidgets_1.6.2      bit_4.0.5             
    ##  [91] processx_3.8.0         cowplot_1.1.1          tidyselect_1.2.0      
    ##  [94] parallelly_1.35.0      RcppAnnoy_0.0.20       plyr_1.8.8            
    ##  [97] magrittr_2.0.3         R6_2.5.1               profvis_0.3.7         
    ## [100] generics_0.1.3         DBI_1.1.3              withr_2.5.0           
    ## [103] pillar_1.9.0           fitdistrplus_1.1-8     survival_3.5-5        
    ## [106] abind_1.4-5            sp_1.6-0               tibble_3.2.1          
    ## [109] future.apply_1.10.0    hdf5r_1.3.8            crayon_1.5.2          
    ## [112] KernSmooth_2.23-20     utf8_1.2.3             spatstat.geom_3.2-9   
    ## [115] plotly_4.10.1          urlchecker_1.0.1       rmarkdown_2.21        
    ## [118] grid_4.2.3             data.table_1.14.8      callr_3.7.3           
    ## [121] digest_0.6.31          xtable_1.8-4           tidyr_1.3.0           
    ## [124] httpuv_1.6.9           R.utils_2.12.2         munsell_0.5.0         
    ## [127] viridisLite_0.4.1      sessioninfo_1.2.2
