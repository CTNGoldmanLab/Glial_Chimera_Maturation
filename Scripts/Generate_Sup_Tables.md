Make Supplemental Tables
================
2024-07-12

### load library

``` r
library(xlsx)
library(data.table)
library(tidyverse)
```

## Sup Table 1A GPC vs PSC DE

``` r
table1A <- read.delim("output/DE/GPC.vs.PSC.sig.txt")
names(table1A) <- c("Gene", "pVal", "FDR", "Log2FC")
write.xlsx(table1A, file = "output/SupTables/SupTable1.xlsx", sheetName = "GPC vs PSC DE", row.names = F)
```

``` r
# 120 vs ESC K4me3

K4me3_120 <- read.delim("output/CUT&Tag/DiffBind/K4me3_CTd120_diff_peak_annotation_promoter.txt")

K4me3_120 <- K4me3_120[,c(1,2,3,7,8,9,10,20,21)]
names(K4me3_120) <- c("Chr", "Start", "End", "Log2FC", "pVal", "FDR", "Annotation", "Gene", "Description")
K4me3_120 <- K4me3_120[order(K4me3_120$Log2FC, decreasing = T),]

write.xlsx(K4me3_120, file = "output/SupTables/SupTable1.xlsx", sheetName = "H3K4me3 Promoter 120 vs 0", row.names = F, append = T)

# 180 vs ESC K4me3
K4me3_180 <- read.delim("output/CUT&Tag/DiffBind/K4me3_CTd180_diff_peak_annotation_promoter.txt")

K4me3_180 <- K4me3_180[,c(1,2,3,7,8,9,10,20,21)]
names(K4me3_180) <- c("Chr", "Start", "End", "Log2FC", "pVal", "FDR", "Annotation", "Gene", "Description")
K4me3_180 <- K4me3_180[order(K4me3_180$Log2FC, decreasing = T),]

write.xlsx(K4me3_180, file = "output/SupTables/SupTable1.xlsx", sheetName = "H3K4me3 Promoter 180 vs 0", row.names = F, append = T)

#120 vs ESC K27me3
K27me3_120 <- read.delim("output/CUT&Tag/DiffBind/K27me3_CTd120_diff_peak_annotation_promoter.txt")

K27me3_120 <- K27me3_120[,c(1,2,3,7,8,9,10,20,21)]
names(K27me3_120) <- c("Chr", "Start", "End", "Log2FC", "pVal", "FDR", "Annotation", "Gene", "Description")
K27me3_120 <- K27me3_120[order(K27me3_120$Log2FC, decreasing = T),]

write.xlsx(K27me3_120, file = "output/SupTables/SupTable1.xlsx", sheetName = "H3K27me3 Promoter 120 vs 0", row.names = F, append = T)

#180 vs ESC K27me3
K27me3_180 <- read.delim("output/CUT&Tag/DiffBind/K27me3_CTd180_diff_peak_annotation_promoter.txt")

K27me3_180 <- K27me3_180[,c(1,2,3,7,8,9,10,20,21)]
names(K27me3_180) <- c("Chr", "Start", "End", "Log2FC", "pVal", "FDR", "Annotation", "Gene", "Description")
K27me3_180 <- K27me3_180[order(K27me3_180$Log2FC, decreasing = T),]

write.xlsx(K27me3_180, file = "output/SupTables/SupTable1.xlsx", sheetName = "H3K27me3 Promoter 180 vs 0", row.names = F, append = T)
```

## Enhancer SupTables

``` r
geneHancer <- read.delim("data_for_import/CUT&Tag/GeneHancer_AnnotSV_gene_association_scores_v5.18_elite.txt", header = F)

#120 vs ESC K27ac
K27ac_120 <- read.delim("output/CUT&Tag/DiffBind/K27ac_CTd120_diff_peak_annotation_enhancer.txt", header = F)
K27ac_120 <- merge(K27ac_120, geneHancer, by.x = 4, by.y =  1)
K27ac_120 <- K27ac_120[,c(2,3,4,11,12,13,1,14)]
names(K27ac_120) <- c("Chr", "Start", "End", "Log2FC", "pVal", "FDR", "Enhancer", "Gene")

write.xlsx(K27ac_120, file = "output/SupTables/SupTable1.xlsx", sheetName = "H3K27ac Enhancer 120 vs 0", row.names = F, append = T)

K27ac_180 <- read.delim("output/CUT&Tag/DiffBind/K27ac_CTd180_diff_peak_annotation_enhancer.txt", header = F)
K27ac_180 <- merge(K27ac_180, geneHancer, by.x = 4, by.y =  1)
K27ac_180 <- K27ac_180[,c(2,3,4,11,12,13,1,14)]
names(K27ac_180) <- c("Chr", "Start", "End", "Log2FC", "pVal", "FDR", "Enhancer", "Gene")

## Writing this last one keeps crashing. Will drag in manually
#write.xlsx(K27ac_180, file = "output/SupTables/SupTable1.xlsx", sheetName = "H3K27ac Enhancer 180 vs 0", row.names = F, append = T)
write.table(K27ac_180, "output/SupTables/K27ac_180.txt", sep = "\t", quote = F, row.names = F)
```

## In Vitro Supp Tables

``` r
invitroComparisons <- c("GPC1.vs.Rest", "GPC2.vs.Rest", "GPC3.vs.Rest", "GPC4.vs.Rest", "NPC.vs.Rest")

for(i in invitroComparisons){
  temp <- assign(i, read.delim(paste0("output/DE/",i,".txt")))
  print(dim(temp))
  temp$comparison <- i
  temp <- temp[order(temp$logFC, decreasing = T),]
  tempSig <- temp[temp$FDR < 0.01 & abs(temp$logFC) > .25,]
  tempSig <- tempSig[,c(4,2,3)]
  names(tempSig) <- c("Gene", "FDR", "Log2FC")
  write.xlsx(tempSig, file = "output/SupTables/SupTable2.xlsx", sheetName = paste0(i, ".DE"), append = T, row.names = F) 
}
```

    ## [1] 11687     4
    ## [1] 11687     4
    ## [1] 11687     4
    ## [1] 11687     4
    ## [1] 11687     4

## In Vivo Supp Tables

# GPC Enriched Genes

``` r
invivoComparisons1 <- c("imOL.vs.GPC", "imAstrocyte.vs.GPC", "cGPC.vs.GPC", "cGPC.vs.cAPC", "cAPC.vs.GPC", "cGPC.vs.cAPC")


for(i in invivoComparisons1){
  temp <- read.delim(paste0("output/DE/",i,".sig.txt"))
  temp <- temp[,c(4,2,3)]
  names(temp) <- c("Gene", "FDR", "Log2FC")
  temp <- temp[order(temp$Log2FC, decreasing = T),]
  assign(paste0(i,".sig"), temp)
  #write.xlsx(temp, file = "output/SupTables/SupTable3.xlsx", sheetName = paste0(i, ".DE"), append = T, row.names = F) 
}

GPCenriched <- rbind(imOL.vs.GPC.sig[imOL.vs.GPC.sig$Log2FC < 0,], imAstrocyte.vs.GPC.sig[imAstrocyte.vs.GPC.sig$Log2FC < 0,])
GPCenriched <- GPCenriched[GPCenriched$Gene %in% GPCenriched[duplicated(GPCenriched$Gene),]$Gene
,]

imOL.vs.GPC.sig.filt <- imOL.vs.GPC.sig[imOL.vs.GPC.sig$Gene %in% GPCenriched$Gene,]
imAstrocyte.vs.GPC.sig.filt <- imAstrocyte.vs.GPC.sig[imAstrocyte.vs.GPC.sig$Gene %in% GPCenriched$Gene,]

GPCenriched <- merge(imOL.vs.GPC.sig.filt, imAstrocyte.vs.GPC.sig.filt, by.x = 1, by.y = 1)
names(GPCenriched) <- c("Gene", "imOL vs GPC FDR", "imOL vs GPC Log2FC", "imAstrocyte vs GPC FDR", "imAstrocyte vs GPC Log2FC")

GPCenriched <- GPCenriched[order(GPCenriched$`imOL vs GPC Log2FC`, decreasing = F),]

write.xlsx(GPCenriched, file = "output/SupTables/SupTable3.xlsx", sheetName = "GPC Enriched Genes", append = T, row.names = F) 
```

# Cycling Enriched Genes

``` r
cyclingEnriched <- rbind(cGPC.vs.GPC.sig[cGPC.vs.GPC.sig$Log2FC > 0,], cAPC.vs.GPC.sig[cAPC.vs.GPC.sig$Log2FC > 0,])
cyclingEnriched <- cyclingEnriched[cyclingEnriched$Gene %in% cyclingEnriched[duplicated(cyclingEnriched$Gene),]$Gene
,]

cGPC.vs.GPC.sig.filt <- cGPC.vs.GPC.sig[cGPC.vs.GPC.sig$Gene %in% cyclingEnriched$Gene,]
cAPC.vs.GPC.sig.filt <- cAPC.vs.GPC.sig[cAPC.vs.GPC.sig$Gene %in% cyclingEnriched$Gene,]

cyclingEnriched <- merge(cGPC.vs.GPC.sig.filt, cAPC.vs.GPC.sig.filt, by.x = 1, by.y = 1)
names(cyclingEnriched) <- c("Gene", "cGPC vs GPC FDR", "cGPC vs GPC Log2FC", "cAPC vs GPC FDR", "cAPC vs GPC Log2FC")


cyclingEnriched <- cyclingEnriched[order(cyclingEnriched$`cGPC vs GPC Log2FC`, decreasing = T),]

write.xlsx(cyclingEnriched, file = "output/SupTables/SupTable3.xlsx", sheetName = "Cycling Enriched Genes", append = T, row.names = F) 
```

# cGPC vs cAPC Genes

``` r
write.xlsx(cGPC.vs.cAPC.sig, file = "output/SupTables/SupTable3.xlsx", sheetName = "cGPC vs cAPC DE", append = T, row.names = F) 
```

## Invitro Invivo Supp Tables

# In Vivo vs In Vitro GPC DE

``` r
Invivo.vs.Invitro.GPC4.sig <- read.delim("output/DE/Invivo.vs.Invitro.GPC4.sig.txt")
Invivo.vs.Invitro.GPC4.sig <- Invivo.vs.Invitro.GPC4.sig[,c(4,2,3)]
names(Invivo.vs.Invitro.GPC4.sig) <- c("Gene", "FDR", "Log2FC")
Invivo.vs.Invitro.GPC4.sig <- Invivo.vs.Invitro.GPC4.sig[order(Invivo.vs.Invitro.GPC4.sig$Log2FC, decreasing = T),]

write.xlsx(Invivo.vs.Invitro.GPC4.sig, file = "output/SupTables/SupTable4.xlsx", sheetName = "In Vivo vs In Vitro GPC4 DE", append = T, row.names = F) 
```

# In Vivo vs In Vitro GPC SCENIC

``` r
Invivo.vs.Invitro.GPC4.SCENIC.sig <- read.csv("output/DE/Invivo.vs.Invitro.SCENICE.sig.csv")
names(Invivo.vs.Invitro.GPC4.SCENIC.sig) <- c("Regulon", "FDR", "Log2FC")

write.xlsx(Invivo.vs.Invitro.GPC4.SCENIC.sig, file = "output/SupTables/SupTable4.xlsx", sheetName = "In Vivo vs In Vitro GPC4 SCENIC", append = T, row.names = F) 
```

# In Vivo and In Vitro SCENIC Targets

``` r
invitroTF.supp <- read.delim("output/Networks/Invitro_Invivo/invitroTF.supp.txt")
invivoTF.supp <- read.delim("output/Networks/Invitro_Invivo/invivoTF.supp.txt")

write.xlsx(invitroTF.supp, file = "output/SupTables/SupTable4.xlsx", sheetName = "In Vitro GPC Regulon Targets", append = T, row.names = F) 
write.xlsx(invivoTF.supp, file = "output/SupTables/SupTable4.xlsx", sheetName = "In Vivo GPC Regulon Targets", append = T, row.names = F) 
```

# In Vivo vs In Vitro GPC IPA terms

``` r
IPA <- read.delim("output/IPA/Invivo.vs.Invitro.GPC4.IPA.sig.txt")

write.xlsx(IPA, file = "output/SupTables/SupTable4.xlsx", sheetName = "In Vivo vs In Vitro IPA", append = T, row.names = F) 
```

## Mouse and Nichenet Supp Tables

# Mouse Markers

``` r
mouseMarkers <- read.table("output/DE/MouseMarkers.txt", header = T)
mouseMarkers <- mouseMarkers[mouseMarkers$p_val_adj < 0.01, ]
mouseMarkers <- mouseMarkers[order(mouseMarkers$cluster, mouseMarkers$avg_log2FC, decreasing = T),]
mouseMarkers <- mouseMarkers[,c(7,5,2,3,4,6)]
names(mouseMarkers) <- c("Gene", "Adj_pVal", "Log2FC", "Pct.1", "Pct.2", "Cluster")

## Also too big, need to use a different package in the future
#write.xlsx(mouseMarkers, file = "output/SupTables/SupTable5.xlsx", sheetName = "Mouse Gene Markers", append = T, row.names = F) 

write.table(mouseMarkers, "output/NicheNet/mouseMarkersForSuppTable.txt", sep = "\t", row.names = F, quote = F)
```

# Ligand Activities and Targets

``` r
ligandActivities <- read.delim("output/NicheNet/ligand_activities_targets.txt")
names(ligandActivities) <- c("Ligand", "AUPR Corrected", "Ligand Activity", "Target", "Ligand Target Weight", "Receiver", "Ligand Activity Normalized", "Scaled Ligand Activity Normalized", "Scaled Ligand Activity")

ligandActivities <- ligandActivities[order(ligandActivities$`Ligand Activity`, ligandActivities$`Ligand Target Weight`, decreasing = T),]

# This is too big... Going to drop it in manually
# write.xlsx(ligandActivities, file = "output/SupTables/SupTable5.xlsx", sheetName = "Ligand Activity and Targets", append = T, row.names = F) 

write.table(ligandActivities, "output/NicheNet/ligandActivitiesForSuppTable.txt", sep = "\t", row.names = F, quote = F)
```

# Ligand Expression

``` r
ligands <- unique(ligandActivities$Ligand)

ligandExpr <- read.delim("output/NicheNet/plotting_tbl.txt")

ligandExpr <- ligandExpr[ligandExpr$gene %in% ligands,]
ligandExpr <- ligandExpr[,c(1:5)]

ligandExpr <- ligandExpr[order(ligandExpr$gene, ligandExpr$celltype),]
names(ligandExpr) <- c("Cell type", "Gene", "Expression", "Scaled Expression", "Fraction Expressed")

write.xlsx(ligandExpr, file = "output/SupTables/SupTable5.xlsx", sheetName = "Ligand Expression", append = T, row.names = F) 
```

# Receptor Expression

``` r
receptorSup <- read.delim("output/NicheNet/receptorSuppTable.txt")
names(receptorSup) <- c("Receiver", "Ligand", "Receptor", "Receptor Expr", "Receptor Fraction")

write.xlsx(receptorSup, file = "output/SupTables/SupTable5.xlsx", sheetName = "Receptor Expression", append = T, row.names = F) 
```

## Invivo Oligo Supp Tables

``` r
# Gene Expression
imOL.vs.GPC.sig <- read.delim("output/DE/imOL.vs.GPC.sig.txt")
imOL.vs.GPC.sig <- imOL.vs.GPC.sig[,c(4,2,3)]
names(imOL.vs.GPC.sig) <- c("Gene", "FDR", "Log2FC")
imOL.vs.GPC.sig <- imOL.vs.GPC.sig[order(imOL.vs.GPC.sig$Log2FC, decreasing = T),]

maOL.vs.imOL.sig <- read.delim("output/DE/maOL.vs.imOL.sig.txt")
maOL.vs.imOL.sig <- maOL.vs.imOL.sig[,c(4,2,3)]
names(maOL.vs.imOL.sig) <- c("Gene", "FDR", "Log2FC")
maOL.vs.imOL.sig <- maOL.vs.imOL.sig[order(maOL.vs.imOL.sig$Log2FC, decreasing = T),]

# SCENIC AUC
imOL.vs.GPC.SCENIC.sig <- read.delim("output/DE/imOL.vs.GPC.SCENIC.sig.txt")
imOL.vs.GPC.SCENIC.sig <- imOL.vs.GPC.SCENIC.sig[,c(6,5,2)]
names(imOL.vs.GPC.SCENIC.sig) <- c("TF", "FDR", "Log2FC")
imOL.vs.GPC.SCENIC.sig <- imOL.vs.GPC.SCENIC.sig[order(imOL.vs.GPC.SCENIC.sig$Log2FC, decreasing = T),]
  
maOL.vs.imOL.SCENIC.sig <- read.delim("output/DE/maOL.vs.imOL.SCENIC.sig.txt")
maOL.vs.imOL.SCENIC.sig <- maOL.vs.imOL.SCENIC.sig[,c(6,5,2)]
names(maOL.vs.imOL.SCENIC.sig) <- c("TF", "FDR", "Log2FC")
maOL.vs.imOL.SCENIC.sig <- maOL.vs.imOL.SCENIC.sig[order(maOL.vs.imOL.SCENIC.sig$Log2FC, decreasing = T),]

#SCENIC Targets
imOL.vs.GPC.TF <- read.delim("output/Networks/Invivo/imOL.vs.GPC.tf.txt")
imOL.vs.GPC.TF <- imOL.vs.GPC.TF[,-9]

maOL.vs.imOL.TF <- read.delim("output/Networks/Invivo/maOL.vs.imOL.tf.txt")
maOL.vs.imOL.TF <- maOL.vs.imOL.TF[,-9]

#Write to SupTable
write.xlsx(imOL.vs.GPC.sig, file = "output/SupTables/SupTable6.xlsx", sheetName = "imOL vs GPC DE", append = T, row.names = F) 
write.xlsx(imOL.vs.GPC.SCENIC.sig, file = "output/SupTables/SupTable6.xlsx", sheetName = "imOL vs GPC SCENIC", append = T, row.names = F) 
write.xlsx(imOL.vs.GPC.TF, file = "output/SupTables/SupTable6.xlsx", sheetName = "imOL vs GPC Regulon Targets", append = T, row.names = F) 


write.xlsx(maOL.vs.imOL.sig, file = "output/SupTables/SupTable6.xlsx", sheetName = "maOL vs imOL DE", append = T, row.names = F) 
write.xlsx(maOL.vs.imOL.SCENIC.sig, file = "output/SupTables/SupTable6.xlsx", sheetName = "maOL vs imOL SCENIC", append = T, row.names = F) 
write.xlsx(maOL.vs.imOL.TF, file = "output/SupTables/SupTable6.xlsx", sheetName = "maOL vs imOL Regulon Targets", append = T, row.names = F) 
```

## Invivo Astro Supp Tables

``` r
# Gene Expression
imAstro.vs.GPC.sig <- read.delim("output/DE/imAstrocyte.vs.GPC.sig.txt")
imAstro.vs.GPC.sig <- imAstro.vs.GPC.sig[,c(4,2,3)]
names(imAstro.vs.GPC.sig) <- c("Gene", "FDR", "Log2FC")
imAstro.vs.GPC.sig <- imAstro.vs.GPC.sig[order(imAstro.vs.GPC.sig$Log2FC, decreasing = T),]

Astro.vs.imAstro.sig <- read.delim("output/DE/Astrocyte.vs.imAstrocyte.sig.txt")
Astro.vs.imAstro.sig <- Astro.vs.imAstro.sig[,c(4,2,3)]
names(Astro.vs.imAstro.sig) <- c("Gene", "FDR", "Log2FC")
Astro.vs.imAstro.sig <- Astro.vs.imAstro.sig[order(Astro.vs.imAstro.sig$Log2FC, decreasing = T),]


# SCENIC AUC
imAstro.vs.GPC.SCENIC.sig <- read.delim("output/DE/imAstro.vs.GPC.SCENIC.sig.txt")
imAstro.vs.GPC.SCENIC.sig <- imAstro.vs.GPC.SCENIC.sig[,c(6,5,2)]
names(imAstro.vs.GPC.SCENIC.sig) <- c("TF", "FDR", "Log2FC")
imAstro.vs.GPC.SCENIC.sig <- imAstro.vs.GPC.SCENIC.sig[order(imAstro.vs.GPC.SCENIC.sig$Log2FC, decreasing = T),]

Astro.vs.imAstro.SCENIC.sig <- read.delim("output/DE/Astro.vs.imAstro.SCENIC.sig.txt")
Astro.vs.imAstro.SCENIC.sig <- Astro.vs.imAstro.SCENIC.sig[,c(6,5,2)]
names(Astro.vs.imAstro.SCENIC.sig) <- c("TF", "FDR", "Log2FC")
Astro.vs.imAstro.SCENIC.sig <- Astro.vs.imAstro.SCENIC.sig[order(Astro.vs.imAstro.SCENIC.sig$Log2FC, decreasing = T),]

#SCENIC Targets
imAstro.vs.GPC.TF <- read.delim("output/Networks/Invivo/imAstro.vs.GPC.tf.txt")
imAstro.vs.GPC.TF <- imAstro.vs.GPC.TF[,-9]

Astro.vs.imAstro.TF <- read.delim("output/Networks/Invivo/Astro.vs.imAstro.tf.txt")
Astro.vs.imAstro.TF <- Astro.vs.imAstro.TF[,-9]

#Write to SupTable
write.xlsx(imAstro.vs.GPC.sig, file = "output/SupTables/SupTable7.xlsx", sheetName = "imAstro vs GPC DE", append = T, row.names = F) 
write.xlsx(imAstro.vs.GPC.SCENIC.sig, file = "output/SupTables/SupTable7.xlsx", sheetName = "imAstro vs GPC SCENIC", append = T, row.names = F) 
write.xlsx(imAstro.vs.GPC.TF, file = "output/SupTables/SupTable7.xlsx", sheetName = "imAstro vs GPC Regulon Targets", append = T, row.names = F) 


write.xlsx(Astro.vs.imAstro.sig, file = "output/SupTables/SupTable7.xlsx", sheetName = "Astro vs imAstro DE", append = T, row.names = F) 
write.xlsx(Astro.vs.imAstro.SCENIC.sig, file = "output/SupTables/SupTable7.xlsx", sheetName = "Astro vs imAstro SCENIC", append = T, row.names = F) 
write.xlsx(Astro.vs.imAstro.TF, file = "output/SupTables/SupTable7.xlsx", sheetName = "Astro vs imAstro Regulon Targets", append = T, row.names = F) 
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
    ##  [1] lubridate_1.9.2   forcats_1.0.0     stringr_1.5.0     dplyr_1.1.1      
    ##  [5] purrr_1.0.1       readr_2.1.4       tidyr_1.3.0       tibble_3.2.1     
    ##  [9] ggplot2_3.4.4     tidyverse_2.0.0   data.table_1.14.8 xlsx_0.6.5       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.2.3   pillar_1.9.0     tools_4.2.3      digest_0.6.31   
    ##  [5] timechange_0.2.0 evaluate_0.20    lifecycle_1.0.3  gtable_0.3.3    
    ##  [9] pkgconfig_2.0.3  rlang_1.1.0      cli_3.6.1        rstudioapi_0.14 
    ## [13] yaml_2.3.7       xfun_0.38        fastmap_1.1.1    rJava_1.0-6     
    ## [17] withr_2.5.0      knitr_1.42       generics_0.1.3   xlsxjars_0.6.1  
    ## [21] vctrs_0.6.1      hms_1.1.3        rprojroot_2.0.3  grid_4.2.3      
    ## [25] tidyselect_1.2.0 glue_1.6.2       R6_2.5.1         fansi_1.0.4     
    ## [29] rmarkdown_2.21   tzdb_0.3.0       magrittr_2.0.3   scales_1.3.0    
    ## [33] htmltools_0.5.5  colorspace_2.1-0 utf8_1.2.3       stringi_1.7.12  
    ## [37] munsell_0.5.0
