Analysis of Cultured hGPCs
================
John Mariani
10/26/2025

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
library(ggridges)

`%not in%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))

axisTitleSize <- 24
axisTextSize <- 18
labelFont = 18
titleFont = 26
tagSize = 26
fig.pt.size = 2

# Because I have a newer function for this but didn't want to remake all the plots this time around to scale
theme_manuscriptOG <-  theme(axis.text = element_text(size = axisTextSize), 
        axis.title = element_text(size = axisTitleSize), 
        title = element_text(size = titleFont), 
        legend.title = element_text(size = titleFont),
        legend.text = element_text(size = axisTitleSize),
        plot.tag = element_text(size = tagSize))

manuscriptPalette <- c("In Vivo" = "red2", 
                       "In Vitro - GPC Stage" = "#2E30FF",
                       "NPC" = "magenta",
                       "GPC1" = "forestgreen",
                       "GPC2" = "darkorange",
                       "GPC3" = "firebrick2",
                       "GPC4" = "turquoise",
                       "Astrocyte" = "dodgerblue2",
                       "imOL" = "gold",
                       "maOL" = "darkorchid4")

source("Scripts/HelperFunctions.R")
source("Scripts/StyleSettings.R")
```

## Make invitro only subset

``` r
invitroInvivo <- readRDS("output/RDS/invitroInvivo.rds")

invitro <- subset(invitroInvivo, subset = stage == "In Vitro - GPC Stage")

DimPlot(invitro, split.by = "cellType")
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
invitroDE <- subset(invitro, subset = cellType %in% c("NPC", "GPC1", "GPC2", "GPC3", "GPC4"))

DimPlot(invitroDE, split.by = "cellType")
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
invitroDE <- NormalizeData(invitroDE)
```

## Determine Gene Expresion Fractions and keep those that are present in 10% of any of the in vitro clusters

``` r
expressionFractions <- DotPlot(invitroDE, assay = "RNA", features = row.names(invitro))$data
names(expressionFractions)
```

    ## [1] "avg.exp"        "pct.exp"        "features.plot"  "id"            
    ## [5] "avg.exp.scaled"

``` r
expressionFractionsFilt <- expressionFractions[expressionFractions$pct.exp > 10,]
highFraction <- unique(expressionFractionsFilt$features.plot)

expressionFractionsDF <- pivot_wider(data = expressionFractions, values_from = pct.exp, names_from = id, id_cols = "features.plot")
```

``` r
# invitro.sca <- makeSCA(invitroDE, highFraction) 
# dim(invitro.sca)
```

## Fit mixed model with fixed cell line, 10x chemistry, and feature complexity effects and a capture random effect

``` r
# options(mc.cores=8)
# getOption("mc.cores")
# 
# modelMAST <- as.formula(object = "~cellType+line+chemistry+ngeneson+(1|orig.ident)")
# 
# ZLM.invitro <-MAST::zlm(formula = modelMAST, sca = invitro.sca, method='glmer',ebayes = F,
#                                             strictConvergence = FALSE, parallel = T)
# 
# #saveRDS(ZLM.invitro, "output/DE/ZLM.invitro.rds")
# ZLM.invitro <- readRDS("output/DE/ZLM.invitro.rds")
# colnames(ZLM.invitro@coefC)
```

## DE for enrichment in each cell state

``` r
# options(mc.cores=8)
# getOption("mc.cores")
# 
# #NPC vs rest 
# temp <- colnames(ZLM.invitro@coefC)
# names(temp) <- temp
# temp[] <- c(0,-1/4,-1/4,-1/4,-1/4,0,0,0)
# temp
# 
# runLR(ZLM.invitro, lrContrast = c(0,-1/4,-1/4,-1/4,-1/4,0,0,0),
#       contrast1 = c(1,0,0,0,0,1/2,1/2,0), 
#       contrast0 = c(1,1/4,1/4,1/4,1/4,1/2,1/2,0), FDR = 0.01, logFC = .25, fileName = "NPC.vs.Rest")
# 
# # GPC4 vs rest
# temp[] <- c(0, -1/4, -1/4, -1/4, 1,0,0,0)
# temp
# 
# runLR(ZLM.invitro, lrContrast = c(0, -1/4, -1/4, -1/4, 1,0,0,0),
#       contrast1 = c(1,0,0,0,1,1/2,1/2,0), 
#       contrast0 = c(1,1/4,1/4,1/4,0,1/2,1/2,0), FDR = 0.01, logFC = .25, fileName = "GPC4.vs.Rest")
# 
# 
# # GPC1 vs rest
# temp[] <- c(0, 1, -1/4, -1/4, -1/4,0,0,0)
# temp
# 
# runLR(ZLM.invitro, lrContrast = c(0, 1, -1/4, -1/4, -1/4,0,0,0),
#       contrast1 = c(1,1,0,0,0,1/2,1/2,0), 
#       contrast0 = c(1, 0, 1/4, 1/4, 1/4,1/2,1/2,0), FDR = 0.01, logFC = .25, fileName = "GPC1.vs.Rest")
#  
# # GPC2 vs rest
# temp[] <- c(0, -1/4, 1, -1/4, -1/4,0,0,0)
# temp
# 
# runLR(ZLM.invitro, lrContrast = c(0, -1/4, 1, -1/4, -1/4,0,0,0),
#       contrast1 = c(1,0,1,0,0,1/2,1/2,0), 
#       contrast0 = c(1, 1/4, 0, 1/4, 1/4,1/2,1/2,0), FDR = 0.01, logFC = .25, fileName = "GPC2.vs.Rest")
#  
# #GPC3 vs rest
# temp[] <- c(0, -1/4, -1/4, 1, -1/4,0,0,0)
# temp
# 
# runLR(ZLM.invitro, lrContrast = c(0, -1/4, -1/4, 1, -1/4,0,0,0),
#       contrast1 = c(1,0,0,1,0,1/2,1/2,0), 
#       contrast0 = c(1, 1/4, 1/4, 0, 1/4,1/2,1/2,0), FDR = 0.01, logFC = .25, fileName = "GPC3.vs.Rest")
```

## Reload DE to save time

``` r
invitroComparisons <- c("GPC4.vs.Rest", "NPC.vs.Rest", "GPC1.vs.Rest", "GPC2.vs.Rest", "GPC3.vs.Rest")

for(i in invitroComparisons){
  temp <- assign(i, read.delim(paste0("output/DE/",i,".txt")))
  print(dim(temp))
  temp$comparison <- i
  temp <- temp[order(temp$logFC, decreasing = T),]
  assign(i, temp)
  temp <- temp[complete.cases(temp),]
  assign(paste0(i,".sig"), temp[temp$FDR < 0.01 & abs(temp$logFC) > .25,])
}
```

    ## [1] 11687     4
    ## [1] 11687     4
    ## [1] 11687     4
    ## [1] 11687     4
    ## [1] 11687     4

## Expression HM

``` r
# How many top genes to show in heatmap
howMany <- 11

topHowMany <- GPC4.vs.Rest.sig$gene[1:howMany]
topHowMany <- c(topHowMany, GPC3.vs.Rest.sig[GPC3.vs.Rest.sig$gene %not in% topHowMany,]$gene[1:howMany])
topHowMany <- c(topHowMany, GPC2.vs.Rest.sig[GPC2.vs.Rest.sig$gene %not in% topHowMany,]$gene[1:howMany])
topHowMany <- c(topHowMany, GPC1.vs.Rest.sig[GPC1.vs.Rest.sig$gene %not in% topHowMany & !grepl(pattern = "ENSG", GPC1.vs.Rest.sig$gene),]$gene[1:howMany])
topHowMany <- c(topHowMany, NPC.vs.Rest.sig[NPC.vs.Rest.sig$gene %not in% topHowMany & !grepl(pattern = "ENSG", NPC.vs.Rest.sig$gene),]$gene[1:howMany])


allDE <- rbindlist(list(GPC4.vs.Rest, GPC3.vs.Rest, GPC2.vs.Rest, GPC1.vs.Rest, NPC.vs.Rest))
allDE <- allDE[allDE$gene %in% topHowMany,]
allDE$comparison <- gsub(x = allDE$comparison, pattern = ".vs.Rest", replacement = "")
allDE$comparison <- factor(allDE$comparison, levels = rev(c("NPC", "GPC1", "GPC2", "GPC3","GPC4")))
allDE$gene <- factor(allDE$gene, levels = rev(topHowMany))

allDE$sig <- symnum(allDE$FDR, cutpoints = c(0, 0.0001,
    0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*"," "))

allDE$sig <- ifelse(abs(allDE$logFC) < .25, " ", allDE$sig)

expressionHMFig <- ggplot(allDE, aes(y = comparison, x = gene, fill = logFC)) + geom_tile(colour = "black") + scale_fill_gradient2(low = "dodgerblue2", mid = "white", high = "red2", midpoint = 0) + scale_y_discrete(expand = c(0,0)) +theme_bw() + theme_manuscriptOG + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.key.height = unit(1, "inch")) + geom_text(aes(label = sig), angle = 90, vjust = 1, size = 12) + guides(fill=guide_colorbar(title="Log2FC Enrichment", title.position = "left", title.theme = element_text(angle = 90, hjust = .5)))

expressionHMFig
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## Integrated Dim Plot

``` r
dim(invitroInvivo)
```

    ## [1] 39605 43142

``` r
integratedDimFig <- DimPlotCustom(invitroInvivo, group.by = "cellType", label = T, pt.size = fig.pt.size) + theme_bw() + theme_manuscriptOG + theme(legend.position = "bottom") + labs(title = "Integrated - 43,142 Cells") + scale_fill_manual(values = manuscriptPalette, aes(legend.title = "CellType")) + guides(fill = guide_legend(override.aes = list(size = 5)))

integratedDimFig
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Dim Plot Split by stage

``` r
stageDimFig <- DimPlotCustom(invitroInvivo, split.by = "stage", group.by = "cellType", ncol = 1, pt.size = fig.pt.size) & theme_bw() & theme_manuscriptOG & NoLegend() & theme(axis.title.y = element_blank())  & scale_fill_manual(values = manuscriptPalette)

table(invitroInvivo$stage)
```

    ## 
    ## In Vitro - GPC Stage              In Vivo 
    ##                37805                 5337

``` r
stageDimFig[[1]] <- stageDimFig[[1]] + labs(title = "In Vivo - 5,337 Cells") 
stageDimFig[[2]] <- stageDimFig[[2]] + xlab("UMAP 1") + labs(title = "In Vitro - 37,805 Cells")
stageDimFig
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
stageDimFig[[2]] | stageDimFig[[1]]
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

## Stacked Plot Cell Type by Stage

``` r
StackedCellType <- as.data.frame(table(invitroInvivo$stage, invitroInvivo$cellType))

unique(StackedCellType$Var1)
```

    ## [1] In Vitro - GPC Stage In Vivo             
    ## Levels: In Vitro - GPC Stage In Vivo

``` r
stackedPlotInvitroInvivo <- ggplot(StackedCellType, aes(fill = Var2, y = Freq, x = Var1))+
  geom_bar(position = "fill", stat = "identity")+
  guides(fill = guide_legend(override.aes = list(size = .6))) + theme_bw() + theme_manuscriptOG + scale_y_continuous(expand = c(0,0))  + labs(y = "Fraction Celltype", x = "Stage", fill = "Cell Type") + theme(legend.key.size = unit(.6, "lines")) + xlab(element_blank()) + NoLegend() & scale_fill_manual(values = manuscriptPalette)

stackedPlotInvitroInvivo 
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Violin Plots

``` r
gpcViolins <- (VlnPlot(invitroDE, c("PDGFRA"), pt.size = 0) & coord_flip() & scale_fill_manual(values = manuscriptPalette)& scale_y_continuous(expand = c(0,0), limits = c(0,3.75)) & theme_bw() & theme_manuscriptOG & NoLegend() +  theme(axis.title.y = element_blank(), axis.title.x = element_blank())) /
(VlnPlot(invitroDE, c("NEUROD1"), pt.size = .01) & coord_flip() & scale_fill_manual(values = manuscriptPalette) & scale_y_continuous(expand = c(0,0), limits = c(0,3.75)) & theme_bw() & theme_manuscriptOG & NoLegend() + theme(axis.title.y = element_blank()))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

## Load in Palantir Data and update seurat metadata

``` r
palantirMetadata <- read.csv("output/Palantir/palantirInvitroInvivoMetadata.csv", row.names = 1)
palantirMetadata <- palantirMetadata[,c(14:21)]

invitroInvivoMeta <- invitroInvivo@meta.data
invitroInvivoMeta <- merge(invitroInvivoMeta, palantirMetadata, by.x = 0, by.y = 0)
row.names(invitroInvivoMeta) <- invitroInvivoMeta$Row.names
invitroInvivoMeta$Row.names <- NULL

identical(row.names(invitroInvivoMeta), Cells(invitroInvivo))
```

    ## [1] FALSE

``` r
invitroInvivoMeta <- invitroInvivoMeta[match(Cells(invitroInvivo), row.names(invitroInvivoMeta)),]

identical(row.names(invitroInvivoMeta), Cells(invitroInvivo))
```

    ## [1] TRUE

``` r
cellLevels <-c("Unselected", "NPC", "GPC1", "GPC2", "GPC3", "GPC4", "imOL", "maOL", "Astrocyte")


invitroInvivoMeta$NPC_branch <- factor(ifelse(invitroInvivoMeta$NPC_branch == T, as.character(invitroInvivoMeta$cellType), "Unselected"), levels = cellLevels)
invitroInvivoMeta$Oligo_branch <- factor(ifelse(invitroInvivoMeta$Oligo_branch == T, as.character(invitroInvivoMeta$cellType), "Unselected"), levels = cellLevels)
invitroInvivoMeta$Astro_branch <- factor(ifelse(invitroInvivoMeta$Astro_branch == T, as.character(invitroInvivoMeta$cellType), "Unselected"), levels = cellLevels)


invitroInvivoMeta$NPC_branch <- factor(invitroInvivoMeta$NPC_branch , levels = cellLevels)
invitroInvivoMeta$Oligo_branch <- factor(invitroInvivoMeta$Oligo_branch , levels = cellLevels)
invitroInvivoMeta$Astro_branch <- factor(invitroInvivoMeta$Astro_branch , levels = cellLevels)


levels(invitroInvivoMeta$NPC_branch)
```

    ## [1] "Unselected" "NPC"        "GPC1"       "GPC2"       "GPC3"      
    ## [6] "GPC4"       "imOL"       "maOL"       "Astrocyte"

``` r
invitroInvivo@meta.data <- invitroInvivoMeta
```

``` r
xmin <- min(invitroInvivo@reductions$umap@cell.embeddings[,1])
xmax <- max(invitroInvivo@reductions$umap@cell.embeddings[,1])

ymin <- min(invitroInvivo@reductions$umap@cell.embeddings[,2])
ymax <- max(invitroInvivo@reductions$umap@cell.embeddings[,2])


third <- (DimPlotCustom(subset(invitroInvivo, subset = Oligo_branch != "Unselected"), group.by = "Oligo_branch", label = T, pt.size = fig.pt.size) & scale_fill_manual(values = manuscriptPalette) & ggtitle("Oligodendrocyte Branch") & labs(tag = "G") | 
DimPlotCustom(subset(invitroInvivo, subset = Astro_branch != "Unselected"), group.by = "Astro_branch", label = T, pt.size = fig.pt.size) & scale_fill_manual(values = manuscriptPalette) & ggtitle("Astrocyte Branch") & labs(tag = "I") |
DimPlotCustom(subset(invitroInvivo, subset = NPC_branch != "Unselected"), group.by = "NPC_branch", label = T, pt.size = fig.pt.size) & scale_fill_manual(values = manuscriptPalette) & ggtitle("NPC Branch") & labs(tag = "K")) & theme_bw() & theme_manuscriptOG & NoLegend() & xlim(xmin, xmax) & ylim(ymin, ymax)
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
third[[2]] <- third[[2]] & theme(axis.title.y = element_blank())

third[[3]] <- third[[3]] & theme(axis.title.y = element_blank())

third
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
levels(invitroInvivoMeta$Oligo_branch)
```

    ## [1] "Unselected" "NPC"        "GPC1"       "GPC2"       "GPC3"      
    ## [6] "GPC4"       "imOL"       "maOL"       "Astrocyte"

``` r
invitroInvivoMeta$Oligo_branch <- factor(invitroInvivoMeta$Oligo_branch, levels = cellLevels)

bottom <- (ggplot(invitroInvivoMeta, aes(x = invitroInvivoMeta$palantir_pseudotime, y = invitroInvivoMeta$Oligo_bp, colour = invitroInvivoMeta$Oligo_branch, group = invitroInvivoMeta$Oligo_branch)) + geom_point(size = 3) + scale_color_manual(values = manuscriptPalette) & labs(tag = "H") |
ggplot(invitroInvivoMeta, aes(x = invitroInvivoMeta$palantir_pseudotime, y = invitroInvivoMeta$Astro_bp, colour = invitroInvivoMeta$Astro_branch)) + geom_point(size = 3) + scale_color_manual(values = manuscriptPalette) & labs(tag = "J") | 
ggplot(invitroInvivoMeta, aes(x = invitroInvivoMeta$palantir_pseudotime, y = invitroInvivoMeta$NPC_bp, colour = invitroInvivoMeta$NPC_branch)) + geom_point(size = 3) + scale_color_manual(values = manuscriptPalette) & labs(tag = "L")) & theme_bw() & theme_manuscriptOG & NoLegend() & labs(x = "Pseudotime", y = "Branch Probability")

bottom[[2]] <- bottom[[2]] & theme(axis.title.y = element_blank())

bottom[[3]] <- bottom[[3]] & theme(axis.title.y = element_blank())

bottom
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Piece together main figure

``` r
integratedDimFig <- integratedDimFig + labs(tag = "A")
stageDimFig[[2]] <- stageDimFig[[2]] + labs(tag = "B")
stageDimFig[[1]] <- stageDimFig[[1]] + labs(tag = "C")
stackedPlotInvitroInvivo <- stackedPlotInvitroInvivo + labs(tag = "D")
gpcViolins[[1]] <- gpcViolins[[1]] + labs(tag = "E")


top <- (integratedDimFig | stageDimFig[[2]] | stageDimFig[[1]] | (stackedPlotInvitroInvivo / gpcViolins))  + plot_layout(widths = c(1,1,1,.45))

second <- (expressionHMFig + labs(tag = "F") + plot_spacer()) + plot_layout(widths = c(1,.001))

(top / second / third / bottom) + plot_layout(heights = c(1,.5,1,.5))
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
#ggsave("output/Figures/Invitro/invitroFig.pdf", width = 30, height = 40)
```

# Supplementary Plots

## Leiden Cluster plot

``` r
dimLeidenFig <- DimPlotCustom(invitroInvivo, group.by = "leiden_clusters", ncol = 1, pt.size = fig.pt.size, label = T) & theme_bw() & theme_manuscriptOG & NoLegend()  & ggtitle("In Vitro and In Vivo Leiden Clusters")

dimLeidenFig
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

## Split by line dim plot

``` r
dimLineFig <- DimPlotCustom(invitro, split.by = "line", group.by = "cellType", ncol = 1, pt.size = fig.pt.size, label = T) & theme_bw() & theme_manuscriptOG & NoLegend() & theme(axis.title.y = element_blank()) & scale_fill_manual(values = manuscriptPalette)

dimLineFig
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
table(invitro$line)
```

    ## 
    ##   C27  WA09 
    ## 16808 20997

``` r
dimLineFig[[1]] <- dimLineFig[[1]] + labs(title = "iPSC (C27) - 16,808 In Vitro Cells") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
dimLineFig[[2]] <- dimLineFig[[2]] + xlab("UMAP 1") + labs(title = "ESC (WA09) - 20,997 In Vitro Cells")
dimLineFig
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

## Sort

``` r
StackedSort <- as.data.frame(table(invitro$sort, invitro$cellType))
StackedSort$Var1 <- relevel(StackedSort$Var1, "DAPI-")

unique(StackedSort$Var1)
```

    ## [1] CD140a+ DAPI-  
    ## Levels: DAPI- CD140a+

``` r
StackedSortPlot <- ggplot(StackedSort, aes(fill = Var2, y = Freq, x = Var1))+
  geom_bar(position = "fill", stat = "identity")+
  guides(fill = guide_legend(override.aes = list(size = .6))) + theme_bw() + theme_manuscriptOG + scale_y_continuous(expand = c(0,0)) + ylab("Percent Identity") + labs(y = "Percent", x = "Stage", fill = "Cell Type", tag = "C") + theme(legend.key.size = unit(.6, "lines")) + xlab(element_blank()) + NoLegend() & scale_fill_manual(values = manuscriptPalette)  & ggtitle("Sort Cell Type Fractions")

StackedSortPlot
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

## Cell Cycle Scoring

``` r
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#MLF1IP
s.genes[s.genes %not in% invitro@assays$RNA@counts@Dimnames[[1]]] <- "CENPU"
s.genes[s.genes %not in% invitro@assays$RNA@counts@Dimnames[[1]]]
```

    ## character(0)

``` r
#FAM64A HN1
g2m.genes[g2m.genes %not in% invitro@assays$RNA@counts@Dimnames[[1]]] <- c("PIMREG", "JPT1")
g2m.genes[g2m.genes %not in% invitro@assays$RNA@counts@Dimnames[[1]]]
```

    ## character(0)

``` r
invitro <- CellCycleScoring(invitro, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


DimPlotCustom(invitro, split.by = "Phase", label = T, group.by = "cellType") 
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
cellCycleVlnGG <- VlnPlot(invitro, c("G2M.Score", "S.Score"), pt.size = 0, idents = c("NPC", "GPC1", "GPC2", "GPC3", "GPC4")) & theme_bw() & theme_manuscriptOG & NoLegend() & theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

cellCycleVlnGG[[1]] <- cellCycleVlnGG[[1]] + ggtitle("G2M Phase Score")
cellCycleVlnGG[[2]] <- cellCycleVlnGG[[2]] + ggtitle("S Phase Score")

cellCycleVlnGG
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
color0 = "grey60"
colPalette = c("dodgerblue2", 
        "gold", "red2")

cellCycleFeatureGG <- FeaturePlot(invitro, c("G2M.Score", "S.Score"), order = T, pt.size = fig.pt.size) & scale_colour_gradientn(colors = c(color0, colorRampPalette(colPalette)(100))) & theme_bw() & theme_manuscriptOG & theme(legend.position = "bottom", legend.key.width = unit(.5, "inch")) & xlab("UMAP 1") & ylab("UMAP 2")
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

``` r
cellCycleFeatureGG[[1]] <- cellCycleFeatureGG[[1]] + ggtitle ("G2M Phase Score")
cellCycleFeatureGG[[2]] <- cellCycleFeatureGG[[2]] + ggtitle ("S Phase Score")


cellCycleFeatureGG
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->

## Mature NPC Feature Plots

``` r
npcFeaturePlots <- FeaturePlotCustom(invitro, c("SYT1", "MEIS2", "BCL11B"), plotLegend = "shared", pt.size = fig.pt.size)   & theme_bw() & theme_manuscriptOG & theme(legend.position = "bottom", legend.key.width = unit(.5, "inch"))

npcFeaturePlots[[1]] <- npcFeaturePlots[[1]] + ylab("UMAP 2")
npcFeaturePlots <- npcFeaturePlots & xlab("UMAP 1")

npcFeaturePlots
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

## Supplemental Figure

``` r
dimLeidenFig <- dimLeidenFig + labs(tag = "A")
dimLineFig[[1]] <- dimLineFig[[1]] + labs(tag = "B")
dimLineFig[[2]] <- dimLineFig[[2]] + labs(tag = "C")
StackedSortPlot <- StackedSortPlot + labs(tag = "D")
cellCycleFeatureGG[[1]] <- cellCycleFeatureGG[[1]] + labs(tag = "E")
cellCycleVlnGG[[1]] <- cellCycleVlnGG[[1]] + labs(tag = "F")
npcFeaturePlots[[1]] <- npcFeaturePlots[[1]] + labs(tag = "G")


top <- dimLeidenFig | dimLineFig[[1]] | dimLineFig[[2]] | StackedSortPlot

middle <- cellCycleFeatureGG | cellCycleVlnGG

bottom <- npcFeaturePlots

top / middle / bottom
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
#ggsave("output/Figures/Invitro/invitroSuppFig.pdf", width = 30, height = 40)
```

## Specificity Analysis

``` r
invivo <- subset(invitroInvivo, subset = stage == "In Vivo")
invivo <- subset(invivo, subset = cellType %in% c("GPC4", "imOL", "maOL", "Astrocyte"))

GPC4expression <- DotPlot(invitroDE, features = GPC4.vs.Rest.sig[GPC4.vs.Rest.sig$logFC > 0,]$gene)$data
GPC4expressionWide <- pivot_wider(GPC4expression,  names_from = id, id_cols = features.plot, values_from = pct.exp)
GPC4expressionWide$spec.GPC4 <- GPC4expressionWide$GPC4 / rowSums(GPC4expressionWide[,2:6])
GPC4expressionWide <- GPC4expressionWide[order(GPC4expressionWide$spec.GPC4, decreasing = T),]
GPC4expressionWide$group <- "GPC4"

GPC3expression <- DotPlot(invitroDE, features = GPC3.vs.Rest.sig[GPC3.vs.Rest.sig$logFC > 0,]$gene)$data
GPC3expressionWide <- pivot_wider(GPC3expression,  names_from = id, id_cols = features.plot, values_from = pct.exp)
GPC3expressionWide$spec.GPC3 <- GPC3expressionWide$GPC3 / rowSums(GPC3expressionWide[,2:6])
GPC3expressionWide <- GPC3expressionWide[order(GPC3expressionWide$spec.GPC3, decreasing = T),]
GPC3expressionWide$group <- "GPC3"

GPC2expression <- DotPlot(invitroDE, features = GPC2.vs.Rest.sig[GPC2.vs.Rest.sig$logFC > 0,]$gene)$data
GPC2expressionWide <- pivot_wider(GPC2expression,  names_from = id, id_cols = features.plot, values_from = pct.exp)
GPC2expressionWide$spec.GPC2 <- GPC2expressionWide$GPC2 / rowSums(GPC2expressionWide[,2:6])
GPC2expressionWide <- GPC2expressionWide[order(GPC2expressionWide$spec.GPC2, decreasing = T),]
GPC2expressionWide$group <- "GPC2"

GPC1expression <- DotPlot(invitroDE, features = GPC1.vs.Rest.sig[GPC1.vs.Rest.sig$logFC > 0,]$gene)$data
GPC1expressionWide <- pivot_wider(GPC1expression,  names_from = id, id_cols = features.plot, values_from = pct.exp)
GPC1expressionWide$spec.GPC1 <- GPC1expressionWide$GPC1 / rowSums(GPC1expressionWide[,2:6])
GPC1expressionWide <- GPC1expressionWide[order(GPC1expressionWide$spec.GPC1, decreasing = T),]
GPC1expressionWide$group <- "GPC1"

NPCexpression <- DotPlot(invitroDE, features = NPC.vs.Rest.sig[NPC.vs.Rest.sig$logFC > 0,]$gene)$data
NPCexpressionWide <- pivot_wider(NPCexpression,  names_from = id, id_cols = features.plot, values_from = pct.exp)
NPCexpressionWide$spec.NPC <- NPCexpressionWide$NPC / rowSums(NPCexpressionWide[,2:6])
NPCexpressionWide <- NPCexpressionWide[order(NPCexpressionWide$spec.NPC, decreasing = T),]
NPCexpressionWide$group <- "NPC"

new_names <- c("Gene", "Specificity", "Population") 
allSpecificity <- bind_rows(
  GPC4expressionWide %>% 
    select(1,7:8) %>% 
    setNames(new_names),
  
  GPC3expressionWide %>% 
    select(1,7:8) %>% 
    setNames(new_names),
  
  GPC2expressionWide %>% 
    select(1,7:8) %>% 
    setNames(new_names),
  
  GPC1expressionWide %>% 
    select(1,7:8) %>% 
    setNames(new_names),
  
  NPCexpressionWide %>% 
    select(1,7:8) %>% 
    setNames(new_names),
)

allSpecificity$Population <- factor(allSpecificity$Population, levels = rev(c("NPC", "GPC1", "GPC2", "GPC3", "GPC4")))

write.csv(allSpecificity, "output/DE/invitroSpecificity.csv", quote = F, row.names = F)
```

## Ridge Plot

``` r
fig.ridgePlot <- ggplot(allSpecificity, aes(x = Specificity, y = Population, fill = Population)) + 
  geom_density_ridges(alpha = 1) + 
  scale_fill_manual(values = manuscriptPalette) + 
  labs(
    title = 'Density of Cell Population Specificity Scores',
    x = 'Cell Population Specificity', tag = "H") +
  theme_manuscript() + geom_point(shape = 21, fill = "white", colour = "black", stroke = .5, size = 2) + geom_vline(xintercept = .7) + theme(legend.position = "none", axis.title.y = element_blank()) 

fig.ridgePlot
```

    ## Picking joint bandwidth of 0.028

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
invitroDE$label <- paste0("In Vitro ", invitroDE$cellType)
invivo$label <- paste0("In Vivo ", invivo$cellType)

dotplotSeurat <- merge(invitroDE, invivo)
dotplotSeurat$label <- factor(dotplotSeurat$label, levels = rev(c("In Vitro NPC", "In Vitro GPC1", "In Vitro GPC2",
                                                              "In Vitro GPC3", "In Vitro GPC4", "In Vivo GPC4",
                                                              "In Vivo imOL", "In Vivo maOL", "In Vivo Astrocyte")))

specificExpression <- DotPlot(dotplotSeurat, group.by = "label", features = allSpecificity[allSpecificity$Specificity > .7,]$Gene)$data


fig.MarkerDotPlot <- ggplot(specificExpression[specificExpression$id %not in% c("In Vivo imOL", "In Vivo maOL"),], aes(size = pct.exp, fill = avg.exp.scaled, y = id, x = features.plot)) + 
  geom_point(color = "black", pch = 21) + 
  scale_size_area(max_size = 5) + 
  scale_fill_gradientn(colors = PurpleAndYellow()) + 
  theme_bw() + 
  theme_manuscript() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "bottom", title.theme = element_text(size = baseSize*axisTitleSize), barwidth = 5, barheight = .5), 
         size = guide_legend(title.position = "bottom", title.theme = element_text(size = baseSize*axisTitleSize))) +
  labs(tag = "I", title = "High Subpopulation Specificity Gene Expression", size = "% Expressed", fill = "Scaled Expression") 

fig.MarkerDotPlot
```

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
bottomS2 <- free(fig.ridgePlot) | fig.MarkerDotPlot

bottomS2
```

    ## Picking joint bandwidth of 0.028

![](14_Invitro_Analysis_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

``` r
#ggsave(bottomS2, file = "output/Figures/Invitro/bottomS2.pdf", width = 11, height = 5, units = "in")
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
    ##  [1] ggridges_0.5.7              EnhancedVolcano_1.16.0     
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
    ## [29] tidyr_1.3.0                 ggplot2_3.5.2              
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
    ##  [64] generics_0.1.3         evaluate_0.20          stringr_1.5.0         
    ##  [67] fastmap_1.1.1          yaml_2.3.7             goftest_1.2-3         
    ##  [70] knitr_1.42             fitdistrplus_1.1-8     purrr_1.0.1           
    ##  [73] RANN_2.6.1             pbapply_1.7-0          future_1.32.0         
    ##  [76] nlme_3.1-162           mime_0.12              ggrastr_1.0.2         
    ##  [79] compiler_4.2.3         rstudioapi_0.14        beeswarm_0.4.0        
    ##  [82] plotly_4.10.1          png_0.1-8              spatstat.utils_3.1-0  
    ##  [85] tibble_3.2.1           stringi_1.7.12         highr_0.10            
    ##  [88] lattice_0.21-8         Matrix_1.6-4           vctrs_0.6.1           
    ##  [91] pillar_1.9.0           lifecycle_1.0.3        spatstat.geom_3.2-9   
    ##  [94] lmtest_0.9-40          RcppAnnoy_0.0.20       cowplot_1.1.1         
    ##  [97] bitops_1.0-7           irlba_2.3.5.1          httpuv_1.6.9          
    ## [100] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20    
    ## [103] gridExtra_2.3          vipor_0.4.7            parallelly_1.35.0     
    ## [106] codetools_0.2-19       MASS_7.3-58.3          xlsxjars_0.6.1        
    ## [109] rprojroot_2.0.3        withr_2.5.0            sctransform_0.3.5     
    ## [112] GenomeInfoDbData_1.2.9 mgcv_1.8-42            parallel_4.2.3        
    ## [115] grid_4.2.3             rmarkdown_2.21         Rtsne_0.16            
    ## [118] spatstat.explore_3.2-7 shiny_1.7.4            ggbeeswarm_0.7.2
