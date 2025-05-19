library(ggplot2)

axisTitleSize <- 24
axisTextSize <- 18
labelFont = 18
titleFont = 22
tagSize = 26

## Make SCA for DE
makeSCA <- function(seurat, cellFractions){
  DefaultAssay(seurat) <- "RNA"
  data.use <- GetAssayData(object = seurat, slot = "data")
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[,1]
  cData <- as.data.frame(seurat@meta.data)
  cData$wellKey <- row.names(cData)
  cData$ngeneson <- scale(cData$nFeature_RNA)
  temp.sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = data.use),
    check_sanity = T,
    cData = cData,
    fData = fdat)
  temp.sca <- temp.sca[as.character(cellFractions),]
  return(temp.sca)
}

## Likelihood Ratio Test for MAST from ZLM
runLR <- function(z, lrContrast, contrast0, contrast1, fileName, FDR = 0.01, logFC = .25){
  temp <- as.data.frame(lrTest(z, as.matrix(lrContrast)))
  temp$FDR <- p.adjust(p = temp$`hurdle.Pr(>Chisq)`, 'fdr')
  temp <- temp[,9:10]
  log2FC <- getLogFC(z, contrast0 = contrast0, contrast1 = contrast1)
  temp$logFC <- log2FC$logFC
  temp <- temp[order(temp$FDR,  decreasing = F),]
  temp$gene <- row.names(temp)
  temp <- temp[,c(4,1,2,3)]
  assign(fileName, value = temp, envir = .GlobalEnv)
  write.table(temp, paste0("output/DE/", fileName, ".txt"), sep = "\t", row.names = F, quote = F)
  temp <- temp[temp$FDR < FDR & abs(temp$logFC) > logFC,]
  write.table(temp, paste0("output/DE/", fileName, ".sig.txt"), sep = "\t", row.names = F, quote = F)
  assign(paste0(fileName, ".sig"), value = temp, envir = .GlobalEnv)
}

## GGPlot Theme for figures
theme_manuscript <- theme_bw() + theme(axis.text = element_text(size = axisTextSize), 
                                       axis.title = element_text(size = axisTitleSize), 
                                       title = element_text(size = titleFont), 
                                       legend.title = element_text(size = titleFont),
                                       legend.text = element_text(size = axisTitleSize),
                                       plot.tag = element_text(size = tagSize))

`%not in%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))

DimPlotPseudotimeKnots <- function(seurat, group.by = "orig.ident", pt.size = 1, curve = sds, models = gamFit){
  p <- DimPlotCustom(seurat, group.by = group.by, pt.size = pt.size)
  for (i in seq_along(slingCurves(curve))) {
    curve_i <- slingCurves(curve)[[i]]
    curve_i <- curve_i$s[curve_i$ord, seq_len(2)]
    #colnames(curve_i) <- c("dim1", "dim2")
    p <- p + geom_path(data = as.data.frame(curve_i), col = "black", size = 1)
    #If I feel like adding an arrow again
    #p <- p + geom_path(data = as.data.frame(curve_i), col = "black", size = 1, arrow = arrow())
  }
  # Adding the knots
  nCurves <- length(slingCurves(curve))
  knots <- S4Vectors::metadata(models)$tradeSeq$knots
  knots_dim <- matrix(ncol = 2, nrow = nCurves * length(knots))
  for (ii in seq_along(slingCurves(curve))) {
    S <- project_to_curve(x = slingCurves(curve)[[ii]]$s,
                          s = slingCurves(curve)[[ii]]$s[slingCurves(curve)[[ii]]$ord, ],
                          stretch = 0)
    for (jj in seq_along(knots)) {
      kn <- knots[jj]
      times <- S$lambda
      knot <- which.min(abs(times - kn))
      knots_dim[jj + (ii-1)*length(knots), ] <- S$s[knot, seq_len(2)]
    }
  }
  knots_dim <- as.data.frame(knots_dim)
  colnames(knots_dim) <- c("UMAP_1", "UMAP_2")
  p <- p +
    geom_point(data = knots_dim, col = "black", size = 2)
  return(p)
}

calcExpressionFractions <- function(seurat, group){
  tempFractions <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(tempFractions) <- c("pct.expr", "features.plot", "id")
  
  for(i in unique(seurat[[group]][[1]])){
    #tempNorm <- as.data.frame(rowMeans((all_GPC@assays$RNA@data[,colnames(all_GPC@assays$RNA@data) %in% row.names(all_GPC@meta.data[all_GPC@meta.data$Group == i,])])))
    tempSparse <- seurat@assays$RNA@counts[,colnames(seurat@assays$RNA@counts) %in% Cells(seurat)[which(seurat[[group]][[1]] %in% i)]]
    tempDF <- as.data.frame(rowSums(tempSparse != 0) * 100 / dim(tempSparse)[2])
    names(tempDF) <- "pct.expr"
    tempDF$features.plot <- row.names(tempDF)
    tempDF$id <- i
    tempFractions <- rbind(tempFractions, tempDF)
  }
  return(tempFractions)
}

## Update of Signac function that doesn't innately filter biotypes
GetGRangesFromEnsDb2 <- function (ensdb, standard.chromosomes = TRUE, verbose = TRUE) {
  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("Please install biovizBase\n", "https://www.bioconductor.org/packages/biovizBase/")
  }
  whole.genome <- as(object = seqinfo(x = ensdb), Class = "GRanges")
  if (standard.chromosomes) {
    whole.genome <- keepStandardChromosomes(whole.genome, 
                                            pruning.mode = "coarse")
  }
  my_lapply <- ifelse(test = verbose, yes = pbapply::pblapply, no = lapply)
  tx <- my_lapply(X = seq_along(whole.genome), FUN = function(x) {
    suppressMessages(expr = biovizBase::crunch(obj = ensdb, 
                                               which = whole.genome[x], columns = c("tx_id", "gene_name", 
                                                                                    "gene_id", "gene_biotype")))
  })
  tx <- do.call(what = c, args = tx)
  return(tx)
}

