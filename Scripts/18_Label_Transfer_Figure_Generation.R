

## Ramos
integratedRamos <- readRDS("output/RDS/integratedRamos.rds")

round((prop.table(table(integratedRamos$predicted.id, integratedRamos$cellType), margin = 2)*100),1)

integratedRamos$cellType <- factor(integratedRamos$cellType, c("NPC", "GPC1", "GPC2", "GPC3", "GPC4", "imOL"))
integratedRamos$predicted.id <- factor(integratedRamos$predicted.id, rev(c("EN", "IN", "nIPC", "mIPC", "gIPC", "gIPC-A", "gIPC-O", "tRG", 'oRG', "RG/AC", "AC", "AC-p", "AC-f", "OPC", "preOL")))

hmDFRamos <- as.data.frame(round((prop.table(table(integratedRamos$predicted.id, integratedRamos$cellType), margin = 2)*100),1))

hm.ramos.gg <- ggplot(hmDFRamos, aes(x = Var2, y = Var1, fill = Freq)) + geom_tile(colour = "black") + labs(x = "In vitro hGPC Cell Type", y = "Fetal scRNA-seq Cell Type - Ramos 2022") + scale_fill_gradient2(low = "white") + scale_x_discrete(expand = c(0,0))  + scale_y_discrete(expand = c(0,0)) + theme_manuscript + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

dimplot.ramos.gg <- ((DimPlotCustom(subset(integratedRamos, subset = stage == "Fetal"), group.by = "celltypes", label = T) + ggtitle("17-41wk GA Fetal - Ramos 2022") + guides(fill = guide_legend(override.aes = list(size = 5)))) |
    (DimPlotCustom(subset(integratedRamos, subset = stage != "Fetal"), group.by = "cellType", label = T) + ggtitle("In vitro - GPC Stage") + scale_fill_manual(values = manuscriptPalette) + guides(fill = guide_legend(override.aes = list(size = 5))))) & theme_manuscript & labs(fill = "Cell Type")


## van Bruggen
integratedVanBruggen <- readRDS("output/RDS/integratedVanBruggen.rds")

round((prop.table(table(integratedVanBruggen$predicted.id, integratedVanBruggen$cellType), margin = 2)*100),1)

integratedVanBruggen$cellType <- factor(integratedVanBruggen$cellType, c("NPC", "GPC1", "GPC2", "GPC3", "GPC4", "imOL"))
integratedVanBruggen$predicted.id <- factor(integratedVanBruggen$predicted.id, levels = rev(levels(integratedVanBruggen$predicted.id)))

hmDFVB <- as.data.frame(round((prop.table(table(integratedVanBruggen$predicted.id, integratedVanBruggen$cellType), margin = 2)*100),1))

hm.vb.gg <- ggplot(hmDFVB, aes(x = Var2, y = Var1, fill = Freq)) + geom_tile(colour = "black") + labs(x = "In vitro hGPC Cell Type", y = "Fetal scRNA-seq Cell Type - van Brueggen 2022") + scale_fill_gradient2(low = "white") + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + theme_manuscript + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

dimplot.vb.gg <- ((DimPlotCustom(subset(integratedVanBruggen, subset = stage == "Fetal"), group.by = "group", label = T) + ggtitle("8-10wk GA Fetal - Van Brueggen 2022") + guides(fill = guide_legend(override.aes = list(size = 5)))) |
    (DimPlotCustom(subset(integratedVanBruggen, subset = stage != "Fetal"), group.by = "cellType", label = T) + ggtitle("In vitro - GPC Stage") + scale_fill_manual(values = manuscriptPalette) +  guides(fill = guide_legend(override.aes = list(size = 5))))) & theme_manuscript & labs(fill = "Cell Type")



## Velmeshev
integratedVelmeshev <- readRDS("output/RDS/integratedVelmeshev.rds")

round((prop.table(table(integratedVelmeshev$predicted.id, integratedVelmeshev$cellType), margin = 2)*100),1)

integratedVelmeshev$cellType <- factor(integratedVelmeshev$cellType, c("NPC", "GPC1", "GPC2", "GPC3", "GPC4", "imOL"))
integratedVelmeshev$predicted.id <- factor(integratedVelmeshev$predicted.id, rev(c("IN", "EN", "Glial Prog", "AC", "OPC")))


hmDFVelm <- as.data.frame(round((prop.table(table(integratedVelmeshev$predicted.id, integratedVelmeshev$cellType), margin = 2)*100),1))

hm.Velm.gg <- ggplot(hmDFVelm, aes(x = Var2, y = Var1, fill = Freq)) + geom_tile(colour = "black") + labs(x = "In vitro hGPC Cell Type", y = "Fetal scRNA-seq Cell Type - Velmeshev 2023") + scale_fill_gradient2(low = "white") + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + theme_manuscript + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

dimplot.velm.gg <- ((DimPlotCustom(subset(integratedVelmeshev, subset = stage == "Fetal"), group.by = "group", label = T) + ggtitle("13-37wk GA Fetal - Velmeshev 2023 ") + guides(fill = guide_legend(override.aes = list(size = 5)))) |
                    (DimPlotCustom(subset(integratedVelmeshev, subset = stage != "Fetal"), group.by = "cellType", label = T) + ggtitle("In vitro - GPC Stage") + scale_fill_manual(values = manuscriptPalette) +  guides(fill = guide_legend(override.aes = list(size = 5))))) & theme_manuscript & labs(fill = "Cell Type")




(dimplot.ramos.gg | hm.ramos.gg) /
  (dimplot.velm.gg | hm.Velm.gg) /
  (dimplot.vb.gg | hm.vb.gg)



###
