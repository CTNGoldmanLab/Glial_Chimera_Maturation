manuscriptPalette <- c("In Vivo" = "red2", 
                       "In Vitro - GPC Stage" = "#2E30FF",
                       "In Vitro - PSC Stage" = "mediumseagreen",
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

# 
# theme_manuscript <-  theme(axis.text = element_text(size = axisTextSize), 
#                            axis.title = element_text(size = axisTitleSize), 
#                            title = element_text(size = titleFont), 
#                            legend.title = element_text(size = axisTitleSize),
#                            legend.text = element_text(size = axisTitleSize),
#                            plot.tag = element_text(size = tagSize),
#                            plot.title = element_text(size = titleFont))



baseSize = 6
axisTextSize <- 1
axisTitleSize <- 1.125
titleSize <- 1.25
legendTitleSize <- 1.125
legendTextSize <- 1
tagSize = 1.4
labelSize <- baseSize * 1


theme_manuscript <- function() {
  theme_bw(base_size = baseSize) + 
    theme(
      axis.text = element_text(size = rel(axisTextSize)), 
      axis.title = element_text(size = rel(axisTitleSize)), 
      title = element_text(size = rel(titleSize)), 
      legend.title = element_text(size = rel(legendTitleSize)),
      legend.text = element_text(size = rel(legendTextSize)),
      plot.tag = element_text(size = rel(tagSize)),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(0,0,0,0,"pt")
    )
}
