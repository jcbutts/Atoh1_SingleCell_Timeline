# Atoh1 Identity

Atoh1_Identity <- function(obj, ptsize) {

  SinglePlotTheme <-  theme(axis.title =element_text(size=30), axis.text=element_text(size=30), plot.title=element_text(size=50),legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=30) ) 
  
  
  
rRLgenes <- c("En1","En2")
cRLgenes <- c("Hoxa2","Hoxb5","Hoxa7")
CellStategenes = c("Mki67", "Nes", "Pax6", "Neurod6", "Rbfox3","Mapt")
ThreeCellState = c("Mki67", "Nes","Mapt")
AtohTargetgenes = c("Barhl1", "Barhl2", "Lhx9","Lhx2")
AtohNTgenes = c("Slc17a7", "Slc17a6", "Sst","Penk","Crh","Nos1","Slc18a3","Th","Tacr1","Slc17a8")
# Penk = Penk1, Slc18A3 = Vacht, Tacr1 = Nk1R
ExcitatoryGenes <- c("Slc17a7","Slc17a6")
InhibitoryGenes <- c("Gad1","Slc6a5")
GlialMarker=c("Apoe","Fabp7","Dbi","Plp1")
Notch = c("Notch1","Notch2","Notch3","Notch4","Dll1","Dll3","Dll4","Jag1","Jag2","Hes1","Hes5","Rbpj")

p <- FeaturePlot(object= obj, features = rRLgenes, pt.size = ptsize, ncol = 2, cols = firemap, combine = FALSE )
p1 <- p[[1]] + SinglePlotTheme + NoAxes()
p2 <- p[[2]] + SinglePlotTheme + NoAxes()
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/rRLIdentity.png", p, base_width = 16, base_height = 8)

p <- FeaturePlot(object= obj, features = cRLgenes, pt.size = ptsize, ncol = 3, cols = firemap, combine = FALSE )
p1 <- p[[1]] + SinglePlotTheme + NoAxes()
p2 <- p[[2]] + SinglePlotTheme + NoAxes()
p3 <- p[[3]] + SinglePlotTheme + NoAxes()
p <-  ggarrange(p1, p2, p3, ncol = 3, common.legend = TRUE, legend = 'right')
save_plot("Plots/cRLIdentity.png", p, base_width = 24, base_height = 8)

p <- DotPlot(obj,features =  c("Hoxc8","Hoxb7","Hoxc6","Hoxa5","Hoxa4","Hoxb3","Hoxd3","Hoxa3","Hoxb2","Hoxa2","En1","En2"))
save_plot("Plots/Hox_DotPlot.png", p, base_width = 20, base_height = 16)

p <- FeaturePlot(object= obj, features = CellStategenes, pt.size = ptsize, ncol = 3, cols = firemap )
save_plot("Plots/CellState.png", p, base_width = 12, base_height = 8)

p <- FeaturePlot(object= obj, features = AtohTargetgenes, pt.size = ptsize, cols = firemap, combine = FALSE )
p1 <- p[[1]] + SinglePlotTheme + NoAxes()
p2 <- p[[2]] + SinglePlotTheme + NoAxes()
p3 <- p[[3]] + SinglePlotTheme + NoAxes()
p4 <- p[[4]] + SinglePlotTheme + NoAxes()
p <-  ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = 'right')
save_plot("Plots/AtohhTargetGenes.png", p, base_width = 32, base_height = 16)

p <- FeaturePlot(object= obj, features = AtohNTgenes, pt.size = ptsize, ncol = 3, cols = firemap )
save_plot("Plots/AtohNTTarget.png", p, base_width = 16, base_height = 20)

p <- FeaturePlot(object= obj, features = c(ExcitatoryGenes, InhibitoryGenes), pt.size = ptsize, ncol = 2, cols = firemap )
save_plot("Plots/ExIn.png", p, base_width = 8, base_height = 8)

p <- FeaturePlot(object= obj, features = GlialMarker, pt.size = ptsize, ncol = 2, cols = firemap )
save_plot("Plots/Glia.png", p, base_width = 8, base_height = 8)

p <- FeaturePlot(object= obj, features = Notch, pt.size = ptsize, ncol = 3, cols = firemap )
save_plot("Plots/NotchGenes.png", p, base_width = 12, base_height = 16)

p <- FeaturePlot(object= obj, features = c("Hbb-bs", "Hbb-y","Hba-a1","Pecam1"), pt.size = ptsize, ncol = 4, cols = firemap )
save_plot("Plots/BloodGenes.png", p, base_width = 16, base_height = 4)

p <- FeaturePlot(object= obj, features = c("Xist"), pt.size = ptsize, ncol = 1, cols = firemap )
save_plot("Plots/Xist.png", p, base_width = 4, base_height = 4)

p <- FeaturePlot(object= obj, features = c("Atoh1"), pt.size = ptsize, ncol = 1, cols = firemap )
save_plot("Plots/Atoh1.png", p, base_width = 4, base_height = 4)

p <- FeaturePlot(object= obj, features = c("Atoh1"), pt.size = ptsize, ncol = 1, cols = SampleCols )
save_plot("Plots/Atoh1_samplecols.png", p, base_width = 4, base_height = 4)

p <- FeaturePlot(object= obj, features = c("WPRE"), pt.size = ptsize, ncol = 1, cols = firemap )
save_plot("Plots/WPRE.png", p, base_width = 4, base_height = 4)

##Violin Plots
p <- VlnPlot(object= obj, features = rRLgenes, pt.size = 0, ncol = 2)
save_plot("Plots/Vln_rRLIdentity.png", p, base_width = 16, base_height = 4)

p <- VlnPlot(object= obj, features = cRLgenes, pt.size = 0, ncol = 3)
save_plot("Plots/Vln_cRLIdentity.png", p, base_width = 24, base_height = 4)

p <- VlnPlot(object= obj, features = CellStategenes, pt.size = 0, ncol = 3)
save_plot("Plots/Vln_CellState.png", p, base_width = 24, base_height = 8)

p <- VlnPlot(object= obj, features = AtohTargetgenes, pt.size = 0, ncol = 2)
save_plot("Plots/Vln_AtohhTargetGenes.png", p, base_width = 16, base_height = 8)

p <- VlnPlot(object= obj, features = AtohNTgenes, pt.size = 0, ncol = 3)
save_plot("Plots/Vln_AtohNTTarget.png", p, base_width = 32, base_height = 20)

p <- VlnPlot(object= obj, features = c(ExcitatoryGenes, InhibitoryGenes), pt.size = 0, ncol = 2)
save_plot("Plots/Vln_ExIn.png", p, base_width = 16, base_height = 8)

p <- VlnPlot(object= obj, features = Notch, pt.size = 0, ncol = 3)
save_plot("Plots/Vln_NotchGenes.png", p, base_width = 24, base_height = 16)

p <- VlnPlot(object= obj, features = c("Hbb-bs", "Hbb-y","Hba-a1","Pecam1"), pt.size = 0, ncol = 4)
save_plot("Plots/Vln_BloodGenes.png", p, base_width = 32, base_height = 4)

p <- VlnPlot(object= obj, features = c("Xist"), pt.size = 0, ncol = 1)
save_plot("Plots/Vln_Xist.png", p, base_width = 8, base_height = 4)

p <-VlnPlot(object= obj, features = c("Atoh1"), pt.size = 0, ncol = 1)
save_plot("Plots/Vln_Atoh1.png", p, base_width = 8, base_height = 4)

p <- VlnPlot(object= obj, features = c("WPRE"), pt.size = 0, ncol = 1)
save_plot("Plots/Vln_WPRE.png", p, base_width = 8, base_height = 4)




}