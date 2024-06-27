##Atoh1 lineage identity

Atoh1_Lineage <- function(obj) {

AESGenes=c("Lhx2", "Hoxb3", "Hoxb4", "Hoxb5","Nhlh1","Nhlh2")
PNGenes=c("Hoxa2","Hoxa4","Pax6","Zic1", "Zic2","Barhl1")
EGLGenes=c("En2", "Pax6","Barhl1","Zic1","Reln","Neurod1")
#OTher potential genes "Sox2","Pax2", "Ntn1", "Zic2"
CNGenes=c("En2", "Tbr1", "Lmx1a","Eomes")
DCNGenes=c("En2", "Pax6", "Tbr1", "Lmx1a", "Eomes", "Pou3f1", "Irx3", "Slc17a7")
PreCBGenes=c("Barhl1","Foxp2","Tlx3","Arhgef2")
CochGenes=c("Hoxa2", "Lhx9", "Lhx2","Atoh7","Mafb","Maf")
RLSGenes=c("En2", "Ntn1", "Robo2","Slit2","Neurod1")
NonRLGenes=c("Lbx1", "Phox2b","Gad2")
PBGenes= c("Atoh1", "En1", "Runx1", "Tacr1", "Calca", "Adcyap1","Th")
RostralPonsGenes = c("En1","Lhx9")

p <- FeaturePlot(object= obj, features = PNGenes, pt.size = 0, ncol = 3, cols = firemap )
save_plot("Plots/PNGenes.png", p, base_width = 12, base_height = 8)

p <- FeaturePlot(object= obj, features = EGLGenes, pt.size = 0, ncol = 3, cols = firemap )
save_plot("Plots/EGLGenes.png", p, base_width = 12, base_height = 8)

p <- FeaturePlot(object= obj, features = CNGenes, pt.size = 0, ncol = 2, cols = firemap )
save_plot("Plots/CNGenes.png", p, base_width = 8, base_height = 8)

p <- FeaturePlot(object= obj, features = DCNGenes, pt.size = 0, ncol = 3, cols = firemap )
save_plot("Plots/DCNGenes.png", p, base_width = 12, base_height = 12)

p <- FeaturePlot(object= obj, features = PreCBGenes, pt.size = 0, ncol = 2, cols = firemap )
save_plot("Plots/PreCBGenes.png", p, base_width = 8, base_height = 8)

p <- FeaturePlot(object= obj, features = CochGenes, pt.size = 0, ncol = 3, cols = firemap )
save_plot("Plots/CochGenes.png", p, base_width = 12, base_height = 8)

p <- FeaturePlot(object= obj, features = RLSGenes, pt.size = 0, ncol = 3, cols = firemap )
save_plot("Plots/RLSGenes.png", p, base_width = 12, base_height = 8)

p <- FeaturePlot(object= obj, features = NonRLGenes, pt.size = 0, ncol = 3, cols = firemap )
save_plot("Plots/NonRLGenes.png", p, base_width = 12, base_height = 4)

p <- FeaturePlot(object= obj, features = PBGenes, pt.size = 0, ncol = 3, cols = firemap )
save_plot("Plots/PBGenes.png", p, base_width = 12, base_height = 12)

p <- FeaturePlot(object= obj, features = AESGenes, pt.size = 0, ncol = 3, cols = firemap )
save_plot("Plots/AESGenes.png", p, base_width = 12, base_height = 8)

p <- FeaturePlot(object= obj, features = RostralPonsGenes, pt.size = 0, ncol = 2, cols = firemap )
save_plot("Plots/RostralPonsGenes.png", p, base_width = 8, base_height = 4)

#Violin
p <- VlnPlot(object= obj, features = PNGenes, pt.size = 0, ncol = 3)
save_plot("Plots/Vln_PNGenes.png", p, base_width = 24, base_height = 8)

p <- VlnPlot(object= obj, features = EGLGenes, pt.size = 0, ncol = 3)
save_plot("Plots/Vln_EGLGenes.png", p, base_width = 24, base_height = 8)

p <- VlnPlot(object= obj, features = CNGenes, pt.size = 0, ncol = 2 )
save_plot("Plots/Vln_CNGenes.png", p, base_width = 16, base_height = 8)

p <- VlnPlot(object= obj, features = DCNGenes, pt.size = 0, ncol = 3 )
save_plot("Plots/Vln_DCNGenes.png", p, base_width = 24, base_height = 12)

p <- VlnPlot(object= obj, features = PreCBGenes, pt.size = 0, ncol = 2 )
save_plot("Plots/Vln_PreCBGenes.png", p, base_width = 16, base_height = 8)

p <- VlnPlot(object= obj, features = CochGenes, pt.size = 0, ncol = 3)
save_plot("Plots/Vln_CochGenes.png", p, base_width = 24, base_height = 8)

p <- VlnPlot(object= obj, features = RLSGenes, pt.size = 0, ncol = 3 )
save_plot("Plots/Vln_RLSGenes.png", p, base_width = 24, base_height = 8)

p <- VlnPlot(object= obj, features = NonRLGenes, pt.size = 0, ncol = 3 )
save_plot("Plots/V;n_NonRLGenes.png", p, base_width = 24, base_height = 4)

p <- VlnPlot(object= obj, features = PBGenes, pt.size = 0, ncol = 3 )
save_plot("Plots/Vln_PBGenes.png", p, base_width = 24, base_height = 12)

p <- VlnPlot(object= obj, features = AESGenes, pt.size = 0, ncol = 3 )
save_plot("Plots/Vln_AESGenes.png", p, base_width = 24, base_height = 8)

p <- VlnPlot(object= obj, features = RostralPonsGenes, pt.size = 0, ncol = 2 )
save_plot("Plots/Vln_RostralPonsGenes.png", p, base_width = 16, base_height = 4)


}
