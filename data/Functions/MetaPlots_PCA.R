MetaPlots_PCA <- function(obj) {

png("Plots/PCA_cluster_withlegend.png", width = 1200, height = 1000)
print(DimPlot(obj, pt.size = 0.5, reduction = "pca", label = TRUE, label.size = 16))
dev.off()

png("Plots/PCA_byTimpoint.png", width = 1200, height = 1000)
print(DimPlot(object = obj, reduction = 'pca', group.by="TimePoint", cols = TimepointCols, pt.size = 0.5, label.size = 3) + theme(axis.text=element_text(size=16),legend.text=element_text(size=20)) + guides(colour = guide_legend(override.aes = list(size=5)))) 
dev.off()

png("Plots/PCA_byColor.png", width = 1200, height = 1000)
print(DimPlot(object = obj, reduction = 'pca', group.by="color", cols = ColorCols, pt.size = 0.5, label.size = 3) + theme(axis.text=element_text(size=16),legend.text=element_text(size=20)) + guides(colour = guide_legend(override.aes = list(size=5))))
dev.off()

png("Plots/PCA_byorigident.png", width = 1800, height = 1000)
print(DimPlot(object = obj, reduction = 'pca', group.by="orig.ident", cols = OrigIdentCols2, pt.size = 0.5, label.size = 3) + theme(axis.text=element_text(size=16),legend.text=element_text(size=15)) + guides(colour = guide_legend(override.aes = list(size=5))))
dev.off()

png("Plots/PCA_bySample.png", width = 1200, height = 1000)
print(DimPlot(object = obj, reduction = 'pca', group.by="Sample", cols = SampleCols, pt.size = 0.5, label.size = 3) + theme(axis.text=element_text(size=16),legend.text=element_text(size=20)) + guides(colour = guide_legend(override.aes = list(size=5))))
dev.off()

png("Plots/PCA_byReplicate.png", width =1200, height = 1000)
print(DimPlot(object = obj, reduction = 'pca', group.by="Replicate", pt.size = 0.5, label.size = 3) + theme(axis.text=element_text(size=16),legend.text=element_text(size=15)) + guides(colour = guide_legend(override.aes = list(size=5))))
dev.off()

png("Plots/PCA_splitbySample.png", width = 3000, height = 400)
print(DimPlot(object = obj, reduction = 'pca', group.by = "orig.ident", split.by="orig.ident", cols = OrigIdentCols2, pt.size = 0.05, label.size = 3)) 
dev.off()

}
