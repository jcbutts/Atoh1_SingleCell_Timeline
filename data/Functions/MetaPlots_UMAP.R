MetaPlots_UMAP <- function(obj, ptsize) {
  #A function to create plots to plot the metadata onto UMAp 
  # Arg: 
    # obj: a Seurat object
  #
  # Returns:
   # png files containing metadata umaps
UMAPPlotTheme <- theme(axis.title =element_text(size=30), axis.text=element_text(size=30),legend.text=element_text(size=25, face = 'bold'))  + NoAxes()

p <- DimPlot(object = obj, reduction = 'umap', group.by="TimePoint", cols = TimepointCols, pt.size = ptsize, label.size = 3) + UMAPPlotTheme + guides(colour = guide_legend(override.aes = list(size=12)))
save_plot("MetaPlots_UMAP/UMAP__byTimpoint.png", p, base_width = 9, base_height = 8)

p <- DimPlot(object = obj, reduction = 'umap', group.by="color", cols = ColorCols, pt.size = ptsize, label.size = 3) + UMAPPlotTheme + guides(colour = guide_legend(override.aes = list(size=12)))
save_plot("MetaPlots_UMAP/UMAP_byColor.png", p, base_width = 9, base_height = 8)

p <- DimPlot(object = obj, reduction = 'umap', group.by="orig.ident", cols = OrigIdentCols2, pt.size = ptsize, label.size = 3) + UMAPPlotTheme + guides(colour = guide_legend(override.aes = list(size=12)))
save_plot("MetaPlots_UMAP/UMAP_byorigident.png", p, base_width = 16, base_height = 8)

p <- DimPlot(object = obj, reduction = 'umap', group.by="Sample", cols = SampleCols, pt.size = ptsize, label.size = 3) + UMAPPlotTheme+ guides(colour = guide_legend(override.aes = list(size=12))) 
save_plot("MetaPlots_UMAP/UMAP_bySample.png", p, base_width = 10, base_height = 8)

p <- DimPlot(object = obj, reduction = 'umap', group.by="Replicate", pt.size = ptsize, label.size = 3) + UMAPPlotTheme + guides(colour = guide_legend(override.aes = list(size=12)))
save_plot("MetaPlots_UMAP/UMAP_byReplicate.png", p, base_width = 8, base_height = 8)

p <- DimPlot(object = obj, reduction = 'umap', group.by = "orig.ident", split.by="orig.ident", cols = OrigIdentCols2, pt.size = ptsize) + theme(text = element_text(size = 25, angle =60, face = 'bold'),legend.text=element_text(size=25, face = 'bold', angle = 0)) + NoAxes() + guides(colour = guide_legend(override.aes = list(size=6)))  
save_plot("MetaPlots_UMAP/UMAP_splitbySample.png", p, base_width = 45, base_height = 8)

p <- VlnPlot(object = obj, feature = "nFeature_RNA", pt.size = 0) + theme(axis.title =element_text(size=30), axis.text=element_text(size=30),legend.text=element_text(size=25), plot.title=element_text(size=50))
save_plot("MetaPlots_UMAP/nFeature_ByCluster.png", p, base_width = 16, base_height = 8)

}