MinDist <- function(obj, dims, res, spread) {
  # Plots varying valudes of min.dist for a set spread value
  #
  # Args:
  #   obj: Seurat Object
  #   dims: number of PC dimensions 
  #   res: cluster resolution
  #   spread: spread value to use for RUNUMAP
  #
  # Returns:
  #   saved UMAP plots
  
  
  counter = 1
  values <- c(0.01, seq(from = 0.05, to = 0.5, by = 0.05))
  for (i in values) {
    
    
    obj <- RunUMAP(obj, spread = spread, min.dist = values[counter],  dims = 1:dims, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
    obj <- FindClusters(obj, verbose = FALSE, resolution= res)
    
    name <- paste("MinDist/spread",as.character(spread),"min.dist", as.character(values[counter]),".png")
    png(name, width = 1200, height = 1000)
    print(DimPlot(obj, pt.size = 0.5, label = TRUE, label.size = 16) + theme(axis.text=element_text(size=16)))
    dev.off()
    
    counter = counter+1
  }
}