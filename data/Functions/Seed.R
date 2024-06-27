Seed<- function(obj, res, dims,num_seed) {
  
  counter = 1
  values <- sample.int(500 ,num_seed)
  for (i in values) {
    
    
    obj <- RunUMAP(obj, dims = 1:dims, seed.use = values, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
    obj <- FindClusters(obj, verbose = FALSE, resolution= res)
    
    name <- paste("seed", as.character(values[counter]),".png")
    png(name, width = 1200, height = 1000)
    print(DimPlot(obj, pt.size = 0.5, label = TRUE, label.size = 16) + theme(axis.text=element_text(size=16)))
    dev.off()
    
    counter = counter+1
  }
}