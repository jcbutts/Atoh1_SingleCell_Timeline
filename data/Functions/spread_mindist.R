
Spread_mindist <- function(obj, dims, res) {
counter = 1
values <- c(0.01, seq(from = 0.05, to = 0.5, by = 0.05))
for (i in values) {
  
  
  obj <- RunUMAP(obj, spread = 1, min.dist = values[counter],  dims = 1:dims, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE, resolution= res)
  
  name <- paste("spread_1","min.dist", as.character(values[counter]),".png")
  png(name, width = 1200, height = 1000)
  print(DimPlot(obj, pt.size = 0.5, label = TRUE, label.size = 16) + theme(axis.text=element_text(size=16)))
  dev.off()
  
  counter = counter+1
}

#Change spread
counter = 1
values <- c(0.1, 0.5, seq(from = 1, to = 5, by =1  ))
for (i in values) {
  
  
  obj <- RunUMAP(obj, spread = values[counter], min.dist = 0.2,  dims = 1:dims, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE, resolution= res)
  
  name <- paste("mindist_0.2","spread", as.character(values[counter]),".png")
  png(name, width = 1200, height = 1000)
  print(DimPlot(obj, pt.size = 0.5, label = TRUE, label.size = 16) + theme(axis.text=element_text(size=16)))
  dev.off()
  
  counter = counter+1
}
}
