Spread<- function(obj, dims, res, mindist) {
counter = 1
values <- c( 0.5, seq(from = 1, to = 5, by =1  ))
for (i in values) {
  
  
  obj <- RunUMAP(obj, spread = values[counter], min.dist = mindist,  dims = 1:dims, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE, resolution= res)
  
  name <- paste("Spread/mindist",as.character(mindist),"spread", as.character(values[counter]),".png")
  png(name, width = 1200, height = 1000)
  print(DimPlot(obj, pt.size = 0.5, label = TRUE, label.size = 16) + theme(axis.text=element_text(size=16)))
  dev.off()
  
  counter = counter+1
}
}
