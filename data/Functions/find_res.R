#obj = 
#values = 



find_res <- function(obj){
values = seq(from = 0.1, to = 1, by =0.1)
counter = 1
for (i in values) {
  
  
  obj <- FindClusters(obj, verbose = FALSE, resolution= values[counter])
  
  name <- paste("FindRes/UMAP_","res", as.character(values[counter]),".png")
  png(name, width = 1200, height = 1000)
  print(DimPlot(obj, pt.size = 0.5, label = TRUE, label.size = 16) + theme(axis.text=element_text(size=16)))
  dev.off()
  
  counter = counter+1
}
}