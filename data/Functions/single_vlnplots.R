single_vlnplots <- function(obj, gene_list){
  
  SinglePlotTheme <-  theme(axis.title =element_text(size=30), axis.text=element_text(size=30), plot.title=element_text(size=50), legend.position = "none" ) 
  counter = 1
  for (i in gene_list) {
    name <- paste("Plots/Single_Plot/","vln_", as.character(gene_list[counter]),".png")
    
    p <-VlnPlot(object= obj, features = gene_list[counter], pt.size = 0, ncol = 1) + SinglePlotTheme
    save_plot(name, p, base_width = 12, base_height = 8)
    counter = counter + 1
  }
  
  
}
