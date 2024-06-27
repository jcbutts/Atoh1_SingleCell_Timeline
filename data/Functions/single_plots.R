single_plots <- function(obj, gene_list, pt.size){
  
  SinglePlotTheme <-  theme(axis.title =element_text(size=30), axis.text=element_text(size=30), plot.title=element_text(size=50), legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=30)) 
counter = 1
  for (i in gene_list) {
    name <- paste("Plots/Single_Plot/", as.character(gene_list[counter]),".png")
    
    p <-   FeaturePlot(object= obj, features = gene_list[counter] , pt.size = pt.size, ncol = 1, cols = firemap )+ SinglePlotTheme + NoAxes()
    save_plot(name, p, base_width = 8, base_height = 8)
    counter = counter + 1
  }
  
  
}
