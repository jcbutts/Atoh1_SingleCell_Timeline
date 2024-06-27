ModScoreRidgePlot <- function(mod.scores.df) {
  # Plots phenotype specific ridgeline plots for celltype module scores
  #
  # Parameters:
  #   mod.scores.df: dataframe with module scores
  #   ###TODO: check with Ryan if downsample is necessary.
  #p.vals <- CalcModPvals(df, downsample.n = 300) #should this be downsample.n = wilcox.downsample?
  #p.vals$`Module score` = 4.5
  
  p <- mod.scores.df %>% 
    ggplot(
      aes(x= Cluster1, 
          y= celltype)) + 
    geom_density_ridges(
      rel_min_height = 0.01, 
      aes(fill = celltype, 
          col = celltype),
      size = 0.25, 
      alpha = 0.75) +
    scale_fill_manual(
      values = rev(Category.cols)) +
    scale_color_manual(
      values = rev(Category.cols)) +
    ggpubr::theme_pubr(base_size = 7, x.text.angle = 45) + 
    ylab("") + 
    xlab("Module score") + 
    NoLegend()
   # geom_text(
    #  data = p.vals,
     # aes(label = p.anno),
    #  size = 7 / .pt,
     # nudge_y = 0.1)
  
  p
}
