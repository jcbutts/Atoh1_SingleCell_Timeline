
################################################################################
# Purpose: Analyze and plot module score                                       #
# Initial Author: Ryan Dhindsa 
# Modified by J. Butts and R. Dhinsda
################################################################################


CalcModPvals <- function(module.df) {
  # Calculates p-values per module
  res <- data.frame( Cluster = character(), 
                     p = numeric())
  
  x <- module.df
  
  for(cluster in unique(x$celltype)){
    ds <- x %>%
      mutate(curr_cluster = ifelse(celltype == cluster, T, F))
      # group_by(curr_cluster) %>%
      #slice_sample(n=500)
    
    wilcox.res <- wilcox.test(`Cluster1` ~ curr_cluster, 
                              data = ds, 
                              alternative = "less")
    new.row <- data.frame(celltype = cluster, 
                          p = wilcox.res$p.value)
    res <- rbind(res, new.row)
  }
  colnames(res) <- c("celltype", "p")
  res$FDR <- p.adjust(res$p, method = "fdr")
  res$p.adj <- p.adjust(res$p, method = "bonferroni")
  
  res$p.adj <- signif(res$p.adj, 2)  # round to 2 digits
  res$p.anno <- ifelse(res$FDR < 0.01, "*", " ")
  
  return(res)
}