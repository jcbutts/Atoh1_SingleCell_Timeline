#Arguments
#DEG_list = columns labeled "DEG" and "Peak"
#CategoryList = 


PlotEnrichment <- function(DEG_list, CategoryList,perc.overlap) {
  #Create empty dataframe
  df.fisher <- as.data.frame(NA)
  df.fisher[,c("Category","pval","conf.interval1","conf.interval2","odds_ratio","Percent_peaks")] <- NA
  df.percent.peaks <- c()  
  
  #Calculate Fishers test and percent peaks for each cell state  
  counter = 1
  for (i in CategoryList) {
    
    df<- DEG_list %>% 
      filter(cluster == i) 
    df.temp.mat <- as.matrix(table(df$DEG, df$Peak))
    df.percent.peaks <- (df.temp.mat[2,2]/(df.temp.mat[2,2]+df.temp.mat[2,1]))*100
    df.temp <- fisher.test(df.temp.mat)
    df.temp2 <- c(0,i, df.temp$p.value,df.temp$conf.int[[1]],df.temp$conf.int[[2]],df.temp$estimate, df.percent.peaks)
    df.fisher[counter,] <- df.temp2
    
    counter = counter+1
  }
  
  df.fisher <- df.fisher %>% mutate(pval_adj = as.numeric(pval) * length(CategoryList))
  
  #Plot the enrichment
  df.fisher$Category <- factor(df.fisher$Category, level = CategoryList)
  p <- ggplot(df.fisher, aes(x = Category, y = as.numeric(odds_ratio))) + geom_point(color = Category.cols, size = 5) +
    geom_errorbar(aes(ymin=as.numeric(conf.interval1), ymax=as.numeric(conf.interval2)), color = Category.cols, width = 0.3) +
    geom_hline(yintercept=1,linetype=2) +
    geom_text(aes(label = paste0("p=", signif(as.numeric(pval_adj), digits = 2)), y = 3), color = Category.cols, 
              position = position_dodge(width = 0.5), size = 16 / .pt, hjust = 0.5, show.legend = F) +
    geom_text(aes(label = paste0( signif(as.numeric(Percent_peaks) ,digits = 2), "%"), y = 2.85), color = Category.cols, 
              position = position_dodge(width = 0.5), size = 16 / .pt, hjust = 0.5, show.legend = F) +
    geom_text(aes(label = paste0( signif(as.numeric(perc.overlap$n.expressed))), y = 2.70), color = Category.cols, 
              position = position_dodge(width = 0.5), size = 16 / .pt, hjust = 0.5, show.legend = F) +
    labs(y = 'Odds Ratio') + 
    theme_bw() + theme(text = element_text(size = 20, face = 'bold'), axis.text.x =  element_text(angle = 45, hjust = 1, vjust = 1) ) 
  return(p)
}