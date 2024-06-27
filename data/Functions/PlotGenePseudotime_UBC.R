library(tidyverse)

# Functions ---------------------------------------------------------------------
GetPseudotime <- function(obj, traj.name, genelist) {
  df <- FetchData(obj, c(traj.name, genelist))
  names(df)[names(df) == traj.name] <- 'pseudotime'
  df <- df %>% 
    drop_na(pseudotime) %>% 
    mutate(pseudotime_perc = percent_rank(pseudotime)) %>% 
    select(-c(pseudotime))
    
  df$trajectory <- traj.name
  df <- df %>% 
    reshape2::melt(id = c("trajectory", "pseudotime_perc"))
  
  return(df)
  
}


GetMultipleTraj <- function(obj, traj.list, genelist) {
  df.list <- as.list(c())
    for (i in 1:length(traj.list)) {
      temp.df <- GetPseudotime(obj, traj.list[i], genelist)
      df.list[[i]] <- temp.df 
    }
  
  full.df <- do.call(rbind, df.list)
  return(full.df)
  
}

PlotGenePseudotime_UBC <- function(traj.dat) {
  p <- ggplot(traj.dat, aes(x = pseudotime_perc, y = value, col = trajectory)) + geom_smooth(method = 'loess') +
    coord_cartesian(ylim = c(0, NA), expand = TRUE) +
    xlab("Pseudotime") + ylab("Expression Level") +
    facet_wrap(~variable, ncol = 1, strip.position = 'right', scale = 'free_y') +
    scale_y_continuous(n.breaks = 3) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(size = 10, face = 'bold'),
          legend.position = 'none',
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20, face = 'bold'))
    
  p
  
}

# Main --------------------------------------------------------------------------


#traj.cols = cols =c("orchid","darkorchid1","darkorchid3","darkorchid4", "turquoise","lightskyblue","lightskyblue2","dodgerblue","mediumblue","navy","midnightblue")

  
#p <- ggplot(traj.dat, aes(x = pseudotime_perc, y = value, col = trajectory)) + geom_smooth() +
 # scale_color_manual(values = traj.cols) +
 # facet_wrap(~variable, ncol = 1, strip.position = 'right') +
 # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 # plot(p)


#df <- df %>% 
 # mutate(cluster2_pt_per = percent_rank(cluster2_pseudotime)) %>% 
# mutate(cluster5_pt_per = percent_rank(cluster5_pseudotime))


#df2 <- df %>% 
 # select(-c('cluster2_pseudotime','cluster5_pseudotime')) %>% 
#  reshape2::melt(value.name = c('cluster2_pt_per','cluster5_pt_per'))

#p <- ggplot(df, aes(x = cluster2_pseudotime, y = Dll3)) + 
 # geom_point() +
  #geom_smooth(method = 'loess')
#plot(p)