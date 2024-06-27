#library(ggplot2)
#library(RColorBrewer)
#library(viridis)
#library(ggpubr)
#library(readxl)
#library(reshape)
#library(cowplot)
#library(plotly)
#library(ggridges)


MakeTimePlot <- function(Timing_df) {

ident_label <- c("E9.5_GFP", "E10.5_GFP", 
                             "E11.5_GFPTd", "E12.5_GFPTd",
                             "E13.5_GFPTd","E14.5_GFPTd", "E16.5_GFPTd")
variable_label <- c("RLS", "CB_nuclei","EGL1", "EGL2", "Early_Pons","CLS", "Caudal_Pons", "PES", "CES", "AES_early", "AES_late")

#Colors to match migration streams
cols =c(Early_RLS_col, CBNuc_cols, EGL1_cols, EGL2_cols, Early_pons_col, CLS_col, Caudal_pons_col, PES_col, CES_col,AES_col, AES_late_col)

label.setup <- theme(axis.text.x = element_text(size = 14, face = 'bold', angle = 45, hjust = 1, vjust = 1),
                     axis.text.y = element_text(size = 14),
                     axis.title = element_text(size = 25, face = "bold"),
                     plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
                     legend.title=element_text(size=20, face = "bold"), legend.text=element_text(size=14))

#Order labels
MigrationStream$variable <- factor(MigrationStream$variable, levels = variable_label)
MigrationStream$orig.ident <- factor(MigrationStream$orig.ident, levels = ident_label)

#Linegraph
p <-  ggplot(MigrationStream, aes(x = orig.ident, y = value, group = variable, colour = variable, fill = variable)) +                                    
  geom_area(alpha = 0.7) + geom_line(size = 1.5) + scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  facet_wrap(~variable, nrow = 11, strip.position = 'right') +
  theme(legend.position = 'none') +
  scale_x_discrete(limits = ident_label) +
  scale_y_continuous(breaks = c(0, 4000)) +
  theme_classic() +
  theme(strip.text = element_blank()) +
  label.setup 
p <- p + labs(x = "Sample", y = "Number of cells in cluster")
return(p)

}