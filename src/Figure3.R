### Figure 3

#Set directory as you see fit
setwd("C:/Jessie/AllSamples/Publication/Figure3")

##Create Folders
dir.create("InitialQC")
dir.create("MinDist")
dir.create("Spread")
dir.create("MetaPlots_UMAP")
dir.create("Plots")
dir.create("FindRes")

###Functions
source('C:/Jessie/Functions/find_res.R')
source('C:/Jessie/Functions/spread.R')
source('C:/Jessie/Functions/mindist.R')
source('C:/Jessie/Functions/MetaPlots_UMAP.R')
source('C:/Jessie/Functions/RenameClusterGroups.R')

###Libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggpubr)
library(data.table)
library(ggplot2)

###Dependencies

PlotTheme <- theme(plot.title=element_text(size=50), legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=30)) 

#Combined SCT Seurat Object
load("C:/Jessie/AllSamples/Publication/Figure1/FullDataset_normalizedRNA.Rdata")

#From the combined Seurat objects folder
#load("C:/Jessie/AllSamples/Publication/Seurat_Objects/FullDataset_normalizedRNA.Rdata")

#Category ID Labels
Category_IDs <- "C:/Jessie/AllSamples/Publication/RenameCluster_Files/Category_IDs.txt"

#####Script#####

###Note, I need the Category IDs for the last part. Would probably be a good idea stash these IDs and load in the SO above.
Idents(CombinedSCT) <- 'seurat_clusters'

#Rename clusters by category
CombinedSCT <- RenameClusterGroups(CombinedSCT, Category_IDs)

##Reorder levels
Category.levels = c("Progenitor", "Intermediate Progenitor", "Migrating_1", "Migrating_2", "Mature","Non RL")

#Add labels to metadata and reorder
CombinedSCT[["FullData.category"]] <- Idents(object = CombinedSCT)
CombinedSCT$FullData.category <- factor(CombinedSCT$FullData.category, levels = Category.levels)

Idents(CombinedSCT) <- 'seurat_clusters'

##Look at cells that express mki67 and where their expression is over1
#Figure3B 
p <- DimPlot(object = CombinedSCT, reduction = 'umap', pt.size = 0.5, cells.highlight = WhichCells(object = CombinedSCT, expression = Mki67 > 1), cols.highlight = 'steelblue4') + 
  ggtitle("Mki67 > 1") +
  theme(plot.title = element_text(size = 40, hjust = 0.5, face = "bold")) +
  NoAxes() + NoLegend()
save_plot("Plots/Figure3_FullData_Mki67_expression_over1.png", p, base_width = 8, base_height = 9)

#stash old identity# 
CombinedSCT[["old.ident"]] <- Idents(object = CombinedSCT)

Prolif_subset <- subset(x = CombinedSCT, subset = Mki67 > 1)

png("InitialQC/VlnQCbySample.png", width = 2000, height = 400)
VlnPlot(object = Prolif_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.blood"), group.by = "orig.ident", cols = OrigIdentCols, ncol = 4, pt.size = 0)
dev.off()

Prolif_subset <- SCTransform(Prolif_subset, vars.to.regress = c("percent.mito","percent.blood"), verbose = FALSE)

Prolif_subset <- RunPCA(Prolif_subset, verbose = FALSE)
png("InitialQC/QC_PCHeatmap.png", width = 500, height = 1000)
DimHeatmap(Prolif_subset, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png("InitialQC/QC_ElbowPlot.png", width = 500, height = 500)
ElbowPlot(Prolif_subset, ndims = 50)
dev.off()

Prolif_subset <- RunUMAP(Prolif_subset, dims = 1:50, verbose = FALSE)
Prolif_subset <- FindNeighbors(Prolif_subset, dims = 1:50, verbose = FALSE)
Prolif_subset <- FindClusters(Prolif_subset, verbose = FALSE, resolution= 0.8)

#Res at 50PC
find_res(Prolif_subset)
MinDist(Prolif_subset, 50, 0.4, 1)
Spread(Prolif_subset, 50, 0.4, 0.3)

#Final Setup
Prolif_subset <- RunUMAP(Prolif_subset, spread = 1, min.dist = 0.3, dims = 1:50, verbose = FALSE)
Prolif_subset <- FindNeighbors(Prolif_subset, dims = 1:50, verbose = FALSE)
Prolif_subset <- FindClusters(Prolif_subset, verbose = FALSE, resolution= 0.3)

MetaPlots_UMAP(Prolif_subset, 1)

#Figure3C
p <- DimPlot(Prolif_subset, group.by = 'Sample', cols = SampleCols, pt.size = 1, label = FALSE) + NoAxes() +
  guides(colour = guide_legend(override.aes = list(size=14))) +
  theme(legend.text = element_text(size = 40, face = 'bold'))
save_plot("Plots/Figure3_UMAP_ProlifSubset_bySample.png", p, base_width = 12, base_height = 8)

#Figure3D
p <- DimPlot(Prolif_subset, pt.size = 1, label = TRUE, label.size = 16) + theme(legend.position = "none") + NoAxes()
save_plot("Plots/Figure3_UMAP_ProlifSubset_byCluster.png", p, base_width = 8, base_height = 8)

save(Prolif_subset, file = "Prolif_subset.Rdata")

##### Genes Expression Prolif subset #####
DefaultAssay(Prolif_subset) <- "RNA"
Prolif_subset <- NormalizeData(Prolif_subset)

save(Prolif_subset, file = "Prolif_subset_normalizedRNA.Rdata")

all_markers_prolif <- FindAllMarkers(object = Prolif_subset,
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)

top10_prolif <- all_markers_prolif %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(all_markers_prolif, "Plots/all_markers_prolif_subset.csv" )
write.csv(top10_prolif, "Plots/top10_prolif_subsset.csv" )
save(all_markers_prolif, file = "Plots/all_markers_prolif_subset.RData")

#Figure3F
p <- FeaturePlot(object= Prolif_subset, features = c('Nes','Pax6'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes() 
p1 <- p[[1]] 
p2 <- p[[2]]
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure3_Feature_prolifsub_Nes_Pax6.png", p, base_width = 16, base_height = 8)

#Figure3Supp
p <- FeaturePlot(object= CombinedSCT, features = c('Nes','Pax6'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]] 
p2 <- p[[2]]
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure3Supp_Feature_Full_Nes_Pax6.png", p, base_width = 16, base_height = 8)

#Figure3G
p <- FeaturePlot(object= Prolif_subset, features = c('Ccnd1','Ccnd2'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]] 
p2 <- p[[2]]
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure3_Feature_prolifsub_Ccnd1_Ccnd2.png", p, base_width = 16, base_height = 8)

#Figure3Supp
p <- FeaturePlot(object= Prolif_subset, features = c('En1','En2'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]]
p2 <- p[[2]]
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure3Supp_Feature_prolifsub_rRL.png", p, base_width = 16, base_height = 8)

#Figure3Supp
p <- FeaturePlot(object= Prolif_subset, features = c('Hoxa2','Hoxa5'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]]
p2 <- p[[2]]
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure3Supp_Feature_prolifsub_cRL.png", p, base_width = 16, base_height = 8)

##Make the heatmap
#Scale data for input to heatmap
Prolif_subset <- ScaleData(object = Prolif_subset, features = rownames(x = Prolif_subset), verbose = TRUE)

#We will rotate the map 90 degrees so I flipped the order of the clusters so that they appear in descending order
my_levels = c('12','11','10','9','8','7','6','5','4','3','2','1','0')
levels(Prolif_subset) <- my_levels

#Get color palette information to then reverse the order of the colors to match the clusters above the heatmap
numclusters <- 13
my_color_palette <- hue_pal()(length(1:numclusters))

#Figure3E 
p <- DoHeatmap(Prolif_subset, assay = 'RNA', features = top10_prolif$gene, angle = 270, size = 12, draw.lines = TRUE, group.colors = rev(my_color_palette)) + theme(axis.text.y = element_text(size = 30, face = 'bold'), legend.position = 'bottom', legend.key.size = unit(1.5, "cm"), legend.title = element_text(size = 30), legend.text = element_text(size=30)) 
save_plot('Plots/Figure3_Prolif_Heatmap_top10.png', p, base_width = 20  , base_height =  40 )

### Differential expression on the progenitor and intermediate progenitor
#UMAP with just the progenitors and intermediate progenitors that are being compared.
Prog_comparison <- WhichCells(object = CombinedSCT, idents = c('0', '3', '29'))
Intermed_prog_comparison <- WhichCells(object = CombinedSCT, idents = c('35','7','8'))

p <- DimPlot(CombinedSCT, label=F, pt.size = 0.5, 
             cells.highlight = list(Prog_comparison, Intermed_prog_comparison), 
             cols.highlight = rev(c(Progenitor_cols, IntermediateProgenitor_cols)), cols = "gray") +
  scale_color_manual(labels = rev(c("Progenitor", "Intermediate Progenitor","Other")),
                     values = rev(c(Progenitor_cols, IntermediateProgenitor_cols, "gray"))) +
  guides(colour = guide_legend(override.aes = list(size=10), reverse = TRUE)) + NoAxes() +
  theme(legend.position = 'top', legend.text = element_text(size=22, face = 'bold'))
save_plot("Plots/Figure3J_FullData_UMAP_colored_progenitor_and_intermediateprog.png", p, base_width = 8, base_height = 8)

#Subset out the progenitor 
Idents(CombinedSCT) <- 'FullData.category'
Prog_and_Intermediate_subset <- subset(CombinedSCT, idents = c("Progenitor" ,"Intermediate Progenitor"))
#Subset out the cb rostral
Idents(object = Prog_and_Intermediate_subset ) <- "seurat_clusters"
Prog_and_Intermediate_subset <- subset(Prog_and_Intermediate_subset, idents = c("6","15","37"), invert = TRUE)
#Back to full.data category metadata 
Idents(object = Prog_and_Intermediate_subset ) <- "FullData.category"

#DEG expression to compare progenitors and intermediate progenitors
ProgenitorvIntermediateProg <- FindMarkers(CombinedSCT, ident.1 = "Progenitor", ident.2 = "Intermediate Progenitor", logfc.threshold = 0.25)
write.csv(ProgenitorvIntermediateProg, "Plots/ProgenitorvIntermediateProg.csv" )

#Pull out top and bottom 20 from DEG list and reorder by avg_logFC
top20_progen <- ProgenitorvIntermediateProg %>% top_n(n = 20, wt = avg_logFC) %>% arrange(desc(avg_logFC))
bottom20_progen <-ProgenitorvIntermediateProg %>% top_n(n = -20, wt = avg_logFC) %>% arrange(avg_logFC)
Combine_top_progn <- rbind(top20_progen, bottom20_progen)

#Scale the data for heatmap plotting
Prog_and_Intermediate_subset<- ScaleData(object = Prog_and_Intermediate_subset, features = rownames(x = Prog_and_Intermediate_subset), verbose = TRUE)

#Figure3K
p <- DoHeatmap(Prog_and_Intermediate_subset, assay = "RNA", features = row.names(Combine_top_progn), size = 8, angle = 0, hjust = 0.5, group.colors = c(Progenitor_cols, IntermediateProgenitor_cols)) +
  theme(axis.text.y = (element_text(size = 16, face = 'bold')), legend.position = "bottom")
save_plot("Plots/Figure3k_Prog v InterProg_heatmap.png", p, base_width = 11, base_height = 12)
save_plot("Plots/Figure3k_Prog v InterProg_heatmap.pdf", p, base_width = 11, base_height = 12)





