### Figure 7

#Set directory as you see fit
setwd("C:/Jessie/AllSamples/Publication/Figure7")

##Create Folders
dir.create("InitialQC")
dir.create("MinDist")
dir.create("Spread")
dir.create("MetaPlots_UMAP")
dir.create("Monocle")
dir.create("Plots")
dir.create("FindRes")

###Functions
source('C:/Jessie/Functions/find_res.R')
source('C:/Jessie/Functions/spread.R')
source('C:/Jessie/Functions/mindist.R')
source('C:/Jessie/Functions/MetaPlots_UMAP.R')
source("C:/Jessie/Functions/PlotGenePseudotime_UBC.R")

###Libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

###Dependencies

#Full Dataset Seurat Object
load("C:/Jessie/AllSamples/Publication/Figure1/FullDataset_normalizedRNA.Rdata")

#From the combined Seurat objects folder
#load("C:/Jessie/AllSamples/Publication/Seurat_Objects/FullDataset_normalizedRNA.Rdata")

#Prolif subset Seurat Object
load("C:/Jessie/AllSamples/Publication/Figure3/Prolif_subset_normalizedRNA.Rdata")

#From the combined Seurat objects folder
#load("C:/Jessie/AllSamples/Publication/Seurat_Objects/Prolif_subset_normalizedRNA.Rdata")

#Peak_List File
load("C:/Jessie/AllSamples/101221/Subset/CUT&RUN/Filtered_Peak_List.Rdata")
##Functions - MinDist, Spread, Monocle, tidyverse, colors, GetMultipleTraj, PlotGenesPseudotime_UBC 

###General Info
PlotTheme <-  theme(plot.title=element_text(size=50), legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=30)) 

#####Script#####
## Figure 7A - Feature plot of UBC progenitor markers in the proiferating subset
p <- FeaturePlot(object= Prolif_subset, features = c('Mki67','Eomes','Pax6','Otx2','Cbfa2t2','Cbfa2t3'), pt.size = 1, cols = firemap, combine = FALSE )
p1 <- p[[1]] + PlotTheme + NoAxes()
p2 <- p[[2]] + PlotTheme + NoAxes()
p3 <- p[[3]] + PlotTheme + NoAxes()
p4 <- p[[4]] + PlotTheme + NoAxes()
p5 <- p[[5]] + PlotTheme + NoAxes()
p6 <- p[[6]] + PlotTheme + NoAxes()
p <-  ggarrange(p1, p2, p3, p4, p5, p6, ncol = 6, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure7_prolifsubset_UBCmarkers.png", p, base_width = 48, base_height = 8, limitsize = FALSE)

rm(Prolif_subset)

#### Subset E16.5 rostral RL population to look at the proliferating UBCs

#stash old identity
CombinedSCT[["old.ident"]] <- Idents(object = CombinedSCT)

#subset the cb lineages
cb_subset <- subset(x = CombinedSCT, idents = c("14","15","23","20"))

png("InitialQC/VlnQCbySample.png", width = 2000, height = 400)
VlnPlot(object = cb_subset , features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.blood"), group.by = "orig.ident", cols = OrigIdentCols, ncol = 4, pt.size = 0)
dev.off()

cb_subset  <- SCTransform(cb_subset , vars.to.regress = c("percent.mito","percent.blood"), verbose = FALSE)

cb_subset  <- RunPCA(cb_subset , verbose = FALSE)
png("InitialQC/QC_PCHeatmap.png", width = 500, height = 1000)
DimHeatmap(cb_subset , dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png("InitialQC/QC_ElbowPlot.png", width = 500, height = 500)
ElbowPlot(cb_subset , ndims = 50)
dev.off()

cb_subset  <- RunUMAP(cb_subset , dims = 1:50, verbose = FALSE)
cb_subset  <- FindNeighbors(cb_subset , dims = 1:50, verbose = FALSE)

#Res at 50PC
find_res(cb_subset)

#Edit sample colors based on samples in subset
SampleCols <- brewer.pal(num_samples,"Paired")
SampleCols <- SampleCols[7:12]

#Find MinDist and Min Spread
MinDist(cb_subset, 50, 0.6, 1)
Spread(cb_subset, 50, 0.6, 0.1)

#Final Settings
cb_subset <- RunUMAP(cb_subset, spread = 3, min.dist = 0.1, dims = 1:50, verbose = FALSE)
cb_subset <- FindNeighbors(cb_subset, dims = 1:50, verbose = FALSE)
cb_subset <- FindClusters(cb_subset, verbose = FALSE, resolution= 0.6)

p <- DimPlot(cb_subset, reduction = "umap", label = TRUE, pt.size = 1, label.size = 16) + NoAxes()
save_plot("Plots/Figure7_UMAP.png", p, base_width = 16, base_height = 12)

MetaPlots_UMAP(cb_subset, 1)

save(cb_subset, file = "C:/Jessie/AllSamples/Publication/Figure7/cb_subset.Rdata")

#### RNA Expression ####

DefaultAssay(cb_subset) <- "RNA"
cb_subset <- NormalizeData(cb_subset)

save(cb_subset, file = "C:/Jessie/AllSamples/Publication/Figure7/cb_subset_normalizedRNA.Rdata")

all_markers_cb_subset <- FindAllMarkers(object = cb_subset,
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)

top20_cb_subset <- all_markers_cb_subset %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

write.csv(all_markers_cb_subset, "all_markers_cb_subset.csv" )
write.csv(top20_cb_subset, "Plots/top20_cb_subset.csv" )
save(all_markers_cb_subset, file = "all_markers_cb_subset.RData")

Cluster0VCluster6and9 <- FindMarkers(cb_subset, ident.1 = '0', ident.2 = c('6','9'), logfc.threshold = 0.25)
write.csv(Cluster0VCluster6and9, "cb_subset_Cluster0vCluster6and9.csv" )

#Figure 7D - FeaturePlot of UBC and EGL markers in the cb_subset
p <- FeaturePlot(object= cb_subset, features = c('Mki67','Pax6','Eomes','Otx2', 'Cdk1','Reln','Tlx3','Robo2'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]]
p4 <- p[[4]] 
p <-  ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure7_FeaturePlot_CbsubsetTrajGenes1.png", p, base_width = 32, base_height = 8)

p <- FeaturePlot(object= cb_subset, features = c('Cdk1','Reln','Tlx3','Robo2'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]]
p4 <- p[[4]] 
p <-  ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure7_FeaturePlot_CbsubsetTrajGenes2.png", p, base_width = 32, base_height = 8)


#Figure 7H - Feature plot of new markers of proliferating UBC in cb_subset
p <- FeaturePlot(object= cb_subset, features = c('Barx2','Npnt','Neurod4', 'Mgarp'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]]
p4 <- p[[4]]
p <-  ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure7_FeaturePlot_cbsubset_UBCProgGenes.png", p, base_width = 16, base_height = 16)

#Figure Supplement - Feature plot of new markers of proliferating UBC in full dataset
p <- FeaturePlot(object= CombinedSCT, features = c('Barx2','Npnt','Neurod4', 'Mgarp'), pt.size = 0.5, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]]
p4 <- p[[4]]
p <-  ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure7sup_FeaturePlot_full_UBCProgGenes.png", p, base_width = 16, base_height = 16)

##Do the new genes contain Atoh1 Peaks?
new_UBC_genes <- c('Barx2','Npnt','Neurod4', 'Mgarp')
UBCprog_genes_with_Peaks <- sum(new_UBC_genes %in% Peak_list)
Percent_UBCprog_genes_with_Peaks <- UBCprog_genes_with_Peaks/length(new_UBC_genes)*100

#Subset just the Unipolar Brush Cell lineage
UBC_subset <- subset(cb_subset, idents = c('13','0','1','2','12'))

#Pull colors from the original UMAP to match for the subset
num_clusters = 16
my_color_palette <- hue_pal()(num_clusters)

clusters_in_subset <- c(0,1,2,12,13)
cluster_pos <- 1+ clusters_in_subset
Eomes_color_palette <- my_color_palette[cluster_pos]

p <- DimPlot(UBC_subset, reduction = "umap", label = TRUE, pt.size = 1, label.size = 16, cols = Eomes_color_palette) + NoAxes() +
  theme(legend.position = 'none')
save_plot("Plots/Figure7sup_UMAP_UBC_subset.png", p, base_width = 14, base_height = 12)

save(UBC_subset, file = "C:/Jessie/AllSamples/Publication/Figure7/UBC_subset_normalizedRNA.Rdata")

#Subset just the cells that express Eomes
Eomes_UBC_subset <- subset(UBC_subset, subset = Eomes >1)

#Figure 7E UMAP of Eomes+ UBC
p <- DimPlot(Eomes_UBC_subset, reduction = "umap", label = TRUE, pt.size = 1, label.size = 16, cols = Eomes_color_palette) + NoAxes() +
     theme(legend.position = 'none')
save_plot("Plots/Figure7_UMAP_Eomes_subset.png", p, base_width = 14, base_height = 12)

save(Eomes_UBC_subset, file = "C:/Jessie/AllSamples/Publication/Figure7/Eomes_UBC_subset_normalizedRNA.Rdata")

##Monocal 

DefaultAssay(Eomes_UBC_subset) <- "SCT"
#This is generally what I followed; http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html
cds <- as.cell_data_set(Eomes_UBC_subset)

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_Eomes_UBC.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

#Learn Graph
cds <- learn_graph(cds, use_partition = FALSE)
png("Monocle/Monocle_Trajectory&Partition_trajectory_label_Eomes_UBC.png", width = 500, height = 500)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE )
dev.off()

#order Cells
cds <- order_cells(cds)
#Figure 7F - Eomes+ UBC in pseudotime
png("Plots/Figure7_Monocle_pseudotime_Eomes_subset.png", width = 600, height = 500)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE) + NoAxes() + 
           theme(legend.key.size = unit(1.0, "cm"), legend.title = element_text(size=20), legend.text = element_text(size = 20))
dev.off()

##Add pseudotime values into the Seurat object and plot gene expression in pseudotime
pseudotime_Eomes_UBC <- as.data.frame(pseudotime(cds))
Eomes_UBC_subset <- AddMetaData(Eomes_UBC_subset, pseudotime_Eomes_UBC, col.name = 'Eomes_UBC_pseudotime')

###Plot genes in pseudotime
UBC_genes <- c('Atoh1','Sox2','Mki67','Eomes','Otx2','Cbfa2t2','Lmx1a','Mapt')
traj.labels <- 'Eomes_UBC_pseudotime'

traj.dat <- GetMultipleTraj(Eomes_UBC_subset, traj.labels,UBC_genes)

#Figure 7G - change in gene expression of Eomes+ UCB in pseudotime
p<- PlotGenePseudotime_UBC(traj.dat)
save_plot('Plots/Figure7_UBC_genes_Lineplot.pdf', p, base_width = 4, base_height = 6)

