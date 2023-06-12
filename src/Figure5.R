##Figure5

#Set directory as you see fit
setwd("C:/Jessie/AllSamples/Publication/Figure5")

##Create Folders
dir.create("Plots")
dir.create("Monocle")

##Functions
source('C:/Jessie/Functions/RenameClusterGroups.R')
source('C:/Jessie/Functions/MakeTimePlot.R')
source("C:/Jessie/Functions/PlotGenePseudotime.R")

## Libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggpubr)
library(data.table)
library(ggplot2)
library(monocle3)

###Dependencies
#Full Dataset Seurat Object
load("C:/Jessie/AllSamples/Publication/Figure1/FullDataset_normalizedRNA.Rdata")

#From the combined Seurat objects folder
#load("C:/Jessie/AllSamples/Publication/Seurat_Objects/FullDataset_normalizedRNA.Rdata")

#DEG list from Full Dataset
load("C:/Jessie/AllSamples/Publication/Figure1/Plots/all_markers_FullDataset.RData")

#Label of migrating clusters
Migration_IDs <- "C:/Jessie/AllSamples/Publication/RenameCluster_Files/Migration_IDs.txt"

#Peak_List File
load("C:/Jessie/AllSamples/101221/Subset/CUT&RUN/Filtered_Peak_List.Rdata")

#####Script#####

Idents(CombinedSCT) <- 'seurat_clusters'

#Add in IDs of Migration streams
CombinedSCT <- RenameClusterGroups(CombinedSCT, Migration_IDs)

Migration.levels <- c("RLS", "CB_nuclei","EGL1", "EGL2", "Early_Pons","CLS", "Caudal_Pons", "PES", "CES", "AES_early", "AES_late", "Other")
Migration_cols =c(Early_RLS_col, CBNuc_cols, EGL1_cols, EGL2_cols, Early_pons_col, CLS_col, Caudal_pons_col, PES_col, CES_col,AES_col, AES_late_col, 'gray')

CombinedSCT[["Migration.category"]] <- Idents(object = CombinedSCT)
CombinedSCT$Migration.category <- factor(CombinedSCT$Migration.category, levels = Migration.levels)

#Figure5B
p <- DimPlot(CombinedSCT, group.by = 'Migration.category',pt.size = 0.5 ,cols = Migration_cols) + NoAxes() + 
  guides(colour = guide_legend(override.aes = list(size=14))) +
  theme(legend.text = element_text(size = 35, face = 'bold')) 
save_plot("Plots/Figure5_UMAP_by_MigrationStream.png", p, base_width = 12, base_height = 8) 

##Make dataframe with sum of cells in different migration streams across samples
#Calculate Number of cells in each sample by migration ID
n_cells <- FetchData(CombinedSCT, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

#Remove all Td only samples because these are too mature
n_cells <- n_cells %>% 
  filter(!grepl('_Td', orig.ident))

#remove the '_rep2' label from second replicate
n_cells$orig.ident <- str_remove(n_cells$orig.ident, pattern = "_rep2")

#Add together rows that have the same label (sum the cells from the two replicates)
n_cells <- n_cells %>%
  group_by(orig.ident) %>%
  summarise_all(list(sum))

#reorganize the dataframe
MigrationStream <- n_cells %>% pivot_longer(!orig.ident, names_to = 'variable', values_to = 'value')

#Remove information from non-migration stream clusters
MigrationStream <- MigrationStream %>% 
  filter(!grepl('Other', variable))

#Replace NA with 0
MigrationStream <- MigrationStream %>% replace(is.na(.), 0)

#Figure5C
p <- MakeTimePlot(MigrationStream)
save_plot("Plots/Figure5_MigrationStream_LinePlot.png", p, base_width = 8, base_height = 8)
save_plot("Plots/Figure5_MigrationStream_LinePlot.pdf", p, base_width = 8, base_height = 8)

##### Pseudotime ####

#Find pseudotime for each migration stream individually 

#### Cluster 2 ####
#Cluster2 and its intermediate prog. cluster8
Cluster2 <- subset(x = CombinedSCT, idents = c('8','2'))
DefaultAssay(Cluster2) <- "SCT"

cds <- as.cell_data_set(Cluster2)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_cluster2.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

#Learn Graph
cds <- learn_graph(cds, use_partition = FALSE)
png("Monocle/Monocle_Trajectory&Partition_trajectory_label_cluster2.png", width = 500, height = 500)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE )
dev.off()

#order Cells
cds <- order_cells(cds)
png("Monocle/Monocle_pseudotime_cluster2_compare.png", width = 600, height = 500)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
dev.off()

pseudotime_cluster2 <- as.data.frame(pseudotime(cds))

CombinedSCT <- AddMetaData(CombinedSCT, pseudotime_cluster2, col.name = 'cluster2_pseudotime')
rm(Cluster2)
rm(cds)

#### Cluster 5 ####
Cluster5 <- subset(x = CombinedSCT, idents = c('5'))

DefaultAssay(Cluster5) <- "SCT"
cds <- as.cell_data_set(Cluster5)

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_cluster5.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

#Learn Graph
cds <- learn_graph(cds, use_partition = FALSE)
png("Monocle/Monocle_Trajectory&Partition_trajectory_label_cluster5.png", width = 500, height = 500)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE )
dev.off()

#Order cells
cds <- order_cells(cds)
png("Monocle/Monocle_pseudotime_cluster5.png", width = 600, height = 500)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
dev.off()

pseudotime_cluster5 <- as.data.frame(pseudotime(cds))
CombinedSCT <- AddMetaData(CombinedSCT, pseudotime_cluster5, col.name = 'cluster5_pseudotime')
rm(Cluster5)
rm(cds)

#### Cluster 10 ####
Cluster10 <- subset(x = CombinedSCT, idents = c('7','10'))

DefaultAssay(Cluster10) <- "SCT"
cds <- as.cell_data_set(Cluster10)

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_cluster10.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

#Learn Graph
cds <- learn_graph(cds, use_partition = FALSE)
png("Monocle/Monocle_Trajectory&Partition_trajectory_label_cluster10.png", width = 500, height = 500)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE )
dev.off()

#Order Cells
cds <- order_cells(cds)
png("Monocle/Monocle_pseudotime_cluster10.png", width = 600, height = 500)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
dev.off()

pseudotime_cluster10 <- as.data.frame(pseudotime(cds))
CombinedSCT <- AddMetaData(CombinedSCT, pseudotime_cluster10, col.name = 'cluster10_pseudotime')

rm(Cluster10)
rm(cds)

#### Cluster 13 ####
Cluster13 <- subset(x = CombinedSCT, idents = c('7','8','13'))

DefaultAssay(Cluster13) <- "SCT"
cds <- as.cell_data_set(Cluster13)

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_cluster13.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

##Subset out part of progenitor cluster
cds.sub <- subset(as.Seurat(cds), idents = c('2','3','4','5','6'))
cds <- as.cell_data_set(cds.sub)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(CombinedSCT[["SCT"]])

#Regraph
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_subset_cluster13.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

#Learn Graph
cds <- learn_graph(cds, use_partition = FALSE)
png("Monocle/Monocle_Trajectory&Partition_trajectory_label_cluster13.png", width = 500, height = 500)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE )
dev.off()

#Pick where to start
cds <- order_cells(cds)
png("Monocle/Monocle_pseudotime_cluster13.png", width = 600, height = 500)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
dev.off()

pseudotime_cluster13 <- as.data.frame(pseudotime(cds))
CombinedSCT <- AddMetaData(CombinedSCT, pseudotime_cluster13, col.name = 'cluster13_pseudotime')

rm(Cluster13)
rm(cds)

#### Cluster 17 ####
Cluster17 <- subset(x = CombinedSCT, idents = c('8','17'))

DefaultAssay(Cluster17) <- "SCT"
cds <- as.cell_data_set(Cluster17)

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_cluster17.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

##Subset out part of progenitor cluster
cds.sub <- subset(as.Seurat(cds), idents = c('2','3','4'))
cds <- as.cell_data_set(cds.sub)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(CombinedSCT[["SCT"]])

#Regraph
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_subset_cluster17.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

#Learn Graph
cds <- learn_graph(cds, use_partition = FALSE)
png("Monocle/Monocle_Trajectory&Partition_trajectory_label_cluster17.png", width = 500, height = 500)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE )
dev.off()

#Pick where to start
cds <- order_cells(cds)
png("Monocle/Monocle_pseudotime_cluster17.png", width = 600, height = 500)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
dev.off()

pseudotime_cluster17 <- as.data.frame(pseudotime(cds))
CombinedSCT <- AddMetaData(CombinedSCT, pseudotime_cluster17, col.name = 'cluster17_pseudotime')

rm(Cluster17)
rm(cds)

### Cluster 21 ####
Cluster21 <- subset(x = CombinedSCT, idents = c('35','21'))

DefaultAssay(Cluster21) <- "SCT"
cds <- as.cell_data_set(Cluster21)

cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
png("Monocle/Monocle_Full_Trajectory&Partition_cluster21.png", width = 1000, height = 500)
wrap_plots(p1, p2)
dev.off()

#Learn Graph
cds <- learn_graph(cds, use_partition = FALSE)
png("Monocle/Monocle_Trajectory&Partition_trajectory_label_cluster21.png", width = 500, height = 500)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE )
dev.off()

#Pick where to start
cds <- order_cells(cds)
png("Monocle/Monocle_pseudotime_cluster21.png", width = 600, height = 500)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
dev.off()

pseudotime_cluster21 <- as.data.frame(pseudotime(cds))
CombinedSCT <- AddMetaData(CombinedSCT, pseudotime_cluster21, col.name = 'cluster21_pseudotime')

rm(Cluster21)
rm(cds)

save(CombinedSCT, file = "FullDataset_pseudotime.Rdata")

##### 
# Plot pseudotime back in full UMAP

PlotTheme <- theme(title = element_text(size = 50), legend.key.size = unit(2.5, "cm"), legend.text = element_text(size=50))

p <- FeaturePlot(CombinedSCT, pt.size = 0.5, features = c('cluster2_pseudotime','cluster5_pseudotime', 'cluster10_pseudotime','cluster13_pseudotime','cluster17_pseudotime','cluster21_pseudotime'), combine = FALSE) 
p1 <- p[[1]] + NoAxes() + ggtitle("Early RLS") + PlotTheme
p2 <- p[[6]] + NoAxes() + ggtitle("Early Pons") + PlotTheme
p3 <- p[[3]] + NoAxes() + ggtitle("CLS") + PlotTheme
p4 <- p[[4]] + NoAxes() + ggtitle("Caudal Pons") + PlotTheme
p5 <- p[[5]] + NoAxes() + ggtitle("CES") + PlotTheme
p6 <- p[[2]] + NoAxes() + ggtitle("AES Early") + PlotTheme
p <-  ggarrange(p1, p2,p3, p4, p5, p6, common.legend = TRUE, legend = 'right') 
save_plot(p, file = 'Plots/Pseudotime_compiled.png', base_width = 24, base_height = 16)

Migrating_Clusters <- c("2", "21", "13", "10", "17", "5")
#To find shared genes:
all_markers_SCT_filter <- filter(all_markers_FUllDataset, cluster %in% Migrating_Clusters)
genes_common_to_all <- all_markers_SCT_filter %>% group_by(gene) %>% filter(n()>5)
genes_common_to_all <- genes_common_to_all %>%  filter (!duplicated(gene))
genes_common_to_all <- genes_common_to_all$gene
write.csv(genes_common_to_all, "Plots/Figure5_genes_common_to_all.csv")
#[1] "Gadd45g" "Sstr2"   "Mfng"    "Srrm4"   "Atoh1"   "Rassf4"  "Nhlh1"   "Lzts1"   "Dll3"    "Akap6"   "Gse1"    "Cntn2"  
#[13] "Fam78b" 

genes_common_to_5orMore <- all_markers_SCT_filter %>% group_by(gene) %>% filter(n()>4)
genes_common_to_5orMore <- genes_common_to_5orMore %>%  filter (!duplicated(gene))
genes_common_to_5orMore <- genes_common_to_5orMore$gene
write.csv(genes_common_to_5orMore, "Plots/Figure5_genes_common_to_5orMore.csv")
#[1] "Gadd45g" "Insm1"   "Nkd1"    "Ebf2"    "Igsf8"   "Sstr2"   "Mfng"    "Srrm4"   "Atoh1"   "Rassf4"  "Nhlh1"   "Selenom"
#[13] "Cbfa2t2" "Lzts1"   "Dll3"    "Rnd2"    "Mtcl1"   "Zeb1"    "Ube2ql1" "Abracl"  "Epha5"   "Ebf3"    "Akap6"   "Bcl7a"  
#[25] "Barhl1"  "Gse1"    "Nav2"    "Tubb3"   "Elavl4"  "Cntn2"   "Fam78b"  "Vat1"    "Ier2"    "Btbd17"  "Igfbp2"  "Hoxa2"  
#[37] "Eya2"    "Enox2" 

##how many of these genes are in the Peak_List defined in Figure 4 
sum_genes_with_peaks <- sum(genes_common_to_5orMore %in% Peak_list)
percent_of_genes_with_peaks <- sum_genes_with_peaks/length(genes_common_to_5orMore)*100
#[1] 92.10526

##To find unique genes
Idents(CombinedSCT) <- 'seurat_clusters'

###Differential gene expression to find unique genes for each migration stream
RLSvMig <- FindMarkers(CombinedSCT, ident.1 = "2", ident.2 = c("5","10", "13", "17", "21"), logfc.threshold = 0.25, only.pos = T)
write.csv(RLSvMig, "Plots/Figure5_Fullset_RLSvMig.csv" )

AES_earlyvMig <- FindMarkers(CombinedSCT, ident.1 = "5", ident.2 = c("2","10", "13", "17", "21"), logfc.threshold = 0.25, only.pos = T)
write.csv(AES_earlyvMig, "Plots/FIgure5_Fullset_AES_earlyvMig.csv" )

CLSvMig <- FindMarkers(CombinedSCT, ident.1 = "10", ident.2 = c("2","5", "13", "17", "21"), logfc.threshold = 0.25, only.pos = T)
write.csv(CLSvMig, "Plots/Figure5_Fullset_CLSvMig.csv" )

Caudal_ponsvMig <- FindMarkers(CombinedSCT, ident.1 = "13", ident.2 = c("2","5", "10", "17", "21"), logfc.threshold = 0.25, only.pos = T)
write.csv(Caudal_ponsvMig, "Plots/Figure5_Fullset_Caudal_ponsvMig.csv" )

CESvMig <- FindMarkers(CombinedSCT, ident.1 = "17", ident.2 = c("2","5", "10", "13", "21"), logfc.threshold = 0.25, only.pos = T)
write.csv(CESvMig, "Plots/Figure5_Fullset_CESvMig.csv" )

Early_ponsvMig <- FindMarkers(CombinedSCT, ident.1 = "21", ident.2 = c("2","5", "10", "13","17"), logfc.threshold = 0.25, only.pos = T)
write.csv(Early_ponsvMig, "Plots/Figure5_Fullset_Early_ponsvMig.csv" )

### Plot gene changes in pseudotime
Shared_Genes_short <- c("Atoh1","Dll3", "Mfng", "Lzts1","Gadd45g","Cntn2")
Unique_Genes_short <- c("Pxylp1","Lgals1","Pou4f1","C1ql3","Pvalb", "S100a6")
CellState_Genes <- c("Mki67", "Sox2", "Mapt")

traj.cols <- c(Early_RLS_col , Early_pons_col, CLS_col, Caudal_pons_col, CES_col, AES_col)
traj.labels <- c("cluster2_pseudotime","cluster21_pseudotime", "cluster10_pseudotime", "cluster13_pseudotime", "cluster17_pseudotime", "cluster5_pseudotime")

traj.dat.unique <- GetMultipleTraj(CombinedSCT, traj.labels, Unique_Genes_short)
traj.dat.unique$trajectory <- factor(traj.dat.unique$trajectory, level = traj.labels)
p <- PlotGenePseudotime(traj.dat.unique)
p1 <- p + scale_color_manual(labels = c('Early RLS', 'Early Pons', 'CLS', 'Caudal Pons', 'CES', 'AES_early'), values = traj.cols) 

traj.dat.shared <- GetMultipleTraj(CombinedSCT, traj.labels, Shared_Genes_short)
traj.dat.shared$trajectory <- factor(traj.dat.shared$trajectory, level = traj.labels)
p <- PlotGenePseudotime(traj.dat.shared)
p2 <- p + scale_color_manual(labels =c('Early RLS', 'Early Pons', 'CLS', 'Caudal Pons', 'CES', 'AES_early'), values = traj.cols )

#Figure5F Combine unique and shared
p <- ggarrange(p2, p1, common.legend = TRUE, legend = 'right')
save_plot('Plots/Figure5_Combined_Migration_PseudotimeLineplot.png', p, base_width = 8, base_height = 6)
save_plot('Plots/Figure5_Combined_Migration_PseudotimeLineplot.pdf', p, base_width = 8, base_height = 6)

#Figure5E
traj.dat.cellstate <- GetMultipleTraj(CombinedSCT, traj.labels, CellState_Genes)
traj.dat.cellstate$trajectory <- factor(traj.dat.cellstate$trajectory, level = traj.labels)
p <- PlotGenePseudotime(traj.dat.cellstate)
p <- p + scale_color_manual(labels = c('Early RLS', 'Early Pons', 'CLS', 'Caudal Pons', 'CES', 'AES_early'), values = traj.cols )
save_plot('Plots/Figure5_Cellstate_Migration_PseudotimeLinePlot.png', p, base_width = 4, base_height = 3)
save_plot('Plots/FIgure5_Cellstate_Migration_PseudotimeLinePlot.pdf', p, base_width = 4, base_height = 3)



