### Figure 6

#Set directory as you see fit
setwd("C:/Jessie/AllSamples/Publication/Figure6")

###Create Folders
dir.create("Subset1")
dir.create("Subset2")
dir.create("Plots")

###Functions
source('C:/Jessie/Functions/find_res.R')
source('C:/Jessie/Functions/spread.R')
source('C:/Jessie/Functions/mindist.R')
source('C:/Jessie/Functions/MetaPlots_UMAP.R')
source('C:/Jessie/Functions/Atoh1_Identity.R')
source('C:/Jessie/Functions/RenameClusterGroups.R')

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

#Mature subset seurat object
load("C:/Jessie/AllSamples/Publication/Figure2/Mature_subset_normalizedRNA.Rdata")

#From the combined Seurat objects folder
#load("C:/Jessie/AllSamples/Publication/Seurat_Objects/Mature_subset_normalizedRNA.Rdata")

#Peak_List File
load("C:/Jessie/AllSamples/101221/Subset/CUT&RUN/Filtered_Peak_List.Rdata")

#ID List
caudal_pons_ID <- "C:/Jessie/AllSamples/Publication/RenameCluster_Files/CaudalPons_IDs.txt"
Mature_IDs_coarse <- "C:/Jessie/AllSamples/Publication/RenameCluster_Files/Mature_IDs_coarse.txt"

###General Info
PlotTheme <-  theme(plot.title=element_text(size=50), legend.key.size = unit(0.75, "cm"), legend.text = element_text(size=30)) 

#####Script#####
setwd("C:/Jessie/AllSamples/Publication/Figure6/Subset1")
dir.create("InitialQC")
dir.create("MinDist")
dir.create("Spread")
dir.create("FindRes")
dir.create("MetaPlots_UMAP")
dir.create("Plots")

#stash old identity# 
CombinedSCT[["old.ident"]] <- Idents(object = CombinedSCT)

#Including both of the intermediate progenitors so nothing is missed. Taking whole "caudal pons" stream - 7,8,13,24,33,34
caudalpons_initial_subset <- subset(x = CombinedSCT, idents = c("7","8","13","24","33","34"))

png("InitialQC/VlnQCbySample.png", width = 2000, height = 400)
VlnPlot(object = caudalpons_initial_subset , features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.blood"), group.by = "orig.ident", cols = OrigIdentCols, ncol = 4, pt.size = 0)
dev.off()

caudalpons_initial_subset  <- SCTransform(caudalpons_initial_subset , vars.to.regress = c("percent.mito","percent.blood"), verbose = FALSE)

caudalpons_initial_subset  <- RunPCA(caudalpons_initial_subset , verbose = FALSE)
png("InitialQC/QC_PCHeatmap.png", width = 500, height = 1000)
DimHeatmap(caudalpons_initial_subset , dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png("InitialQC/QC_ElbowPlot.png", width = 500, height = 500)
ElbowPlot(caudalpons_initial_subset , ndims = 50)
dev.off()

caudalpons_initial_subset  <- RunUMAP(caudalpons_initial_subset , dims = 1:50, verbose = FALSE)
caudalpons_initial_subset  <- FindNeighbors(caudalpons_initial_subset , dims = 1:50, verbose = FALSE)

#Res at 50PC
find_res(caudalpons_initial_subset)

caudalpons_initial_subset <- RunUMAP(caudalpons_initial_subset, dims = 1:50, verbose = FALSE)
caudalpons_initial_subset <- FindNeighbors(caudalpons_initial_subset, dims = 1:50, verbose = FALSE)
caudalpons_initial_subset <- FindClusters(caudalpons_initial_subset, verbose = FALSE, resolution= 0.4)

MetaPlots_UMAP(caudalpons_initial_subset, 1)

p <- DimPlot(caudalpons_initial_subset, reduction = "umap", label = TRUE, pt.size = 1, label.size = 16) + NoAxes()
save_plot("Plots/Figure6_sup_Caudalpons_initial_UMAP.png", p, base_width = 16, base_height = 12)

save(caudalpons_initial_subset, file = "caudalpons_initial_subset_50PC_04res.Rdata")

#### RNA Expression ####
DefaultAssay(caudalpons_initial_subset) <- "RNA"
caudalpons_initial_subset <- NormalizeData(caudalpons_initial_subset)

Atoh1_Identity(caudalpons_initial_subset, 0.5)

save(caudalpons_initial_subset, file = "caudalpons_initial_subset_normalized RNA.Rdata")

#####
##Subset again to clean up 
setwd("C:/Jessie/AllSamples/Publication/Figure6/Subset2")
dir.create("InitialQC")
dir.create("MinDist")
dir.create("Spread")
dir.create("FindRes")
dir.create("MetaPlots_UMAP")

caudalpons_initial_subset[["old.ident2"]] <- Idents(object = caudalpons_initial_subset)

#To pick this subset, I looked at the 3D projection of the initial subset to pick the best progenitor populations. 
#I also looked at Hoxa5 expression and saw that cluster 12 had high Hoxa5 expression and may be a population coming from the CLS so I took it out 
caudalpons_subset2 <- subset(x = caudalpons_initial_subset, idents = c("1",'2',"7","8","9", "10","11","13","15"))

png("InitialQC/VlnQCbySample.png", width = 2000, height = 400)
VlnPlot(object = caudalpons_subset2 , features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.blood"), group.by = "orig.ident", cols = OrigIdentCols, ncol = 4, pt.size = 0)
dev.off()

caudalpons_subset2  <- SCTransform(caudalpons_subset2 , vars.to.regress = c("percent.mito","percent.blood"), verbose = FALSE)

caudalpons_subset2  <- RunPCA(caudalpons_subset2 , verbose = FALSE)
png("InitialQC/QC_PCHeatmap.png", width = 500, height = 1000)
DimHeatmap(caudalpons_subset2 , dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png("InitialQC/QC_ElbowPlot.png", width = 500, height = 500)
ElbowPlot(caudalpons_subset2 , ndims = 50)
dev.off()

caudalpons_subset2  <- RunUMAP(caudalpons_subset2 , dims = 1:50, verbose = FALSE)
caudalpons_subset2  <- FindNeighbors(caudalpons_subset2 , dims = 1:50, verbose = FALSE)

#Res at 50PC
SampleCols <- brewer.pal(num_samples,"Paired")
SampleCols <- SampleCols[2:12]

find_res(caudalpons_subset2)
MinDist(caudalpons_subset2, 50, 0.4, 1)
Spread(caudalpons_subset2, 50, 0.4, 0.15)

#Set res you want to keep
caudalpons_subset2 <- RunUMAP(caudalpons_subset2, dims = 1:50, spread = 0.5, min.dist = 0.15, verbose = FALSE)
caudalpons_subset2 <- FindNeighbors(caudalpons_subset2, dims = 1:50, verbose = FALSE)
caudalpons_subset2 <- FindClusters(caudalpons_subset2, verbose = FALSE, resolution= 0.3)

MetaPlots_UMAP(caudalpons_subset2, 1)

### Generation of plots for publication
setwd("C:/Jessie/AllSamples/Publication/Figure6")
UMAPPlotTheme <- theme(axis.title =element_text(size=30), axis.text=element_text(size=30),legend.text=element_text(size=25, face = 'bold')) + NoAxes()

#Figure6b
p <- DimPlot(caudalpons_subset2, pt.size = 1, label = TRUE, label.size = 16) + NoAxes() + scale_x_reverse() +
  theme(legend.position = 'none')
save_plot("Plots/Figure6_caudalpons_subset2_UMAP.png", p, base_width = 13, base_height = 12)

#Figure6c
p <- DimPlot(caudalpons_subset2, group.by = "Sample", pt.size = 1, cols =SampleCols, label = F) + 
  guides(colour = guide_legend(override.aes = list(size=14))) +
  theme(legend.text = element_text(size = 50, face = 'bold')) + NoAxes() + scale_x_reverse()
save_plot("Plots/Figure6_caudalpons_subset2_Sample.png", p, base_width = 18, base_height = 12)

#Figure 6 Supp
p <- DimPlot(object = caudalpons_subset2, group.by="orig.ident", cols = OrigIdentCols2, pt.size = 1, label.size = 3) + scale_x_reverse() +  
  UMAPPlotTheme + guides(colour = guide_legend(override.aes = list(size=12)))
save_plot("Plots/Figure6Supp_UMAP_byorigident.png", p, base_width = 16, base_height = 8)

p <- DimPlot(object = caudalpons_subset2, group.by="Replicate", pt.size = 1, label.size = 3) + scale_x_reverse() +
  UMAPPlotTheme + guides(colour = guide_legend(override.aes = list(size=12)))
save_plot("Plots/FIgure6Supp_UMAP_byReplicate.png", p, base_width = 8, base_height = 8)

###Label caudal pons cells back on FullDataset
##Build data frame of Cell ID and cluster ID from caudal pons
Idents(caudalpons_subset2) <- 'Caudal_Pons_subset'
caudalpons_subset2[["Caudal_Pons"]] <- Idents(object = caudalpons_subset2)

caudalpons_subset2_df <- data.frame(Cells(caudalpons_subset2),caudalpons_subset2$Caudal_Pons)
colnames(caudalpons_subset2_df) <- c('CellID','Caudal_Pons')

CombinedSCT_Cells <- data.frame(Cells(CombinedSCT))
colnames(CombinedSCT_Cells) <- 'CellID'

# Make dataframe that Merges Labels onto the list of ALL cell IDs
combined_labels <- CombinedSCT_Cells %>% full_join(caudalpons_subset2_df)

#Add metadata into full dataset
CombinedSCT <- AddMetaData(CombinedSCT, combined_labels$Caudal_Pons, col.name = 'CaudalPons_Subset')

#Figure6A
p <- DimPlot(CombinedSCT, group.by = 'CaudalPons_Subset', pt.size = 1, cols = Caudal_pons_col) + NoAxes() +
  theme(legend.position = 'none')
save_plot("Plots/Figure6_CaudalPons_FullDataSet.png", p, base_width = 14, base_height = 12)

rm(CombinedSCT)
Idents(caudalpons_subset2) <- 'seurat_clusters'

###Label mature cells on to caudal pons

#Relabel Mature_subset with Mature IDs_coarse
Idents(Mature_subset) <- 'seurat_clusters'
Mature_subset <- RenameClusterGroups(Mature_subset, Mature_IDs_coarse)
Mature_subset[["Cluster.Labels"]] <- Idents(object = Mature_subset)

##Build data frame of Cell ID and cluster ID from Mature Subset
mature_nuclei_df <- data.frame(Cells(Mature_subset), Mature_subset$Cluster.Labels)
colnames(mature_nuclei_df) <- c('CellID','MatureClusterID')

CaudalPons_Cells <- data.frame(Cells(caudalpons_subset2))
colnames(CaudalPons_Cells) <- 'CellID'

# Make dataframe that Merges Labels onto the list of ALL cell IDs
mature_combined_labels <- CaudalPons_Cells %>% full_join(mature_nuclei_df)
mature_combined_labels <- mature_combined_labels %>% filter(CellID %in% CaudalPons_Cells$CellID)

#Add metadata into full dataset
caudalpons_subset2 <- AddMetaData(caudalpons_subset2, mature_combined_labels$MatureClusterID, col.name = 'MatureID')

p <- DimPlot(caudalpons_subset2, group.by = 'MatureID', pt.size = 0.5) +
  guides(colour = guide_legend(override.aes = list(size=14))) +
  theme(legend.text = element_text(size = 30, face = 'bold')) + NoAxes() + scale_x_reverse()
save_plot("Plots/Figure6Supp_caudalpons_MatureLabels.png", p, base_width = 24, base_height = 12)

n_cells <- FetchData(caudalpons_subset2, vars = c("ident", "MatureID")) %>%
  dplyr::count(ident, MatureID) %>%
  tidyr::spread(ident, n)
write.csv(n_cells, "Plots/MatureID_onCaudalPons.csv" )

rm(Mature_subset)

save(caudalpons_subset2, file = "CaudalPons_subset2.Rdata")

### RNA Expression ####
DefaultAssay(caudalpons_subset2) <- "RNA"
caudalpons_subset2 <- NormalizeData(caudalpons_subset2)

save(caudalpons_subset2, file = "CaudalPons_subset2_normalized RNA.Rdata")

all_markers_caudalPons_subset2 <- FindAllMarkers(object = caudalpons_subset2,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

top20_caudalPons_subset2 <- all_markers_caudalPons_subset2 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

write.csv(all_markers_caudalPons_subset2, "Plots/all_markers_caudalPons_subset2.csv" )
write.csv(top20_caudalPons_subset2, "Plots/top20_caudalPons_subset2.csv" )
save(all_markers_caudalPons_subset2, file = "Plots/all_markers_caudalPons_subset2.RData")

#Figure6D
p  <- FeaturePlot(object= caudalpons_subset2, features = c('Nes','Sox2','Atoh1','Nhlh1','Vcan', 'C1ql3', 'Slc17a6', 'Mapt'), pt.size = 0.5, cols = firemap ) & PlotTheme & NoAxes() & scale_x_reverse()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]] 
p4 <- p[[4]] 
p5 <- p[[5]] 
p6 <- p[[6]] 
p7 <- p[[7]]
p8 <- p[[8]]
p <-  ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 8, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure6_FeaturePlot_NeuronalMaturation.png", p, base_width = 32, base_height = 4) 

#Figure6E Feature Plots of final mature populations
p  <- FeaturePlot(object= caudalpons_subset2, features = c('Tcf24','Sox7','Kcnip2',"Sncg", "Mafb"),pt.size = 0.5, cols = firemap ) & NoAxes() & scale_x_reverse() &
  theme(plot.title=element_text(size=35), legend.key.size = unit(0.5, "cm"), legend.text = element_text(size=20))
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]] 
p4 <- p[[4]] 
p5 <- p[[5]] 
p <-  ggarrange(p1, p2, p3, p4, p5, ncol = 5, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figur6_FeaturePlot_FinalPop.png", p, base_width = 20, base_height = 4)

#Figure6 supplemental 
p <- FeaturePlot(object= caudalpons_subset2, features = c('En1','En2'), pt.size = 0.5, cols = firemap) & PlotTheme & NoAxes() & scale_x_reverse()
p1 <- p[[1]]
p2 <- p[[2]]
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure6supp_Feature_caudalponssub_rRL.png", p, base_width = 8, base_height = 4)

#Figure6 supplemental
p <- FeaturePlot(object= caudalpons_subset2, features = c('Hoxa2','Hoxa5'), pt.size = 0.5, cols = firemap) & PlotTheme & NoAxes() & scale_x_reverse()
p1 <- p[[1]]
p2 <- p[[2]]
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure6supp_Feature_caudalponssub_cRL.png", p, base_width = 8, base_height = 4)


#Label final populations
Pr5 <- WhichCells(object = caudalpons_subset2, idents = '5')
ILL <- WhichCells(object = caudalpons_subset2, idents = '10')
SuVe <- WhichCells(object = caudalpons_subset2, idents = '11')
SupOlive <- WhichCells(object = caudalpons_subset2, idents = '9')

#Figure6O
p <- DimPlot(caudalpons_subset2, label=F, pt.size = 0.5, 
             cells.highlight = list(Pr5,ILL,SuVe,SupOlive)) + 
  scale_color_manual(labels = rev(c("Pr5", "ILL","SuVe","SON/Ve","Other")),
                     values = rev(c('darkgreen', 'firebrick3', 'magenta2', 'purple4', 'gray'))) +
  guides(colour = guide_legend(override.aes = list(size=16), reverse = T)) +
  theme(legend.text = element_text(size = 55, face = 'bold')) +
  NoAxes() + scale_x_reverse()
save_plot("Plots/Figure6_UMAP_FinalPopLabels.png", p, base_width = 16, base_height = 12)

#Rename clusters by lineage
caudalpons_subset2 <- RenameClusterGroups(caudalpons_subset2, caudal_pons_ID)

#Figure6P
p <- DimPlot(caudalpons_subset2, pt.size = 0.5, label = FALSE, 
             order = rev(c('Progenitor', 'Lineage1', 'Lineage2', 'Lineage3', 'NA'))) + 
  scale_color_manual(labels = c("Progenitor", "Lineage_1","Lineage_2","Lineage_3","Other"),
                     values = c('#00BFC4', '#7CAE00', '#F8766D', '#C77CFF')) +
  guides(colour = guide_legend(override.aes = list(size=16))) +
  theme(axis.text=element_text(size=16), legend.text = element_text(size = 55, face = 'bold')) + NoAxes() + scale_x_reverse()
save_plot("Plots/Figure6_UMAP_LineageIDLabels.png", p, base_width = 17, base_height = 12)

##Differential gene expression on the lineage identities
all_markers_LineageID <- FindAllMarkers(object = caudalpons_subset2,
                                        only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25)

top20_LineageID <- all_markers_LineageID %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

write.csv(all_markers_LineageID, "Plots/all_markers_LineageID.csv" )
write.csv(top20_LineageID, "Plots/top20_LineageID.csv" )
save(all_markers_LineageID, file = "Plots/all_markers_LineageID.RData")

#Figure6Q - Progenitor Markers
p  <- FeaturePlot(object= caudalpons_subset2, features = c('Id1','Nes','Hes5'), pt.size = 0.5, cols = firemap ) & PlotTheme & NoAxes() & scale_x_reverse()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]] 
p <-  ggarrange(p1, p2, p3, ncol = 3, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure6_Feature_Progenitor.png", p, base_width = 12, base_height = 4)

#Figure6R - Lineage 1 Markers
p  <- FeaturePlot(object= caudalpons_subset2, features = c('Neurod6','Nrgn','Foxb1'), pt.size = 0.5,cols = firemap ) & PlotTheme & NoAxes() & scale_x_reverse()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]] 
p <-  ggarrange(p1, p2, p3, ncol = 3, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure6_Feature_lineage1.png", p, base_width = 12, base_height = 4)

#Figure6S - Lineage 2 Markers
p  <- FeaturePlot(object= caudalpons_subset2, features = c('Lmo3','Cartpt','C1ql4'), pt.size = 0.5, cols = firemap ) & PlotTheme & NoAxes() & scale_x_reverse() 
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]] 
p <-  ggarrange(p1, p2, p3, ncol = 3, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure6_Feature_lineage2.png", p, base_width = 12, base_height = 4)

#Figure6T - Lineage 3 Markers
p  <- FeaturePlot(object= caudalpons_subset2, features = c('Hoxa3', 'Maf','Lamp5'), pt.size = 0.5, cols = firemap ) & PlotTheme & NoAxes() & scale_x_reverse()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]] 
p <-  ggarrange(p1, p2, p3, ncol = 3, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure6_Feature_lineage3.png", p, base_width = 12, base_height = 4)

Prog_genes <- c('Id1','Nes','Hes5')
Prog_genes %in% Peak_list

Lin1_genes <- c('Neurod6','Nrgn', 'Foxb1')
Lin1_genes %in% Peak_list

Lin2_genes <- c('Lmo3','Cartpt','C1ql4')
Lin2_genes %in% Peak_list

Lin3_genes <- c('Hoxa3', 'Maf','Lamp5')
Lin3_genes %in% Peak_list

