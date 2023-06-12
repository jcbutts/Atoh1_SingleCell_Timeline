### Figure 2

#Set directory as you see fit
setwd("C:/Jessie/AllSamples/Publication/Figure2")

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

##Dependencies
Initial_IDs <- "C:/Jessie/AllSamples/Publication/RenameCluster_Files/Initial_IDs.txt"
Mature_IDs_coarse <- "C:/Jessie/AllSamples/Publication/RenameCluster_Files/Mature_IDs_coarse.txt"
Mature_IDs_fine <- "C:/Jessie/AllSamples/Publication/RenameCluster_Files/Mature_IDs_fine.txt"

#Combined SCT Seurat Object -- FullDataset_normalizedRNA
load("C:/Jessie/AllSamples/Publication/Figure1/FullDataset_normalizedRNA.Rdata")

#From the combined Seurat objects folder
#load("C:/Jessie/AllSamples/Publication/Seurat_Objects/FullDataset_normalizedRNA.Rdata")

PlotTheme <- theme(plot.title=element_text(size=50), legend.key.size = unit(0.75, "cm"), legend.text = element_text(size=30)) 

#### Script ####
#Figure2A - Feature Plot Cell State
p <- FeaturePlot(object= CombinedSCT, features = c("Mki67","Sox2","Nhlh1", 'Mapt'), pt.size = 0.5, cols = firemap) & NoAxes() & 
  theme(plot.title=element_text(size=80), legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=50))
p1 <- p[[1]]
p2 <- p[[2]] 
p3 <- p[[3]]
p4 <- p[[4]]
p <-  ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure2_Feature_FullData_CellState.png", p, base_width = 32, base_height = 8)

#Figure2B - Rostral Caudal identity
#Rostral RL genes
p <- FeaturePlot(object= CombinedSCT, features = c("En1","En2","Hoxa2", 'Hoxa5'), pt.size = 0.5, cols = firemap) & NoAxes() & 
  theme(plot.title=element_text(size=80), legend.key.size = unit(1.5, "cm"), legend.text = element_text(size=50))
p1 <- p[[1]]
p2 <- p[[2]] 
p3 <- p[[3]]
p4 <- p[[4]]
p <-  ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure2_Feature_FullData_RostralCaudal.png", p, base_width = 32, base_height = 8)

#Figure2Supp
p <- FeaturePlot(object= CombinedSCT, features = c("Barhl1","Barhl2", "Lhx2", "Lhx9"), pt.size = 0.5, cols = firemap) & NoAxes() & PlotTheme
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]]
p4 <- p[[4]]
p <-  ggarrange(p1, p2, p3, p4, ncol = 4, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure2Supp_AtohTargets.png", p, base_width = 32, base_height = 8)

##DotPlot of genes to determine Initial ID
my_levels = c("38" ,"36" ,"31" ,"25", "10"  ,"22" ,"1" ,"32" ,"30" ,"5" ,"16"  ,"9" ,"26"  ,"17" ,"28" ,"13", "33", "34", "24", "21",'7','8', '35','0','3','29','37','6','15','2','11','4','19','18','39','27','12','14','23','20')
levels(CombinedSCT) <- my_levels

genes <- rev(c("Mki67", "Nes","Sox2","En1", "En2","Hoxa2", "Hoxa5",
               "Hoxa7","Pax6","Barhl1","Lhx9", "Slc17a6","Sst"))

#Figure2Supp
p <- DotPlot(CombinedSCT, features = genes, col.max = 4, col.min = 0)
save_plot("Plots/Figure2Supp_InitialID_DotPlot.pdf", p, base_width = 11, base_height = 8) 

##Rename Clusters by initial ID
Idents(CombinedSCT) <- 'seurat_clusters'

CombinedSCT <- RenameClusterGroups(CombinedSCT, Initial_IDs)

Initial.ID.levels <- c("Progenitor", "Early RLS",  "Cb RLS", "Early Pons", "Caudal Pons", "CES", "AES", "PES","CLS","Spinal Cord" ,"SST","Non RL")
Initial.ID.cols <- c(Progenitor_cols, Early_RLS_col, CBNuc_cols, Early_pons_col, Caudal_pons_col , CES_col, AES_col, PES_col, CLS_col,"navy","slategray", "gray")

CombinedSCT[["Initial.ID"]] <- Idents(object = CombinedSCT)
CombinedSCT$Initial.ID <- factor(CombinedSCT$Initial.ID, levels = Initial.ID.levels)

#Figure2C
p <- DimPlot(CombinedSCT, pt.size = 0.5, label = FALSE, group.by = 'Initial.ID' , cols = Initial.ID.cols)  + NoAxes() + PlotTheme +
  guides(colour = guide_legend(override.aes = list(size=10)))
save_plot("Plots/Figure2_UMAP_by_InitialID.png", p, base_width = 11, base_height = 8) 

#Find markers by initial ID
CombinedSCT.small <- subset(CombinedSCT, downsample = 5000)

all_markers_InitialID <- FindAllMarkers(object = CombinedSCT.small,
                                        only.pos = TRUE,
                                        min.pct = 0.25,
                                        logfc.threshold = 0.25)

top10_InitialID <- all_markers_InitialID %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(all_markers_InitialID, "Plots/all_markers_InitialID.csv" )
write.csv(top10_InitialID, "Plots/top10_InitialID.csv" )
save(all_markers_InitialID, file = "Plots/all_markers_InitialID.RData")

Initial.ID.genes <- rev(c("Id3", "Ube2c", "C1ql4", "Irx3", "Neurod1", "Reln", "Lgals1","Onecut2",
                          "C1ql3", "Pou3f1" ,"Mafb", "Atoh7", "Nfix", "Nfib", "Pid1", "Zfhx3", "Pdzrn4", "Rbp1",
                          "Lhx1", "Snhg11", "Sst", "Syt4", "Phox2b", "Ttr"))
CombinedSCT$Initial.ID <- factor(CombinedSCT$Initial.ID, levels = rev(Initial.ID.levels))

#Figure2D
p <- DotPlot(CombinedSCT, features = Initial.ID.genes, group.by = "Initial.ID", col.max = 4, col.min = 0) + 
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = 'bold'),
        legend.key.size = unit(0.5, "cm"), legend.text = element_text(size=10), legend.title = element_text(size = 10))
save_plot("Plots/Figure2_InitialID_genesofinterest_DotPlot.pdf", p, base_width = 7, base_height = 4) 

rm(CombinedSCT.small)

#### Mature Subset ####
p <- DimPlot(object = CombinedSCT, pt.size = 0.5, cells.highlight = WhichCells(object = CombinedSCT, expression = Mapt > 2.5), cols.highlight = 'steelblue4') +
  ggtitle("Mapt > 2.5") +
  theme(plot.title = element_text(size = 40, hjust = 0.5, face = "bold")) +
  NoLegend() + NoAxes()
save_plot("Plots/Figure2_Mapt_expression_over2pt5.png", p, base_width = 8, base_height = 8)

CombinedSCT[["old.ident"]] <- Idents(object = CombinedSCT)

Mature_subset <- subset(x = CombinedSCT, subset = Mapt > 2.5)

rm(CombinedSCT)

#redo colors so that they are the same from original data set
OrigIdentCols <- OrigIdentCols
SampleCols <- SampleCols[2:12]
TimepointCols <- TimepointCols[2:7]
OrigIdentCols <- OrigIdentCols[c(3,5:24)]

png("InitialQC/VlnQCbySample.png", width = 2000, height = 400)
VlnPlot(object = Mature_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.blood"), group.by = "orig.ident", cols = OrigIdentCols, ncol = 4, pt.size = 0)
dev.off()

Mature_subset <- SCTransform(Mature_subset, vars.to.regress = c("percent.mito","percent.blood"), verbose = FALSE)

Mature_subset <- RunPCA(Mature_subset, verbose = FALSE)
png("InitialQC/QC_PCHeatmap.png", width = 500, height = 1000)
DimHeatmap(Mature_subset, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png("InitialQC/QC_ElbowPlot.png", width = 500, height = 500)
ElbowPlot(Mature_subset, ndims = 50)
dev.off()

Mature_subset <- RunUMAP(Mature_subset, dims = 1:50, verbose = FALSE)
Mature_subset <- FindNeighbors(Mature_subset, dims = 1:50, verbose = FALSE)

#Res at 50PC
find_res(Mature_subset)
MinDist(Mature_subset, 50, 0.8, 1)
Spread(Mature_subset, 50, 0.8, 0.2)

Mature_subset <- RunUMAP(Mature_subset, spread =3, min.dist = 0.2, dims = 1:50, verbose = FALSE)
Mature_subset <- FindNeighbors(Mature_subset, dims = 1:50, verbose = FALSE)
Mature_subset <- FindClusters(Mature_subset, verbose = FALSE, resolution= 0.8)

MetaPlots_UMAP(Mature_subset, 0.5)

p <- DimPlot(Mature_subset, group.by = 'Sample', cols = SampleCols, pt.size = 0.5, label = FALSE) + NoAxes() +
  guides(colour = guide_legend(override.aes = list(size=14))) +
  theme(legend.text = element_text(size = 40, face = 'bold'))
save_plot("Plots/Figure2_UMAP_MatureSubset_bySample.png", p, base_width = 13, base_height = 8)

p <- DimPlot(Mature_subset, pt.size = 0.5, label = TRUE, label.size = 12) + theme(legend.position = "none") + NoAxes()
save_plot("Plots/Fiugre2Supp_UMAP_MatureSubset_byCluster.png", p, base_width = 8, base_height = 8)

save(Mature_subset, file = "Mature_subset.Rdata")

#RNA Expression
DefaultAssay(Mature_subset) <- "RNA"
Mature_subset <- NormalizeData(Mature_subset)

save(Mature_subset, file = "Mature_subset_normalizedRNA.Rdata")

all_markers_Mature <- FindAllMarkers(object = Mature_subset,
                                               only.pos = TRUE,
                                               min.pct = 0.25,
                                               logfc.threshold = 0.25)

top10_Mature_subset <- all_markers_Mature %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(all_markers_Mature, "Plots/all_markers_Mature_subset.csv" )
write.csv(top10_Mature_subset, "Plots/top10_Mature_subset.csv" )
save(all_markers_Mature, file = "Plots/all_markers_Mature_subset8.RData")

n_cells <- FetchData(Mature_subset, vars = c("ident", "Sample")) %>%
  dplyr::count(ident, Sample) %>%
  tidyr::spread(ident, n)
write.csv(n_cells, "Plots/nCells_byCluster.csv" )

p <- DotPlot(Mature_subset,features =  c("Hoxc8","Hoxb7","Hoxc6","Hoxa5","Hoxa4","Hoxb4","Hoxa3","Hoxb3","Hoxb2","Hoxa2","En1","En2")) + 
  theme(axis.text=element_text(size=24), axis.title = element_text(size=30),legend.text=element_text(size=16)) 
save_plot("Plots/Supp_Dotplot_Hox.png", p, base_width = 16, base_height = 12)

p <- DotPlot(Mature_subset,features =  c("Lhx2","Lhx9","Barhl2","Barhl1")) +
  theme(axis.text=element_text(size=24), axis.title = element_text(size=30),legend.text=element_text(size=16)) 
save_plot("Plots/Supp_Dotplot_Atoh1Targets.png", p, base_width = 8, base_height = 12)

p <- DotPlot(Mature_subset,features =  c("Slc17a6","Slc17a7","Sst","Th","Tacr1","Slc18a3")) +
  theme(axis.text=element_text(size=24), axis.title = element_text(size=30),legend.text=element_text(size=16)) 
save_plot("Plots/Supp_Dotplot_Neurotransmitters.png", p, base_width = 10, base_height = 12)


#Rename Clusters by coarse ID
Idents(Mature_subset) <- 'seurat_clusters'
Mature_subset <- RenameClusterGroups(Mature_subset, Mature_IDs_coarse)

p <- DimPlot(Mature_subset, label = TRUE, pt.size = 0.5, label.size = 6, repel = TRUE)  + NoLegend() + NoAxes()
save_plot("Plots/Figure2_UMAP_MatureLabels_coarse.png", p, base_width = 10, base_height = 8) 

all_markers_Mature_coarse <- FindAllMarkers(object = Mature_subset,
                                     only.pos = TRUE,
                                     min.pct = 0.25,
                                     logfc.threshold = 0.25)

top10_Mature_coarse <- all_markers_Mature_coarse %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(all_markers_Mature_coarse, "Plots/all_markers_Mature_coarse.csv" )
write.csv(top10_Mature_coarse, "Plots/top10_Mature_coarse.csv" )
save(all_markers_Mature_coarse, file = "Plots/all_markers_Mature_coarse.csv")


##Rename Clusters by fine ID
Idents(Mature_subset) <- 'seurat_clusters'
Mature_subset <- RenameClusterGroups(Mature_subset, Mature_IDs_fine)

p <- DimPlot(Mature_subset, label = TRUE, pt.size = 0.5, label.size = 6, repel = TRUE)  + NoLegend() + NoAxes()
save_plot("Plots/Figure2Supp_UMAP_MatureLabels_fine.png", p, base_width = 10, base_height = 8) 

p <- FeaturePlot(object= Mature_subset, features = c("Tcf24","Kcnip2", "Sox7", "Arid5a", "Spx") , pt.size = 0.5, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]] 
p2 <- p[[2]] 
p3 <- p[[3]] 
p4 <- p[[4]] 
p5 <- p[[5]] 
p <-  ggarrange(p1, p2, p3, p4, p5, ncol = 5, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure2_NewMarkers_Mature_Feature.png", p, base_width = 40, base_height = 8)

p <- FeaturePlot(object= Mature_subset, features = "Rora", pt.size = 0.5, cols = firemap) + PlotTheme + NoAxes()
save_plot("Plots/Figure2Supp_Feature_Mature_Rora.png", p, base_width = 9, base_height = 8)

p <- FeaturePlot(object= Mature_subset, features = "Rgs4", pt.size = 0.5, cols = firemap) + PlotTheme + NoAxes()
save_plot("Plots/Figure2Supp_Feature_Mature_Rgs4.png", p, base_width = 9, base_height = 8)

p <- FeaturePlot(object= Mature_subset, features = "Pou4f1", pt.size = 0.5, cols = firemap) + PlotTheme + NoAxes()
save_plot("Plots/Figure2Supp_Feature_Mature_Pou4f1.png", p, base_width = 9, base_height = 8)








