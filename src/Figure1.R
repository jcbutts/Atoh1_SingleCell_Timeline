### Figure 1

#Set directory as you see fit
setwd("C:/Jessie/AllSamples/Publication/Figure1")

options(future.globals.maxSize= 500000000000)
memory.limit(size = 175000)

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

###Libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggpubr)
library(data.table)

### Dependencies
#Combined SCT Seurat Object before subset -- Initial_FullDataset_normalizedRNA
load("C:/Jessie/AllSamples/Publication/Upstream_analysis/Initial_FullDataset_normalizedRNA.Rdata")

#From the combined Seurat objects folder
#load("C:/Jessie/AllSamples/Publication/Seurat_Objects/Initial_FullDataset_normalizedRNA.Rdata")

PlotTheme <- theme(plot.title=element_text(size=50), legend.key.size = unit(0.75, "cm"), legend.text = element_text(size=30)) 

#### Script #####

### Subset out spurious clusters
CombinedSCT <- subset(x = CombinedSCT, idents = c("19","43","32","34","42","44","41"), invert = TRUE)

png("InitialQC/VlnQCbySample.png", width = 2000, height = 400)
VlnPlot(object = CombinedSCT, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.blood"), group.by = "orig.ident", cols = OrigIdentCols, ncol = 4, pt.size = 0)
dev.off()

png("InitialQC/DotplotCbySample_Prefilter.png", width = 1000, height = 400)
plot1 <- FeatureScatter(object = CombinedSCT, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 0.5)
plot2 <- FeatureScatter(object = CombinedSCT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.5)
ggarrange(plot1, plot2, legend = 'right', common.legend = TRUE)
dev.off()

CombinedSCT <- SCTransform(CombinedSCT, vars.to.regress = c("percent.mito","percent.blood"), verbose = FALSE)

##PC QC
CombinedSCT <- RunPCA(CombinedSCT, verbose = FALSE)
png("InitialQC/QC_PCHeatmap.png", width = 500, height = 1000)
DimHeatmap(CombinedSCT, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png("InitialQC/QC_ElbowPlot.png", width = 500, height = 500)
ElbowPlot(CombinedSCT, ndims = 50)
dev.off()

#Test Different PC and res
CombinedSCT <- RunUMAP(CombinedSCT, dims = 1:50, verbose = FALSE)
CombinedSCT <- FindNeighbors(CombinedSCT, dims = 1:50, verbose = FALSE)

find_res(CombinedSCT)

##Spread and MinDist at 50PC
MinDist(CombinedSCT, 50, 0.8,1)
Spread(CombinedSCT, 50, 0.8,0.2)


######Final Setup ######
CombinedSCT <- RunUMAP(CombinedSCT, spread = 3, min.dist = 0.2,  dims = 1:50, verbose = FALSE)
CombinedSCT <- FindNeighbors(CombinedSCT, dims = 1:50, verbose = FALSE)
CombinedSCT <- FindClusters(CombinedSCT, verbose = FALSE, resolution= 0.8)

### Characterization of metadata
MetaPlots_UMAP(CombinedSCT, 0.5)

p <- DimPlot(CombinedSCT, pt.size = 0.5, label = TRUE, label.size = 12) + theme(legend.position = "none") + NoAxes()
save_plot("Plots/Figure1_UMAP_FullDataset_clusters.png", p, base_width = 8, base_height = 8)

p <- DimPlot(CombinedSCT, group.by = 'Sample', cols = SampleCols, pt.size = 0.5, label = FALSE) + NoAxes() +
  guides(colour = guide_legend(override.aes = list(size=14))) +
  theme(legend.text = element_text(size = 40, face = 'bold'))
save_plot("Plots/Figure1_UMAP_FullDataset_bySample.png", p, base_width = 13, base_height = 8)

save(CombinedSCT, file = "FullDataset.Rdata")

##### RNA Expression
DefaultAssay(CombinedSCT) <- "RNA"
CombinedSCT <- NormalizeData(CombinedSCT)

save(CombinedSCT, file = "FullDataset_normalizedRNA.Rdata")

all_markers_FUllDataset <- FindAllMarkers(object = CombinedSCT,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

top10_FUllDataset <- all_markers_FUllDataset %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(all_markers_FUllDataset, "Plots/all_markers_FullDataset.csv" )
write.csv(top10_FUllDataset, "Plots/top10_FullDataset.csv" )
save(all_markers_FUllDataset, file = "Plots/all_markers_FullDataset.RData")

###Figure 1F 
p <- FeaturePlot(object= CombinedSCT, features = c('Atoh1','WPRE'), pt.size = 1, cols = firemap) & PlotTheme & NoAxes()
p1 <- p[[1]] 
p2 <- p[[2]]
p <-  ggarrange(p1, p2, common.legend = TRUE, legend = 'right')
save_plot("Plots/Figure1_Feature_Atoh1_WPRE.png", p, base_width = 16, base_height = 8)



