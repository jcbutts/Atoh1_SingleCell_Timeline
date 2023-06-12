setwd("C:/Jessie/AllSamples/Publication/Upstream_analysis")

#Uses the regressblood and filter blood code at the beginning

#Scale up to 50GB
options(future.globals.maxSize= 500000000000)
memory.limit(size = 175000)

library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(magrittr)

###Directories
dir.create("Plots")
dir.create("SpuriousClusters")

##Load in samples that have already had doublets removed
##Samples already filtered by: E12.5_Td <- subset(x = E12.5_Td, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & nCount_RNA  > 500 & percent.mito < 0.4)

###### Load Metadata #####
#E9.5
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E9.5_rep1_rmdp.Rdata")

load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E9.5_rep2_rmdp.Rdata")

#10.5
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E10.5_rep1_rmdp.Rdata")

load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E10.5_rep2_rmdp.Rdata")

#11.5
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E11.5_TdOnly_rep1_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E11.5_GFPTd_rep1_rmdp.Rdata")

load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E11.5_TdOnly_rep2_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E11.5_GFPTd_rep2_rmdp.Rdata")

#12.5
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E12.5_TdOnly_rep1_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E12.5_GFPTd_rep1_rmdp.Rdata")

load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E12.5_TdOnly_rep2_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E12.5_GFPTd_rep2_rmdp.Rdata")

#13.5
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E13.5_TdOnly_rep1_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E13.5_GFPTd_rep1_rmdp.Rdata")

load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E13.5_TdOnly_rep2_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E13.5_GFPTd_rep2_rmdp.Rdata")

#14.5
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E14.5_TdOnly_rep1_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E14.5_GFPTd_rep1_rmdp.Rdata")

load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E14.5_TdOnly_rep2_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E14.5_GFPTd_rep2_rmdp.Rdata")

#16.5
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E16.5_TdOnly_rep1_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_RemapSamples/E16.5_GFPTd_rep1_rmdp.Rdata")

load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E16.5_TdOnly_rep2_rmdp.Rdata")
load("C:/Jessie/CombinedAnalysis/DoubletAnalysis/DoubletRemoval_03_Deep_GARP/E16.5_GFPTd_rep2_rmdp.Rdata")



#Set some metadata to pull out later 
#9.5
E9.5_GFP@meta.data[,"color"] <- "GFP"
E9.5_GFP@meta.data[,"TimePoint"] <- "E9.5"
E9.5_GFP@meta.data[,"Replicate"] <- "1"
E9.5_GFP@meta.data[,"Sample"] <- "E9.5 GFP"

E9.5_GFP_rep2@meta.data[,"color"] <- "GFP"
E9.5_GFP_rep2@meta.data[,"TimePoint"] <- "E9.5"
E9.5_GFP_rep2@meta.data[,"Replicate"] <- "2"
E9.5_GFP_rep2@meta.data[,"Sample"] <- "E9.5 GFP"

#10.5
E10.5_GFP@meta.data[,"color"] <- "GFP"
E10.5_GFP@meta.data[,"TimePoint"] <- "E10.5"
E10.5_GFP@meta.data[,"Replicate"] <- "1"
E10.5_GFP@meta.data[,"Sample"] <- "E10.5 GFP"

E10.5_GFP_rep2@meta.data[,"color"] <- "GFP"
E10.5_GFP_rep2@meta.data[,"TimePoint"] <- "E10.5"
E10.5_GFP_rep2@meta.data[,"Replicate"] <- "2"
E10.5_GFP_rep2@meta.data[,"Sample"] <- "E10.5 GFP"

#11.5
E11.5_Td@meta.data[,"color"] <- "Td"
E11.5_Td@meta.data[,"TimePoint"] <- "E11.5"
E11.5_Td@meta.data[,"Replicate"] <- "1"
E11.5_Td@meta.data[,"Sample"] <- "E11.5 Td"

E11.5_GFPTd@meta.data[,"color"] <- "GFPTd"
E11.5_GFPTd@meta.data[,"TimePoint"] <- "E11.5"
E11.5_GFPTd@meta.data[,"Replicate"] <- "1"
E11.5_GFPTd@meta.data[,"Sample"] <- "E11.5 GFPTd"

E11.5_Td_rep2@meta.data[,"color"] <- "Td"
E11.5_Td_rep2@meta.data[,"TimePoint"] <- "E11.5"
E11.5_Td_rep2@meta.data[,"Replicate"] <- "2"
E11.5_Td_rep2@meta.data[,"Sample"] <- "E11.5 Td"

E11.5_GFPTd_rep2@meta.data[,"color"] <- "GFPTd"
E11.5_GFPTd_rep2@meta.data[,"TimePoint"] <- "E11.5"
E11.5_GFPTd_rep2@meta.data[,"Replicate"] <- "2"
E11.5_GFPTd_rep2@meta.data[,"Sample"] <- "E11.5 GFPTd"

#12.5
E12.5_Td@meta.data[,"color"] <- "Td"
E12.5_Td@meta.data[,"TimePoint"] <- "E12.5"
E12.5_Td@meta.data[,"Replicate"] <- "1"
E12.5_Td@meta.data[,"Sample"] <- "E12.5 Td"

E12.5_GFPTd@meta.data[,"color"] <- "GFPTd"
E12.5_GFPTd@meta.data[,"TimePoint"] <- "E12.5"
E12.5_GFPTd@meta.data[,"Replicate"] <- "1"
E12.5_GFPTd@meta.data[,"Sample"] <- "E12.5 GFPTd"

E12.5_Td_rep2@meta.data[,"color"] <- "Td"
E12.5_Td_rep2@meta.data[,"TimePoint"] <- "E12.5"
E12.5_Td_rep2@meta.data[,"Replicate"] <- "2"
E12.5_Td_rep2@meta.data[,"Sample"] <- "E12.5 Td"

E12.5_GFPTd_rep2@meta.data[,"color"] <- "GFPTd"
E12.5_GFPTd_rep2@meta.data[,"TimePoint"] <- "E12.5"
E12.5_GFPTd_rep2@meta.data[,"Replicate"] <- "2"
E12.5_GFPTd_rep2@meta.data[,"Sample"] <- "E12.5 GFPTd"

#13.5
E13.5_Td@meta.data[,"color"] <- "Td"
E13.5_Td@meta.data[,"TimePoint"] <- "E13.5"
E13.5_Td@meta.data[,"Replicate"] <- "1"
E13.5_Td@meta.data[,"Sample"] <- "E13.5 Td"

E13.5_GFPTd@meta.data[,"color"] <- "GFPTd"
E13.5_GFPTd@meta.data[,"TimePoint"] <- "E13.5"
E13.5_GFPTd@meta.data[,"Replicate"] <- "1"
E13.5_GFPTd@meta.data[,"Sample"] <- "E13.5 GFPTd"

E13.5_Td_rep2@meta.data[,"color"] <- "Td"
E13.5_Td_rep2@meta.data[,"TimePoint"] <- "E13.5"
E13.5_Td_rep2@meta.data[,"Replicate"] <- "2"
E13.5_Td_rep2@meta.data[,"Sample"] <- "E13.5 Td"

E13.5_GFPTd_rep2@meta.data[,"color"] <- "GFPTd"
E13.5_GFPTd_rep2@meta.data[,"TimePoint"] <- "E13.5"
E13.5_GFPTd_rep2@meta.data[,"Replicate"] <- "2"
E13.5_GFPTd_rep2@meta.data[,"Sample"] <- "E13.5 GFPTd"

#E14.5
E14.5_Td@meta.data[,"color"] <- "Td"
E14.5_Td@meta.data[,"TimePoint"] <- "E14.5"
E14.5_Td@meta.data[,"Replicate"] <- "1"
E14.5_Td@meta.data[,"Sample"] <- "E14.5 Td"

E14.5_GFPTd@meta.data[,"color"] <- "GFPTd"
E14.5_GFPTd@meta.data[,"TimePoint"] <- "E14.5"
E14.5_GFPTd@meta.data[,"Replicate"] <- "1"
E14.5_GFPTd@meta.data[,"Sample"] <- "E14.5 GFPTd"

E14.5_Td_rep2@meta.data[,"color"] <- "Td"
E14.5_Td_rep2@meta.data[,"TimePoint"] <- "E14.5"
E14.5_Td_rep2@meta.data[,"Replicate"] <- "2"
E14.5_Td_rep2@meta.data[,"Sample"] <- "E14.5 Td"

E14.5_GFPTd_rep2@meta.data[,"color"] <- "GFPTd"
E14.5_GFPTd_rep2@meta.data[,"TimePoint"] <- "E14.5"
E14.5_GFPTd_rep2@meta.data[,"Replicate"] <- "2"
E14.5_GFPTd_rep2@meta.data[,"Sample"] <- "E14.5 GFPTd"

#E16.5
E16.5_Td@meta.data[,"color"] <- "Td"
E16.5_Td@meta.data[,"TimePoint"] <- "E16.5"
E16.5_Td@meta.data[,"Replicate"] <- "1"
E16.5_Td@meta.data[,"Sample"] <- "E16.5 Td"

E16.5_GFPTd@meta.data[,"color"] <- "GFPTd"
E16.5_GFPTd@meta.data[,"TimePoint"] <- "E16.5"
E16.5_GFPTd@meta.data[,"Replicate"] <- "1"
E16.5_GFPTd@meta.data[,"Sample"] <- "E16.5 GFPTd"

E16.5_Td_rep2@meta.data[,"color"] <- "Td"
E16.5_Td_rep2@meta.data[,"TimePoint"] <- "E16.5"
E16.5_Td_rep2@meta.data[,"Replicate"] <- "2"
E16.5_Td_rep2@meta.data[,"Sample"] <- "E16.5 Td"

E16.5_GFPTd_rep2@meta.data[,"color"] <- "GFPTd"

E16.5_GFPTd_rep2@meta.data[,"TimePoint"] <- "E16.5"
E16.5_GFPTd_rep2@meta.data[,"Replicate"] <- "2"
E16.5_GFPTd_rep2@meta.data[,"Sample"] <- "E16.5 GFPTd"

#Merge Objects

CombinedSCT <- merge(x = E9.5_GFP , y = c(E9.5_GFP_rep2, E10.5_GFP, E10.5_GFP_rep2,
                                          E11.5_Td, E11.5_GFPTd, E11.5_Td_rep2, E11.5_GFPTd_rep2,
                                          E12.5_Td, E12.5_GFPTd, E12.5_Td_rep2, E12.5_GFPTd_rep2,
                                          E13.5_Td, E13.5_GFPTd, E13.5_Td_rep2, E13.5_GFPTd_rep2,
                                          E14.5_Td, E14.5_GFPTd, E14.5_Td_rep2, E14.5_GFPTd_rep2,
                                          E16.5_Td, E16.5_GFPTd, E16.5_Td_rep2, E16.5_GFPTd_rep2), 
                     add.cell.ids=c("E9.5_GFP.data", "E9.5_GFP_rep2.data", "E10.5_GFP.data", "E10.5_GFP_rep2.data",
                                    "E11.5_Td.data", "E11.5_GFPTd.data", "E11.5_Td_rep2.data", "E11.5_GFPTd_rep2.data",  
                                    "E12.5_Td.data", "E12.5_GFPTd.data", "E12.5_Td_rep2.data", "E12.5_GFPTd_rep2.data",
                                    "E13.5_Td.data", "E13.5_GFPTd.data", "E13.5_Td_rep2.data", "E13.5_GFPTd_rep2.data",
                                    "E14.5_Td.data", "E14.5_GFPTd.data", "E14.5_Td_rep2.data", "E14.5_GFPTd_rep2.data",
                                    "E16.5_Td.data", "E16.5_GFPTd.data", "E16.5_Td_rep2.data", "E16.5_GFPTd_rep2.data"))



rm(E9.5_GFP, E9.5_GFP_rep2, E10.5_GFP, E10.5_GFP_rep2, E11.5_GFPTd, E11.5_Td, E11.5_GFPTd_rep2, E11.5_Td_rep2,
   E12.5_GFPTd, E12.5_Td, E12.5_GFPTd_rep2, E12.5_Td_rep2, E13.5_GFPTd, E13.5_Td, E13.5_GFPTd_rep2, E13.5_Td_rep2, 
   E14.5_GFPTd, E14.5_Td, E14.5_GFPTd_rep2, E14.5_Td_rep2, E16.5_GFPTd, E16.5_Td, E16.5_GFPTd_rep2, E16.5_Td_rep2)

#Remove Blood Genes
#counts <- GetAssayData(CombinedSCT, assay = "RNA")
#counts <- counts[-(which(rownames(counts) %in% c('Hbb-bs','Hbb-bt','Hba-a1','Hbb-y','Hbb-a2'))),]
#CombinedSCT <- subset(CombinedSCT, features = rownames(counts))


#Reorder objects for graphing 
CombinedSCT$orig.ident <- factor(CombinedSCT$orig.ident, levels = c("E9.5_GFP", "E9.5_GFP_rep2", "E10.5_GFP", "E10.5_GFP_rep2",
                                                                    "E11.5_GFPTd", "E11.5_GFPTd_rep2", "E11.5_Td", "E11.5_Td_rep2",
                                                                    "E12.5_GFPTd", "E12.5_GFPTd_rep2", "E12.5_Td", "E12.5_Td_rep2",
                                                                    "E13.5_GFPTd", "E13.5_GFPTd_rep2", "E13.5_Td", "E13.5_Td_rep2",  
                                                                    "E14.5_GFPTd", "E14.5_GFPTd_rep2", "E14.5_Td", "E14.5_Td_rep2",
                                                                    "E16.5_GFPTd", "E16.5_GFPTd_rep2", "E16.5_Td", "E16.5_Td_rep2"))
CombinedSCT$Sample <- factor(CombinedSCT$Sample, levels = c("E9.5 GFP", "E10.5 GFP", "E11.5 GFPTd", "E11.5 Td", "E12.5 GFPTd", "E12.5 Td","E13.5 GFPTd", "E13.5 Td","E14.5 GFPTd", "E14.5 Td", "E16.5 GFPTd", "E16.5 Td" ))
CombinedSCT$TimePoint <- factor(CombinedSCT$TimePoint, levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5","E14.5","E16.5") )

#color options
library(RColorBrewer)
library(viridis)
num_timepoints <- 7 
num_origident <- 24 
num_samples <- 12

TimepointCols <- brewer.pal(num_timepoints+1,'PuBu')
TimepointCols <- TimepointCols[2:length(TimepointCols)]
SampleCols <- brewer.pal(num_samples,"Paired")
OrigIdentCols <- viridis(num_origident, option="viridis")
OrigIdentCols2 <- viridis(num_origident, option="inferno")
ColorCols <- c("springgreen4","gold2","tomato3")
firemap <- c("gray80","orange","firebrick2","firebrick3","firebrick4")

##QC 
mito.features_CombinedSCT <- grep(pattern = "^mt.", x = rownames(x = CombinedSCT), value = TRUE)
percent.mito_CombinedSCT <- Matrix::colSums(x = GetAssayData(object = CombinedSCT, slot = 'counts')[mito.features_CombinedSCT, ]) / Matrix::colSums(x = GetAssayData(object = CombinedSCT, slot = 'counts'))

CombinedSCT[['percent.mito']] <- percent.mito_CombinedSCT

#Pull out blood genes for regression
hbb.features_CombinedSCT <- grep(pattern = "^Hbb.", x = rownames(x = CombinedSCT), value = TRUE)
hba.features_CombinedSCT <- grep(pattern = "^Hba.", x = rownames(x = CombinedSCT), value = TRUE)
blood.features_CombinedSCT <- c(hbb.features_CombinedSCT, hba.features_CombinedSCT)
percent.blood_CombinedSCT <- Matrix::colSums(x = GetAssayData(object = CombinedSCT, slot = 'counts')[blood.features_CombinedSCT, ])/ Matrix::colSums(x = GetAssayData(object = CombinedSCT, slot = 'counts'))

CombinedSCT[['percent.blood']] <- percent.blood_CombinedSCT

png("Plots/VlnQCbySample_Prefilter.png", width = 2000, height = 400)
VlnPlot(object = CombinedSCT, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.blood"), group.by = "orig.ident", cols = OrigIdentCols, ncol = 4, pt.size = 0)
dev.off()

png("Plots/PercentBlood_PreFilter_ByOrigIdent.png", width = 800, height = 400)
VlnPlot(object = CombinedSCT, feature = "percent.blood", group.by = "orig.ident", pt.size = 0.5, cols = OrigIdentCols)
dev.off()

png("Plots/PercentBlood_PreFilter_Zoom_ByOrigIdent.png", width = 800, height = 400)
VlnPlot(object = CombinedSCT, feature = "percent.blood", group.by = "orig.ident", pt.size = 0.5, cols = OrigIdentCols, y.max = 0.15)
dev.off()

png("Plots/DotplotCbySample_Prefilter.png", width = 1000, height = 400)
plot1 <- FeatureScatter(object = CombinedSCT, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 0.5)
plot2 <- FeatureScatter(object = CombinedSCT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.5)
ggarrange(plot1, plot2, legend = 'right', common.legend = TRUE)
dev.off()

#Filtering percent mito a little more, was at 0.4 with doublet removal
CombinedSCT <- subset(x = CombinedSCT, subset =  percent.mito < 0.1 & percent.blood < 0.1)

png("Plots/VlnQCbySample_PostFilter.png", width =2000, height = 400)
VlnPlot(object = CombinedSCT, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.blood"), group.by = "orig.ident", cols = OrigIdentCols, ncol = 4, pt.size = 0)
dev.off()

png("Plots/DotplotCbySample_PostFilter.png", width = 1000, height = 400)
plot1 <- FeatureScatter(object = CombinedSCT, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size = 0.5,)
plot2 <- FeatureScatter(object = CombinedSCT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.5,)
ggarrange(plot1, plot2, legend = 'right', common.legend = TRUE)
dev.off()

gc()

##SCT Transform on the merged object 
CombinedSCT <- SCTransform(CombinedSCT, vars.to.regress = c("percent.mito","percent.blood"), verbose = FALSE)

##PC QC
CombinedSCT <- RunPCA(CombinedSCT, verbose = FALSE)
png("Plots/QC_PCHeatmap.png", width = 500, height = 1000)
DimHeatmap(CombinedSCT, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png("Plots/QC_ElbowPlot.png", width = 500, height = 500)
ElbowPlot(CombinedSCT, ndims = 50)
dev.off()

#Final Setup
CombinedSCT <- RunUMAP(CombinedSCT, dims = 1:50, verbose = FALSE)
CombinedSCT <- FindNeighbors(CombinedSCT, dims = 1:50, verbose = FALSE)
CombinedSCT <- FindClusters(CombinedSCT, verbose = FALSE, resolution= 0.8)

png("Plots/SCT_UMAP_res08_50PC_withlegend.png", width = 1200, height = 1000)
DimPlot(CombinedSCT, repel = TRUE, pt.size = 0.5, label = TRUE, label.size = 16) 
dev.off()

#save(new_CombinedSCT, file = "Initial_FullDataset_newSCT.Rdata")
save(CombinedSCT, file = "Initial_FullDataset.Rdata")

n_cells <- FetchData(CombinedSCT, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)
write.csv(n_cells, file= "Plots/cell number in each cluster.csv")

##### GENE EXPRESSION #####
DefaultAssay(CombinedSCT) <- "RNA"
CombinedSCT <- NormalizeData(CombinedSCT)

save(CombinedSCT, file = "Initial_FullDataset_normalizedRNA.Rdata")

all_markers_Initial_fulldataset <- FindAllMarkers(object = CombinedSCT,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

top10_Initial_fulldataset <- all_markers_Initial_fulldataset %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(all_markers_Initial_fulldataset, "Plots/all_markers_Initial_fulldataset.csv" )
write.csv(top10_Initial_fulldataset, "Plots/top10_Initial_fulldataset.csv" )
save(all_markers_Initial_fulldataset, file = "Plots/all_markers__Initial_fulldataset.RData")

###Spurious Clusters

#Cluster 19 - Mature inhibitory Cells
png("SpuriousClusters/Cluster19.png", width = 2400, height = 800)
VlnPlot(object= CombinedSCT, features = c("Tfap2b","Pax2","Gad2", "Gad1","Slc6a5", "Slc32a1"), pt.size = 0, ncol = 3 ) 
dev.off()

png("SpuriousClusters/Cluster19_Localization.png", width = 400, height = 400)
DimPlot(object = CombinedSCT, reduction = 'umap', pt.size = 0.75, cells.highlight = WhichCells(object = CombinedSCT, idents = "19")) + NoLegend()
dev.off()

#Need to find better genes for the last 3
png("SpuriousClusters/Cluster32.png", width = 2400, height = 400)
VlnPlot(object= CombinedSCT, features = c("Ascl1","Ednrb","Ptprz1"), pt.size = 0, ncol = 3 ) 
dev.off()

png("SpuriousClusters/Cluster32_Localization.png", width = 400, height = 400)
DimPlot(object = CombinedSCT, reduction = 'umap', pt.size = 0.75, cells.highlight = WhichCells(object = CombinedSCT, idents = "32")) + NoLegend()
dev.off()

#Cluster 34 - Meninges 
png("SpuriousClusters/Cluster34.png", width = 2400, height = 800)
VlnPlot(object= CombinedSCT, features = c("Col1a2","Col1a1","Dcn", "Dkk2","Foxc2", "Ogn"), pt.size = 0, ncol = 3 ) 
dev.off()

png("SpuriousClusters/Cluster34_Localization.png", width = 400, height = 400)
DimPlot(object = CombinedSCT, reduction = 'umap', pt.size = 0.75, cells.highlight = WhichCells(object = CombinedSCT, idents = "34")) + NoLegend()
dev.off()

#Cluster 41 - Endothelial Cells
png("SpuriousClusters/Cluster41.png", width = 2400, height = 800)
VlnPlot(object= CombinedSCT, features = c("Klf2","Igfbp7","Ramp2", "Pecam1","Foxl2", "Foxq1"), pt.size = 0, ncol = 3 ) 
dev.off()

png("SpuriousClusters/Cluster41_Localization.png", width = 400, height = 400)
DimPlot(object = CombinedSCT, reduction = 'umap', pt.size = 0.75, cells.highlight = WhichCells(object = CombinedSCT, idents = "41")) + NoLegend()
dev.off()

#Cluster 42 - white blood cells
png("SpuriousClusters/Cluster42.png", width = 2400, height = 800)
VlnPlot(object= CombinedSCT, features = c("C1qb","Ccl4","Tyrobp", "Ccr5","Ptprc", "Irf5"), pt.size = 0, ncol = 3 ) 
dev.off()

png("SpuriousClusters/Cluster42_Localization.png", width = 400, height = 400)
DimPlot(object = CombinedSCT, reduction = 'umap', pt.size = 0.75, cells.highlight = WhichCells(object = CombinedSCT, idents = "42")) + NoLegend()
dev.off()

#Cluster 43 - Smooth muscle cells
png("SpuriousClusters/Cluster43.png", width = 2400, height = 800)
VlnPlot(object= CombinedSCT, features = c("Rgs5","Vtn","Acta2", "Ednra","Tagln", "Tbx2"), pt.size = 0, ncol = 3 )
dev.off()

png("SpuriousClusters/Cluster43_Localization.png", width = 400, height = 400)
DimPlot(object = CombinedSCT, reduction = 'umap', pt.size = 0.75, cells.highlight = WhichCells(object = CombinedSCT, idents = "43")) + NoLegend()
dev.off()

#Cluster 44 - platelets 
png("SpuriousClusters/Cluster44.png", width = 2400, height = 800)
VlnPlot(object= CombinedSCT, features = c("Ppbp","Clec1b","Fermt3", "F10","Vwf", "Gp6"), pt.size = 0, ncol = 3 )
dev.off()

png("SpuriousClusters/Cluster44_Localization.png", width = 400, height = 400)
DimPlot(object = CombinedSCT, reduction = 'umap', pt.size = 0.75, cells.highlight = WhichCells(object = CombinedSCT, idents = "44")) + NoLegend()
dev.off()

