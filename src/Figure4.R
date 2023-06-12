### Figure 4

#Set directory as you see fit
setwd("C:/Jessie/AllSamples/Publication/Figure4")

##Create Folders
dir.create("Plots")

###Functions
source('C:/Jessie/Functions/RenameClusterGroups.R')
source("C:/Jessie/Functions/StackedVlnPlot.R")
source("C:/Jessie/Functions/CalcModPVal.R")
source("C:/Jessie/Functions/ModScoreRidgePlot.R")
source("C:/Jessie/Functions/PlotEnrichment.R")

## Libraries
library(Seurat)
library(lme4)
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggpubr)
library(data.table)
library(ggplot2)
library(ggridges)

###Dependencies
#Full Dataset Seurat Object
load("C:/Jessie/AllSamples/Publication/Figure1/FullDataset_normalizedRNA.Rdata")

#From the combined Seurat objects folder
#load("C:/Jessie/AllSamples/Publication/Seurat_Objects/FullDataset_normalizedRNA.Rdata")

#Category ID Labels
Category_IDs <- "C:/Jessie/AllSamples/Publication/RenameCluster_Files/Category_IDs.txt"

#CUT&RUN peak files
E14_closest <-read.table('C:/Jessie/AllSamples/Publication/CUT&RUN_Data/E14_Subset.score.narrowPeak_Closest', sep = '\t', 
                         col.names = c("chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue","peak","local_chrom","local_chromStart","local_chromEnd","local_GeneName","GeneID","plus","Distance"))
E12_closest <-read.table('C:/Jessie/AllSamples/Publication/CUT&RUN_Data/E12_Subset.score.narrowPeak_Closest', sep = '\t', 
                         col.names = c("chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue","peak","local_chrom","local_chromStart","local_chromEnd","local_GeneName","GeneID","plus","Distance"))

#Fiji Quantification
Fiji_Quantification <- read.csv("C:/Jessie/AllSamples/Publication/Fiji_Quantification/Fiji_Quantification.csv")

#####Script#####
###Rename clusters and add to metadata by category and plot vln plots
Idents(CombinedSCT) <- 'seurat_clusters'

#Rename clusters by category
CombinedSCT <- RenameClusterGroups(CombinedSCT, Category_IDs)

##Reorder levels
Category.levels = c("Progenitor", "Intermediate Progenitor", "Migrating_1", "Migrating_2", "Mature","Non RL")
#Assign Colors
Category.cols =c(Progenitor_cols, IntermediateProgenitor_cols, Migrating_1_cols, Migrating_2_cols, Mature_cols, NonRL_cols)

#Add labels to metadata and reorder
CombinedSCT[["FullData.category"]] <- Idents(object = CombinedSCT)
CombinedSCT$FullData.category <- factor(CombinedSCT$FullData.category, levels = Category.levels)

#Figure4A - UMAP labeled by category
p <- DimPlot(CombinedSCT, pt.size = 0.5, label = FALSE, group.by = 'FullData.category', cols = Category.cols)  + NoAxes() +
  guides(colour = guide_legend(override.aes = list(size=14))) +
  theme(legend.text = element_text(size = 50, face = 'bold')) 
save_plot("Plots/Figure4_UMAP_by_category.png", p, base_width = 16, base_height = 8) 

##Vln Plot
CategoryGenes = c("Mki67", "Ube2c","Sox2","Nes","Atoh1", "Nhlh1", "Slc17a6", "Mapt")

#Figure4B - Stacked VlnPlot by category
p <- StackedVlnPlot(CombinedSCT, features = CategoryGenes, cols = Category.cols, group.by = 'FullData.category',  pt.size = 0)
save_plot("Plots/Figure4_Vln_by_category.png", p, base_width = 8 , base_height = 15) 
save_plot("Plots/Figure4_Vln_by_category.pdf", p, base_width = 8 , base_height = 15) 

###CUT&RUN Analysis
#Combine E12 and E14 peaks into 1 file
closest_combined <- rbind(E12_closest, E14_closest)

#filter out peaks that are beyond + or - 5000 bp from TSS
closest_combined <- closest_combined %>% filter(Distance < abs(5000)) 
#Remove duplicated genes
closest_combined <- closest_combined[!duplicated(closest_combined$local_GeneName), ]

#Make list of just peak genes
Peak_list <- closest_combined$local_GeneName
#Intersect the list with genes that are present in the Seurat Object
Peak_list <- intersect(Peak_list, row.names(CombinedSCT))

save(Peak_list, file = "Filtered_Peak_List.Rdata")
write.csv(Peak_list, "Plots/Peak_list.csv" )

##### Calculate module score across different cell states ####
#Subset Object
CombinedSCT.small.mod <- subset(CombinedSCT, downsample = 10000)

#Calculate module score
CombinedSCT.small.mod <- AddModuleScore(CombinedSCT.small.mod, list(Peak_list))

#plot module scores
p <- FeaturePlot(CombinedSCT.small.mod, features = "Cluster1", pt.size = 0.5) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + NoAxes() +
  ggtitle("Module Score of \n ATOH1 targets") +
  theme(plot.title = element_text(size = 40, face = "bold"), legend.key.size = unit(1, "cm"), legend.text = element_text(size=16), legend.position = 'bottom')
save_plot("Plots/Figure4_ModuleScore.png", p, base_width = 8, base_height = 9)
save_plot("Plots/Figure4_ModuleScore.pdf", p, base_width = 8, base_height = 9)

##Determine significance of module score across cell state

#pull out data that contains module scores by cell type
module.df <- FetchData(CombinedSCT.small.mod, vars = c("Cluster1","celltype"))
module.df$celltype <- factor(module.df$celltype, levels = rev(c("Progenitor", "Intermediate Progenitor", "Migrating_1","Migrating_2","Mature","Non RL")))
#Run function to calculate Module Pvalues
ModPVals <- CalcModPvals(module.df)
write.csv(ModPVals, file= "Plots/ModulePvalues.csv")

p <- ModScoreRidgePlot(module.df) + 
  theme(axis.text = element_text(size = 20, face = 'bold'), axis.title = element_text(size = 30, face = 'bold'))

save_plot("Plots/Figure4_ModuleScore_RidgePlot.png", p, base_width = 8, base_height = 8)
save_plot("Plots/Figure4_ModuleScore_RidgePlot.pdf", p, base_width = 8, base_height = 8)

rm(CombinedSCT.small.mod)

##### Perform Enrichment on DEGs #####
CombinedSCT.small <- subset(CombinedSCT, downsample = 2000)
rm(CombinedSCT)
all_markers_enrichment_test <- FindAllMarkers(object = CombinedSCT.small,
                                                                   only.pos = FALSE,
                                                                   min.pct = 0.05,
                                                                   logfc.threshold = 0,
                                                                   return.thresh = 1)

save(all_markers_enrichment_test, file = "Plots/all_markers_enrichment_test.RData")

DEG_list <- all_markers_enrichment_test
CategoryList = c("Progenitor", "Intermediate Progenitor", "Migrating_1", "Migrating_2", "Mature", "Non RL")

cutoff <- 0.01

# We will only consider the upregulated genes
DEG_list <- DEG_list %>% 
  filter(avg_logFC > 0)

DEG_list[ , 'DEG'] <- NA

for (i in 1:length(DEG_list$p_val_adj)) {
  
  if (DEG_list$p_val_adj[i] < cutoff) {
    DEG_list$DEG[i] <- 1
  } else {
    DEG_list$DEG[i] <- 0
  }
}

DEG_list[ , 'Peak'] <- NA
for (i in 1:length(DEG_list$p_val_adj)) {
  
  if (DEG_list$gene[i] %in% Peak_list) {
    DEG_list$Peak[i] <- 1
  } else {
    DEG_list$Peak[i] <- 0
  }
}

p <- PlotEnrichment(DEG_list, CategoryList)
#Odd's Ratio will be printed
#[1] "0.905200877503554" "1.01129968814563"  "1.90189484506598"  "1.97411017643719"  "1.2808226301987"  
#[6] "1.14452646897198"

save_plot("Plots/Figure4_Odds_ratio_DEG_filterpos.png", p, base_width = 8, base_height = 8)
save_plot("Plots/Figure4_Odds_ratio_DEG_filterpos.pdf", p, base_width = 8, base_height = 8)

#### Fiji Statistical Analysis ####
##Run Mix Model Anova for Figure4H
#Quantification file
df <- Fiji_Quantification
#There are two random variables: ROI and replicate
#Different ROIs captured different part of the RL so the perscent overlap is variable.
#There could be variability due to different developmental timing between litters (replicates)
#Build the mixed model with genotype as fixed variable and ROI+replicate as random variables
model = lmer(number ~ genotype + (1|ROI) + (1|replicate), data = df, REML = FALSE)
#Null hypothesis is that genotype does not affect the number
null.hypothesis= lmer(number ~ (1|ROI) + (1|replicate), data = df, REML = FALSE)
#Test our model with null hypothesis
anova(model, null.hypothesis)
#Results
#Data: df
#Models:
#  null.hypothesis: number ~ (1 | ROI) + (1 | replicate)
#model: number ~ genotype + (1 | ROI) + (1 | replicate)
#npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
#null.hypothesis    4 158.72 162.28 -75.362   150.72                    
#model              5 160.07 164.52 -75.036   150.07 0.652  1     0.4194


