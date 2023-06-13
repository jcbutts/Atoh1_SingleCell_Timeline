library(DESeq2)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

setwd("C:/Jessie/AllSamples/Publication/Figure7")

#DEGs on full dataset
load("C:/Jessie/AllSamples/101221/Subset/Plots/all_markers_50PC_08res_SCT.RData")

#Needed to switch from human -> mouse genes or mouse -> human genes
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "dec2021.archive.ensembl.org")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "dec2021.archive.ensembl.org" )

#MB Bulk datasets
G3_bulk <- readRDS("C:/Jessie/AllSamples/Publication/MB_test/Data/dds_g3_vs_all_subgroup.RDS")
G4_bulk <- readRDS("C:/Jessie/AllSamples/Publication/MB_test/Data/dds_g4_vs_all_subgroup.RDS")
Shh_bulk <- readRDS("C:/Jessie/AllSamples/Publication/MB_test/Data/dds_shh_vs_all_subgroup.RDS")

G3_results <- results(G3_bulk)
G4_results <- results(G4_bulk)
Shh_results <- results(Shh_bulk)

###Prepare bulk dataset gene lists for module score
#Remove ensemble names
rownames(G3_results) <- sub("__.*", "", rownames(G3_results))
rownames(G4_results) <- sub("__.*", "", rownames(G4_results))
rownames(Shh_results) <- sub("__.*", "", rownames(Shh_results))

### G3
G3_results <- as.data.frame(G3_results)
G3_results <- rownames_to_column(G3_results, var = "Human_gene")

# Define mouse genes for human genes
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = G3_results$Human_gene, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
colnames(genesV2)[1] <- "Human_gene"

# New dataframe with mouse genes annotated
G3_results_anno <- G3_results %>% left_join(genesV2, by = join_by(Human_gene))

### G4
G4_results <- as.data.frame(G4_results)
G4_results <- rownames_to_column(G4_results, var = "Human_gene")

# Define mouse genes for human genes
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = G4_results$Human_gene, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
colnames(genesV2)[1] <- "Human_gene"

# New dataframe with mouse genes annotated
G4_results_anno <- G4_results %>% left_join(genesV2, by = join_by(Human_gene))

### SHH
Shh_results <- as.data.frame(Shh_results)
Shh_results <- rownames_to_column(Shh_results, var = "Human_gene")

# Define mouse genes for human genes
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Shh_results$Human_gene, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
colnames(genesV2)[1] <- "Human_gene"

# New dataframe with mouse genes annotated
Shh_results_anno <- Shh_results %>% left_join(genesV2, by = join_by(Human_gene))

##### Enrichment of cb genes in MB samples #####
# Same code was repeated with MB_groups set for "G3_results_anno" and "Shh_results_anno"

MB_group <- Shh_results_anno
mouse_DEGs <- all_markers_SCT

clusters <- c("14","15","20" ,"23")

df.fisher <- as.data.frame(NA)
df.fisher[,c("Category","pval","conf.interval1","conf.interval2","odds_ratio")] <- NA
counter = 1

for (x in clusters) {
  
mouse_input <- filter(mouse_DEGs, cluster == x)

MB_enrichment <- MB_group %>% mutate(mouse_markers = ifelse(MGI.symbol %in% mouse_input$gene, "Yes", "No"))
MB_enrichment <- MB_enrichment %>% mutate(sig_gene = ifelse(padj < 0.01 & log2FoldChange > 0.25, "Yes", "No"))

MB_enrichment_filtered <- subset(MB_enrichment, log2FoldChange > 0)
MB_enrichment_filtered <- subset(MB_enrichment_filtered, !grepl("^RP", Human_gene))
MB_enrichment_filtered <- subset(MB_enrichment_filtered, !MGI.symbol == "NA")

contingency_table <- table(MB_enrichment_filtered$sig_gene, MB_enrichment_filtered$mouse_markers)
df.temp <- fisher.test(contingency_table)

df.temp2 <- c(0, x, df.temp$p.value,df.temp$conf.int[[1]],df.temp$conf.int[[2]],df.temp$estimate)
df.fisher[counter,] <- df.temp2

counter = counter + 1
}

df.fisher <- df.fisher %>% mutate(pval_adj = as.numeric(pval) * length(clusters))
OR <- df.fisher$odds_ratio
print(OR)
#Plot the enrichment
cluster_col = c( "seagreen3", "gray67", "darkviolet", "mediumorchid2")
cluster_order <- c("15", "14", "23", "20")

df.fisher$Category <- factor(df.fisher$Category, level = cluster_order)

p3 <- ggplot(df.fisher, aes(x = Category, y = as.numeric(odds_ratio))) + geom_point(color = cluster_col, size = 5) +
  geom_errorbar(aes(ymin=as.numeric(conf.interval1), ymax=as.numeric(conf.interval2)), color = cluster_col, width = 0.3) +
  geom_hline(yintercept=1,linetype=2) +
  geom_text(aes(label = paste0("p=", signif(as.numeric(pval_adj), digits = 2)), y = 6), color = cluster_col,
            position = position_dodge(width = 0.5), size = 16 / .pt, hjust = 0.5, show.legend = F) +
  labs(y = 'Enrichment of SHH MB sig.\n genes in cerebellar clusters') + 
  theme_bw() + theme(text = element_text(size = 20, face = 'bold') ) + 
  scale_x_discrete(labels = c("Cb Prog", "UBC", "Early EGL", "Late EGL")) 
                        
save_plot("Plots/SHH_enrichment.pdf", p3, base_width = 7, base_height = 8) 

p <- ggarrange(p3, p1, p2, ncol = 3)
save_plot("Plots/Figure7_combo_enrichment.pdf", p, base_width = 17, base_height = 4) 

#OR, order is 14, 15, 20, 23

#SHH OR
# "1.38981529204799" "1.34734275536393" "1.70613295492597" "2.37917661538458"

#G4 OR
#"3.14295050337209"  "0.968308822616244" "2.81596386575252"  "1.6661908126669" 

#G3 OR
#"1.49799833813303"  "0.446044990060667" "1.25017204999785"  "0.705022338369173"
