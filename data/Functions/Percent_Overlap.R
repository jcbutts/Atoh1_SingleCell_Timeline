

# Function --------------------------------------------------------------------------
PercentPeakOverlap <- function(DEG_list){
E14_closest <-read.table('E14_Subset.score.narrowPeak_Closest', sep = '\t', 
                         col.names = c("chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue","peak","local_chrom","local_chromStart","local_chromEnd","local_GeneName","GeneID","plus","Distance"))
E12_closest <-read.table('E12_Subset.score.narrowPeak_Closest', sep = '\t', 
                         col.names = c("chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue","peak","local_chrom","local_chromStart","local_chromEnd","local_GeneName","GeneID","plus","Distance"))

#Combine E12 and E14 peaks into 1 file
closest_combined <- rbind(E12_closest, E14_closest)

#filter out absolute value
closest_combined <- closest_combined %>% 
  #filter(Distance <= 0) %>%  filter(Distance > -5000)
  filter(Distance < abs(5000)) 
#Remove duplicated genes
closest_combined <- closest_combined[!duplicated(closest_combined$local_GeneName), ]

DEG_list[ , 'Peak'] <- NA

for (i in 1:length(DEG_list$Genes)) {
  
  if (DEG_list$Genes[i] %in% closest_combined$local_GeneName  ) {
    DEG_list$Peak[i] <- 1
  } else {
    DEG_list$Peak[i] <- 0
  }
}


percent_overlap <- sum(DEG_list$Peak)/length(DEG_list$Gene)*100
print(paste('number of genes =', as.numeric(length(DEG_list$Gene))))
print(paste('Number of genes with peaks =', sum(DEG_list$Peak)))
print(paste('percent of genes with peaks = ', as.numeric(percent_overlap)))
}
