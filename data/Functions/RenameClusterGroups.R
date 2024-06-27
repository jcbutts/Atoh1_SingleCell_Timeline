################################################################################
# Seurat-based pipeline for integrated single-cell dat from mult. genotypes    #
# Author: Ryan Dhindsa                                                         #
################################################################################

# Imports ----------------------------------------------------------------------
#library(Seurat)
#library(cowplot)
#library(data.table)
#library(wesanderson)


# Globals ----------------------------------------------------------------------
#config <- yaml.load_file("~/conf/config.yml")
#obj <- readRDS(config$output$CombinedSCT)

# Functions --------------------------------------------------------------------
RenameClusterGroups <- function(obj, cluster.names.file) {
  # Renames clusters
  # Args:
  #   obj: Seurat object
  #   cluster.names.file: path to txt file containing old and new cluster ids
  # Returns:
  #   obj: seurat object with new cluster identities
  
  # stash current identities
  obj[["orig.clusters"]] <- Idents(object = obj)
  
  cluster.names <- fread(cluster.names.file, 
                         header = T)
  
  current.cluster.ids <- cluster.names$old
  new.cluster.ids <- cluster.names$new
  
  names(new.cluster.ids) <- current.cluster.ids
  
  obj <- RenameIdents(obj, new.cluster.ids)
  
  # stash new identities in meta data
  obj[["celltype"]] <- obj@active.ident
 # obj[["celltype.genotype"]] <- paste0(obj@active.ident, "_", obj$genotype)
  
  #saveRDS(obj, obj.path)
  
  return(obj)
}

# Main -------------------------------------------------------------------------

