library(reticulate)
library(scater)
library(Seurat)
library(cowplot)
library(tidyverse)
library(zellkonverter)
library(dplyr)
library(tidyr)
library(MAST)

sessionInfo()

args <- commandArgs(trailingOnly=TRUE)

setwd("/Users/noamibrahimhayat/Documents/Aachen/projects/hf_single_cell/")

indir   = "/Users/noamibrahimhayat/Documents/Aachen/projects/hf_single_cell/data/"
outdir  = "/Users/noamibrahimhayat/Documents/Aachen/projects/hf_single_cell/results/"

filename = "main_cluster_T-cells.h5ad"
filename = "main_cluster_cardiomyocyte.h5ad"
filename = "main_cluster_fibroblast.h5ad"	
filename = "main_cluster_macrophages.h5ad"		
filename = "main_cluster_neuronal.h5ad"		
filename = "main_cluster_vSMCs.h5ad"
filename = "main_cluster_adipocytes.h5ad"		
filename = "main_cluster_endothelial.h5ad"		
filename = "main_cluster_lymphatic_endo.h5ad"	
filename = "main_cluster_mast_cells.h5ad"		
filename = "main_cluster_pericyte.h5ad"


S_name <- paste0(outdir, filename)
S <- readH5AD(S_name)

S <- as.Seurat(S, counts = "counts", data = "counts") ## CHECK THIS

S <- NormalizeData(S, normalization.method = "LogNormalize", scale.factor = 10000)
S <- FindVariableFeatures(S, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(S)
S <- ScaleData(S, features = all.genes)

S$subcluster_name %>% unique()
S$treatment %>% unique()
#S$subcluster_name <- S$subcluster_name %>% as.character() %>% strsplit(., ": ") %>% lapply(., function(x) x[2]) %>% unlist()

## make subset of object 
Idents(object = S) <-  "subcluster_name"
S$subcluster_name <- S$subcluster_name %>% as.character()

## get DEG by MAST
celltypes <- S$subcluster_name %>% unique()
cond1     <- c("ADKPD", "control")

#print(S$subcluster_name)
Idents(object = S) <-  "subcluster_name"

for (celltypes in unique(S$subcluster_name) ) {
  S_sub <- subset(x=S, idents=celltypes)
  
  Idents(object = S_sub) <-  "treatment"
  print(celltypes)
  
  try(S_sub.markers1 <- FindMarkers(S_sub, ident.1 = cond1[1], ident.2 = cond1[2], min.pct = 0.0, min.cells.group = 3, test.use="MAST"))
  
  celltypes <- str_replace(celltypes, "/", "_")
  
  write.table(S_sub.markers1, file=paste0( "./DEG_", celltypes, "_", cond1[1], ".vs.", cond1[2], ".txt"), row.names=TRUE, col.names=TRUE, sep="\t")
}


