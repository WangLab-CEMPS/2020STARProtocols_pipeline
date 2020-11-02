# Info ---------------------------
## Script name: Enrichment_clusterProfiler.R

## experiment_name:
experiment_name <- "E50h_GM3D"

## Purpose of script:
## 1. XX
## 2. XX

## Author: Guandong Shang
## Date Created: 2020-10-27
## Copyright (c) Guandong Shang, 2020
## Email: shangguandong1996@163.com

# Notes ---------------------------
## 1. XX
## 2. XX   
##

# Prepare -----------------------------------------------------------------

# load up the packages
library(dplyr)

library(clusterProfiler)
library(org.At.tair.db)

# Set Options
options(stringsAsFactors = F)

# load up the data
DiffResult <- readr::read_csv("result/E50h_GM3D_DiffPeakAnno.csv")

## set cutoff
Fold_cutoff <- 1
FDR_cutoff <- 0.05



# Data preparation -------------------------------------------------------------

DiffPeakGene <- list()

DiffResult %>% 
  filter(Fold > Fold_cutoff, 
         FDR < FDR_cutoff) %>% 
  pull(geneId) %>% 
  unique() -> DiffPeakGene[["up"]]

DiffResult %>% 
  filter(Fold < -Fold_cutoff, 
         FDR < FDR_cutoff) %>% 
  pull(geneId) %>% 
  unique() -> DiffPeakGene[["down"]]

str(DiffPeakGene)


# GO Analysis -------------------------------------------------------------

ego <- compareCluster(geneClusters = DiffPeakGene,
                      fun="enrichGO", 
                      OrgDb = org.At.tair.db,
                      keyType = "TAIR",
                      ont = "BP")
readr::write_csv(as_tibble(ego@compareClusterResult),
                 path = paste0("result/",experiment_name,"_DiffPeakGO.csv"))

pdf(paste0("plot/",experiment_name,"_GODotplot.pdf"))
dotplot(ego)
dev.off()
