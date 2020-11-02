# Info ---------------------------
## Script name: DA_DiffBind.R

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
library(DiffBind)
library(BiocParallel)

library(dplyr)

library(TxDb.Athaliana.BioMart.plantsmart28)
library(ChIPseeker)

# Set Options
options(stringsAsFactors = F)

## set cutoff
Fold_cutoff <- 1
FDR_cutoff <- 0.05



# prepare sample info -----------------------------------------------------

tissue <- rep(c("E50h","GM3D"),each = 2)

sample_info <- data.frame(
  SampleID = paste(tissue, c("R1", "R2"), sep = "_"),
  
  Tissue = tissue,
  
  Replicate = 1:2,
  
  bamReads = list.files("rawdata/bam", pattern = "bam$", full.names = T),
  
  Peaks = list.files("rawdata/peak", full.names = T),
  
  PeakCaller = "narrow",
  
  stringsAsFactors = F
  
)  


# Step1 importing data ----------------------------------------------------
dba_meta <- dba(minOverlap = 1, sampleSheet = sample_info)
dba_meta

# Step2 counting reads ----------------------------------------------------
dba_count <- dba.count(dba_meta,minOverlap = 1)
dba_count

pdf(paste0("plot/",experiment_name,"_Sample_Cor.pdf"),
    width = 8,height = 7)
dba.plotHeatmap(dba_count,
                RowAttributes = DBA_TISSUE,
                ColAttributes = F)

dba.plotPCA(dba_count)
dev.off()



# DA ----------------------------------------------------------------------

dba_contrast <- dba.contrast(dba_count, 
                             group1 = dba_count$masks[["E50h"]],
                             group2 = dba_count$masks[["GM3D"]],
                             name1 = "E50h",
                             name2 = "GM3D")
dba_diff <- dba.analyze(dba_contrast)

pdf(paste0("plot/",experiment_name,"_MAplot.pdf"))
dba.plotMA(dba_diff,fold = Fold_cutoff,cex.main=0.8)
abline(h = c(-Fold_cutoff,Fold_cutoff),col = "#ec008c", lty = 5)
dev.off()


# Report&PeakAnno ---------------------------------------------------------

dba_report_all <- dba.report(dba_diff,th = 1)
dba_report_all$feature_id <- paste0(experiment_name,"_",names(dba_report_all))

seqlevels(dba_report_all) <- gsub("Chr","",seqlevels(dba_report_all))

peakAnno <- annotatePeak(dba_report_all,
                         TxDb = TxDb.Athaliana.BioMart.plantsmart28,
                         level = "gene")
peakAnno_tb <- as_tibble(peakAnno@anno)
peakAnno_tb %>% 
  mutate(seqnames = paste0("Chr", seqnames),
         geneChr = paste0("Chr", geneChr)) %>% 
  readr::write_csv(.,path = paste0("result/",experiment_name,"_DiffPeakAnno.csv"))
