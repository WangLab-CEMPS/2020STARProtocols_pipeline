# Info ---------------------------
## Script name: Volcanoplot.R

## experiment_name:
experiment_name <- "E50h_GM3D"

## Purpose of script:
## 1. XX
## 2. XX

## Author: Guandong Shang
## Date Created: 2020-10-28
## Copyright (c) Guandong Shang, 2020
## Email: shangguandong1996@163.com

# Notes ---------------------------
## 1. XX
## 2. XX   
##

# Prepare -----------------------------------------------------------------

# load up the packages
library(tidyverse)

# load up the data
DiffPeakResult <- read_csv("result/E50h_GM3D_DiffPeakAnno.csv") %>% 
  select(Fold, FDR, feature_id, geneId) 

# prepare Point Size & Color ------------------------------------------------------

DiffPeakResult %>% 
  mutate(Color = case_when(
    abs(Fold) <= 1 | FDR >= 0.05 ~ "No_Sig",
    Fold > 1 & FDR < 0.05 ~ "Up",
    Fold < -1 & FDR < 0.05 ~ "down"
  )) %>% 
  mutate(PointSize = case_when(
    abs(Fold) <=  2 | FDR >= 0.05 ~ "A",
    abs(Fold) >= 4 &  FDR < 0.05 ~ "C",
    abs(Fold) > 2 & FDR < 0.05 ~ "B"
  )) -> DiffPeakResult
  


# plot Volcano ------------------------------------------------------------

ggplot(DiffPeakResult, aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = Color, 
                 size = PointSize),alpha = 0.6) +
  scale_size_manual(values=c("A" = 1,
                             "B" = 3,
                             "C" = 5),
                    guide = F) +
  scale_color_manual(values = c("No_Sig" = "#a6a6a6",
                                "Up" = "#d73c31",
                                "down" = "#669cc7")) +
  geom_vline(xintercept = c(-1, 1), 
             color="grey40",
             linetype="longdash", lwd = 0.5) +
  geom_hline(yintercept = -log10(0.05),color="grey40",
             linetype="longdash", lwd = 0.5) +
  theme_bw() -> p

pdf(paste0("plot/",experiment_name,"_Volcano.pdf"))
p
dev.off()
