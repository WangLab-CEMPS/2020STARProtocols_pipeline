# Info ---------------------------
## Script name: DotplotSelect.R

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
GOResult <- read_csv("result/E50h_GM3D_DiffPeakGO.csv")


# filter and plot ---------------------------------------------------------

set.seed(19960203)
GOResult %>% 
  slice(sample(nrow(GOResult), 10)) %>% 
  group_by(Cluster) %>% 
  select(c(1,3,4,7)) %>% 
  mutate(Description = factor(Description, levels = Description),
         GeneRatio = eval(parse(text = GeneRatio))) %>% 
  ggplot(aes(x = Cluster, y = Description)) +
  theme_bw() +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_color_continuous(low="red", high="blue", 
                         guide=guide_colorbar(reverse=TRUE)) +
  scale_size(range=c(3, 8)) +
  ylab(NULL) -> p

pdf(paste0("plot/",experiment_name,"_GODotplot_Select.pdf"))
p
dev.off()
