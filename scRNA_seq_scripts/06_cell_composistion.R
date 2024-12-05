
### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(colorspace)
library(cowplot)
library(RColorBrewer)




#####################
## healthy only
#####################
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))

cell_label_healthy<-cell_label[which(cell_label$Treatment=="Healthy"),]

cell_label_healthy$index<-rownames(cell_label_healthy)

cell_counts<-cell_label_healthy %>% 
  group_by(CellType_refined) %>% 
  summarise(count=length(unique(cell))) %>% 
  mutate(countT= sum(count)) %>%
  group_by(CellType_refined, add=TRUE) %>%
  mutate(per=100*count/countT)

healthy_count<-ggplot(cell_counts, aes(CellType_refined, count, fill=CellType_refined))+geom_bar(stat = "identity", color="black")+
  theme_bw()+th+fillscale_cellType+xlab("")+ylab("Cell Count") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
save_plts(healthy_count, "healthy_count", w=6,h=4)

healthy_percent<-ggplot(cell_counts, aes(CellType_refined, per, fill=CellType_refined))+geom_bar(stat = "identity", color="black")+
  theme_bw()+fillscale_cellType+xlab("")+ylab("Percent of Cells") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
save_plts(healthy_percent, "healthy_percent", w=6,h=4)

