## Load Libraries
library(here)
library(Seurat)

library(SCINA)
library(reshape2)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(tiff)

library(sp)
library(rgeos)
library(scales)
library(viridis)


source("xenium_scripts/00_pretty_plots.R")
source("xenium_scripts/00_long_functions.R")



load(file=here("data/cell_type_labels_BIDCell.RData"))


## broad cell type
plt_umap_xenium$CellType_broad<-as.factor(plt_umap_xenium$CellType)
levels(plt_umap_xenium$CellType_broad)<-c("Immune","Mesenchymal & Cholangiocytes","Mesenchymal & Cholangiocytes","Immune","Erythrocytes",
                                          "Immune","Hepatocytes","Hepatocytes","Hepatocytes","Mesenchymal & Cholangiocytes",
                                          "Mesenchymal & Cholangiocytes","Immune","Endothelial","Immune","Immune",
                                          "Immune","Immune","Immune","Endothelial")





count_plt<-as.data.frame(plt_umap_xenium %>% dplyr::select(CellType,CellType_broad, sample) %>% group_by(sample) %>% count(CellType,CellType_broad))

count_plt<-count_plt[which(count_plt$sample!="C94_2"),]

ggplot(count_plt, aes(CellType,n,  fill=CellType))+
  geom_bar(stat="identity", color="black")+
  fillscale_cellType+theme_bw()+facet_grid(.~CellType_broad, scales="free_x", space="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Cell Count")+guides(fill=guide_legend(ncol=3))+
  theme(legend.position="bottom")




#################
## as a percent
#################


count_plt_percentage <- count_plt %>%
  group_by(sample, CellType,CellType_broad) %>%
  summarise(count = sum(n), .groups = 'drop') %>%
  group_by(sample) %>%
  mutate(total_count = sum(count),
         percentage = (count / total_count) * 100) %>%
  ungroup()




ggplot(count_plt_percentage, aes(sample,percentage))+
  geom_bar(aes(fill=CellType),stat="identity", color="black")+
  facet_grid(.~CellType_broad, scales="free_y")+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Percent of Cells")+guides(fill=guide_legend(ncol=3))+
  theme(legend.position="bottom")+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=15),
    axis.title.y = element_text(size=15))



count_plot<-plot_grid(ggplot(count_plt_percentage[which(count_plt_percentage$CellType_broad=="Hepatocytes"),], aes(sample,percentage))+
            geom_bar(aes(fill=CellType),stat="identity", color="black")+guides(fill=guide_legend(ncol=1))+
            facet_grid(.~CellType_broad, scales="free_y")+fillscale_cellType+theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            xlab("")+ylab("Percent of Cells")+
            theme(legend.position="bottom")+ylim(0,65)+
            theme(
              axis.text = element_text(size=12),
              strip.text = element_text(size=15),
              axis.title.y = element_text(size=15)),
          ggplot(count_plt_percentage[which(count_plt_percentage$CellType_broad=="Mesenchymal & Cholangiocytes"),], aes(sample,percentage))+
            geom_bar(aes(fill=CellType),stat="identity", color="black")+guides(fill=guide_legend(ncol=2))+
            facet_grid(.~CellType_broad, scales="free_y")+fillscale_cellType+theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            xlab("")+ylab("Percent of Cells")+
            theme(legend.position="bottom")+ylim(0,65)+
            theme(
              axis.text = element_text(size=12),
              strip.text = element_text(size=15),
              axis.title.y = element_text(size=15)),
          ggplot(count_plt_percentage[which(count_plt_percentage$CellType_broad=="Endothelial"),], aes(sample,percentage))+
            geom_bar(aes(fill=CellType),stat="identity", color="black")+guides(fill=guide_legend(ncol=1))+
            facet_grid(.~CellType_broad, scales="free_y")+fillscale_cellType+theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            xlab("")+ylab("Percent of Cells")+
            theme(legend.position="bottom")+ylim(0,65)+
            theme(
              axis.text = element_text(size=12),
              strip.text = element_text(size=15),
              axis.title.y = element_text(size=15)),
          ggplot(count_plt_percentage[which(count_plt_percentage$CellType_broad=="Immune"),], aes(sample,percentage))+
            geom_bar(aes(fill=CellType),stat="identity", color="black")+guides(fill=guide_legend(ncol=2))+
            facet_grid(.~CellType_broad, scales="free_y")+fillscale_cellType+theme_bw()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            xlab("")+ylab("Percent of Cells")+
            theme(legend.position="bottom")+ylim(0,65)+
            theme(
              axis.text = element_text(size=12),
              strip.text = element_text(size=15),
              axis.title.y = element_text(size=15)), ncol=4, align = "h")

save_plts(count_plot, "cell_counts", h=6,w=16)
