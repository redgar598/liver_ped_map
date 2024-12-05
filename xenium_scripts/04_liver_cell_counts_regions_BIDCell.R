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

source("scripts/00_pretty_plots.R")

## the panel comes with cell type labels
load(file=here("data/cell_type_labels_BIDCell.RData"))






############
## cell counts
############

cell_counts<-plt_umap_xenium %>% 
  group_by(CellType) %>% 
  summarise(count=length(unique(cell))) %>% 
  mutate(countT= sum(count)) %>%
  group_by(CellType, add=TRUE) %>%
  mutate(per=100*count/countT)

## broad cell type
cell_counts$CellType_broad<-as.factor(cell_counts$CellType)
levels(cell_counts$CellType_broad)<-c("Immune\n&\nErythrocytes","Mesenchymal\n&\nCholangiocytes","Mesenchymal\n&\nCholangiocytes","Immune\n&\nErythrocytes","Immune\n&\nErythrocytes",
                                          "Immune\n&\nErythrocytes","Hepatocytes","Hepatocytes","Hepatocytes","Mesenchymal\n&\nCholangiocytes",
                                          "Mesenchymal\n&\nCholangiocytes","Immune\n&\nErythrocytes","Endothelial","Immune\n&\nErythrocytes","Immune\n&\nErythrocytes",
                                          "Immune\n&\nErythrocytes","Immune\n&\nErythrocytes","Immune\n&\nErythrocytes","Endothelial")

cell_counts$CellType<-factor(cell_counts$CellType, levels=c("Hepatocyte (Pericentral)", "Hepatocyte (Periportal)", "Hepatocyte (Cycling)",
                                                            "LSEC","VEC",
                                                            "HSC (Activated)","HSC (Periportal)","HSC (Quiescent)",
                                                            "Cholangiocytes" ,  "Cholangiocyte (Biliary)" ,
                                                            "CD3+ T-cells", "gd T-cells", "Mature B-cells","Plasma Cells",  "Neutrophil", 
                                                            "Mono-Mac","Macrophage MHCII High","KC Like", "Cycling Myeloid","Erythrocytes"    ))

xenium_count<-ggplot(cell_counts, aes(CellType, count, fill=CellType))+geom_bar(stat = "identity", color="black")+
  theme_bw()+fillscale_cellType+xlab("")+ylab("Cell Count") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
save_plts(xenium_count, "xenium_count", w=6,h=4)

xenium_percent<-ggplot(cell_counts, aes(CellType,  per, fill=CellType))+geom_bar(stat = "identity", color="black")+
  theme_bw()+fillscale_cellType+xlab("")+ylab("Percent of Cells") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
xenium_percent
save_plts(xenium_percent, "xenium_percent", w=6,h=4)


xenium_percent<-ggplot(cell_counts, aes(CellType,  per, fill=CellType))+geom_bar(stat = "identity", color="black")+
  theme_bw()+fillscale_cellType+xlab("")+ylab("Percent of Cells") + facet_grid(~CellType_broad, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
xenium_percent
save_plts(xenium_percent, "xenium_percent_groups", w=8,h=4)

############
## regions
############


samples<-c("output-XETG00082__0011329__tJWH049__20231003__204217",
           "output-XETG00082__0011329__tJWH054__20231003__204217",
           "output-XETG00082__0011329__tJWH050__20231003__204217",
           "output-XETG00082__0011333__tJWH052__20231003__204217",
           "output-XETG00082__0011333__tJWH051__20231003__204217",
           "output-XETG00082__0011333__tJWH055__20231003__204217")
models<-c("2024_02_20_15_44_00",
          "2024_02_25_19_43_04",
          "2024_02_25_14_51_00",
          "2024_02_25_22_53_17",
          "2024_02_25_19_18_28",
          "2024_02_26_01_58_17")





count_all_sample<-lapply(1:6, function(x){
  
  tiff_path<-paste("/media/redgar/Seagate Portable Drive/stroke_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
  smpl<-strsplit(samples[x],"__")[[1]][3]
  
  ##### centroids
  tiff_res <- readTIFF(tiff_path, as.is = TRUE)
  tiff_res <- reshape2::melt(tiff_res)
  tiff_res <- tiff_res[tiff_res$value != 0, ]
  
  colnames(tiff_res) <- c("coord_y", "coord_x", "cell_id")
  tiff_res <- tiff_res[, c("coord_x", "coord_y", "cell_id")]
  tiff_res$cell_id <- as.numeric(tiff_res$cell_id)
  tiff_res$coord_x <- as.numeric(tiff_res$coord_x)
  tiff_res$coord_y <- as.numeric(tiff_res$coord_y)
  tiff_res <- tiff_res[order(tiff_res$cell_id), ]
  rownames(tiff_res) <- NULL
  
  centroids<-as.data.frame(
    tiff_res %>%
      group_by(cell_id) %>%
      summarise(centroid_x = sum(coord_x) / length(coord_x), centroid_y = sum(coord_y) / length(coord_y)))
  
  plt_umap_xenium_sample<-plt_umap_xenium[which(plt_umap_xenium$sample == smpl),]
  plt_umap_xenium_sample$cell<-sapply(1:nrow(plt_umap_xenium_sample), function(x) strsplit(plt_umap_xenium_sample$cell[x],"_")[[1]][1])
  
  cell_centroid<-merge(plt_umap_xenium_sample, centroids, by.x="cell",by.y="cell_id")
  
  ## load shapes
  structure_annotation<-read.csv(paste(here("data/"),smpl, " Coordinates.csv",sep=""), skip=2)
  structure_annotation$Selection<-as.factor(structure_annotation$Selection)

  cell_centroid_stroke<-grab_cells("stroke_site",structure_annotation,cell_centroid)
  cell_centroid_not_stroke<-grab_cells("contralateral_stroke",structure_annotation,cell_centroid)
  
  cell_centroid_stroke$region<-"stroke_site"
  cell_centroid_not_stroke$region<-"contralateral_stroke"
  
  cell_centroid_regions_combined<-rbind(cell_centroid_not_stroke, cell_centroid_stroke)
  cell_centroid_regions_combined$sample<-smpl
  cell_centroid_regions_combined
  
})

count_all_sample<-do.call(rbind,count_all_sample)
 


count_plt<-as.data.frame(count_all_sample %>% dplyr::select(CellType,CellType_broad,region, sample) %>% group_by(region,sample) %>% count(CellType,CellType_broad))

# count_plt$CellType<-factor(count_plt$CellType, levels=c(
#   "VLMC","Astrocytes","TBD (astro CA2)" ,"Endothelial cells","Microglia / perivascular macrophages","Oligodendrocytes","Pericytes / smooth muscle cells", 
#   "L2 IT ENTl" ,"L2/3 IT CTX","L4/5 IT CTX",
#   "L5 NP CTX","L6 CT CTX", "L6 IT ENTl" ,
#   "Lamp5 interneurons","Meis2",
#   "Pvalb interneurons","Sst Chodl inteneurons",
#   "general_interneuron"))

control<-ggplot(count_plt[which((count_plt$sample%in%c("tJWH049","tJWH052"))),], aes(CellType_broad,n,  fill=CellType))+
  geom_bar(stat="identity", color="black")+
  facet_grid(sample~region)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Cell Count")+guides(fill=guide_legend(ncol=3))+
  theme(legend.position="bottom")

reprogrammed<-ggplot(count_plt[which((count_plt$sample%in%c("tJWH050","tJWH051"))),], aes(CellType_broad,n,  fill=CellType))+
  geom_bar(stat="identity", color="black")+
  facet_grid(sample~region)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Cell Count")+guides(fill=guide_legend(ncol=3))+
  theme(legend.position="bottom")

reprogrammed2<-ggplot(count_plt[which((count_plt$sample%in%c("tJWH054","tJWH055"))),], aes(CellType_broad,n,  fill=CellType))+
  geom_bar(stat="identity", color="black")+
  facet_grid(sample~region)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Cell Count")+guides(fill=guide_legend(ncol=3))+
  theme(legend.position="bottom")


gaba<-get_legend(ggplot(count_plt[which(count_plt$CellType_broad=="GABAergic"),], aes(CellType_broad,n,  fill=CellType))+
             geom_bar(stat="identity", color="black")+fillscale_cellType)

glut<-get_legend(ggplot(count_plt[which(count_plt$CellType_broad=="Glutamatergic"),], aes(CellType_broad,n,  fill=CellType))+
                   geom_bar(stat="identity", color="black")+fillscale_cellType)

glial<-get_legend(ggplot(count_plt[which(count_plt$CellType_broad=="Glial"),], aes(CellType_broad,n,  fill=CellType))+
                   geom_bar(stat="identity", color="black")+fillscale_cellType)



plot_grid(control + theme(legend.position = "none"),
          reprogrammed + theme(legend.position = "none"),
          plot_grid(glial, gaba, glut, ncol=3),
          reprogrammed2 + theme(legend.position = "none"))

save_plts(plot_grid(control + theme(legend.position = "none"),
                    reprogrammed + theme(legend.position = "none"),
                    plot_grid(glial, gaba, glut, ncol=3),
                    reprogrammed2 + theme(legend.position = "none")), "cell_count", w=15, h=15)







#################
## as a percent
#################


count_plt_percent<-as.data.frame(count_all_sample %>% 
                           dplyr::select(CellType,CellType_broad,region, sample) %>% 
                           group_by(region,sample) %>% 
                           count(CellType,CellType_broad) %>%
                           mutate(per =  100 *n/sum(n)))

count_plt_label<-as.data.frame(count_all_sample %>% 
                                   dplyr::select(CellType_broad,region, sample) %>% 
                                   group_by(region,sample) %>% 
                                   count(CellType_broad) %>%
                                   mutate(per =  100 *n/sum(n)))

# 
# count_plt_percent$CellType<-factor(count_plt_percent$CellType, levels=c(
#   "VLMC","Astrocytes","TBD (astro CA2)" ,"Endothelial cells","Microglia / perivascular macrophages","Oligodendrocytes","Pericytes / smooth muscle cells", 
#   "L2 IT ENTl" ,"L2/3 IT CTX","L4/5 IT CTX",
#   "L5 NP CTX","L6 CT CTX", "L6 IT ENTl" ,
#   "Lamp5 interneurons","Meis2",
#   "Pvalb interneurons","Sst Chodl inteneurons",
#   "general_interneuron"))

count_plt_percent$region<-as.factor(count_plt_percent$region)
levels(count_plt_percent$region)<-c("Contralateral site", "Stroke site" )
count_plt_label$region<-as.factor(count_plt_label$region)
levels(count_plt_label$region)<-c("Contralateral site", "Stroke site" )

count_plt_label$n<-sapply(1:nrow(count_plt_label), function(x) prettyNum(count_plt_label$n[x], big.mark = ",", scientific = FALSE))

control<-ggplot(count_plt_percent[which((count_plt_percent$sample%in%c("tJWH049","tJWH052"))),], aes(CellType_broad,per))+
  geom_bar(aes(fill=CellType),stat="identity", color="black")+
  facet_grid(sample~region)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Percent of Cells")+guides(fill=guide_legend(ncol=3))+
  geom_text(aes(label=n),count_plt_label[which((count_plt_label$sample%in%c("tJWH049","tJWH052"))),], vjust=-0.5)+
  theme(legend.position="bottom")+ylim(0,75)+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=15),
    axis.title.y = element_text(size=15))+ggtitle("Controls")
control

reprogrammed<-ggplot(count_plt_percent[which((count_plt_percent$sample%in%c("tJWH050","tJWH051"))),], aes(CellType_broad,per))+
  geom_bar(aes(fill=CellType),stat="identity", color="black")+
  facet_grid(sample~region)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Percent of Cells")+guides(fill=guide_legend(ncol=3))+
  geom_text(aes(label=n),count_plt_label[which((count_plt_label$sample%in%c("tJWH050","tJWH051"))),], vjust=-0.5)+
  theme(legend.position="bottom")+ylim(0,75)+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=15),
    axis.title.y = element_text(size=15))+ggtitle("Reprogrammed")

reprogrammed2<-ggplot(count_plt_percent[which((count_plt_percent$sample%in%c("tJWH054","tJWH055"))),], aes(CellType_broad,per))+
  geom_bar(aes(fill=CellType), stat="identity", color="black")+
  facet_grid(sample~region)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Percent of Cells")+guides(fill=guide_legend(ncol=3))+
  geom_text(aes(label=n),count_plt_label[which((count_plt_label$sample%in%c("tJWH054","tJWH055"))),], vjust=-0.5)+
  theme(legend.position="bottom")+ylim(0,75)+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=15),
    axis.title.y = element_text(size=15))


gaba<-get_legend(ggplot(count_plt_percent[which(count_plt_percent$CellType_broad=="GABAergic"),], aes(CellType_broad,n,  fill=CellType))+
                   geom_bar(stat="identity", color="black")+fillscale_cellType+guides(fill = guide_legend(title = "GABAergic\nCell Types"))+
                   theme(
                     legend.title=element_text(size=15),
                     legend.text=element_text(size=12)))

glut<-get_legend(ggplot(count_plt_percent[which(count_plt_percent$CellType_broad=="Glutamatergic"),], aes(CellType_broad,n,  fill=CellType))+
                   geom_bar(stat="identity", color="black")+fillscale_cellType+ guides(fill = guide_legend(title = "Glutamatergic\nCell Types"))+
                   theme(
                     legend.title=element_text(size=15),
                     legend.text=element_text(size=12))  )

glial<-get_legend(ggplot(count_plt_percent[which(count_plt_percent$CellType_broad=="Glial"),], aes(CellType_broad,n,  fill=CellType))+
                    geom_bar(stat="identity", color="black")+fillscale_cellType+ guides(fill = guide_legend(title = "Glial\nCell Types"))+
                    theme(
                      legend.title=element_text(size=15),
                      legend.text=element_text(size=12))  )



plot_grid(control + theme(legend.position = "none"),
          reprogrammed + theme(legend.position = "none"),
          plot_grid(glial, gaba, glut, ncol=3),
          reprogrammed2 + theme(legend.position = "none"))

save_plts(plot_grid(control + theme(legend.position = "none"),
                    reprogrammed + theme(legend.position = "none"),
                    plot_grid(glial,plot_grid(gaba, glut,  ncol=2),  ncol=1),
                    reprogrammed2 + theme(legend.position = "none")), "cell_count_percent", w=12, h=12)



################
## All Count Detail
################

all_counts<-ggplot(count_plt_percent, aes(reorder(CellType, CellType_broad),n))+
  geom_bar(aes(fill=CellType), stat="identity", color="black")+
  facet_grid(sample~region)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Cell Count")+
  geom_text(aes(label=n), vjust=-0.5)+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=15),
    axis.title.y = element_text(size=15))


save_plts(plot_grid(all_counts + theme(legend.position = "none"),
          plot_grid(glial, gaba, glut, ncol=1, align="l", axis="h"), rel_widths = c(3,1)), "cell_count_all_brokendown", w=15, h=12)


################
## percent difference
################
count_plt_percent$treatment<-as.factor(count_plt_percent$sample)
levels(count_plt_percent$treatment)<-c("Control" ,"Reprogrammed" ,"Reprogrammed", "Control" ,"Reprogrammed", "Reprogrammed")

glial_lines<-ggplot(count_plt_percent[which(count_plt_percent$CellType_broad=="Glial"),], aes(region,per))+
  geom_line(aes(group=sample), color="grey60")+
  geom_point(aes(fill=CellType), shape=21, size=2)+
  facet_grid(CellType~treatment)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Percent of Cells")+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=10),
    axis.title.y = element_text(size=15))+ggtitle("Glial")

Glutamatergic_lines<-ggplot(count_plt_percent[which(count_plt_percent$CellType_broad=="Glutamatergic"),], aes(region,per))+
  geom_line(aes(group=sample), color="grey60")+
  geom_point(aes(fill=CellType), shape=21, size=2)+
  facet_grid(CellType~treatment)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Percent of Cells")+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=10),
    axis.title.y = element_text(size=15))+ggtitle("Glutamatergic")

GABAergic_lines<-ggplot(count_plt_percent[which(count_plt_percent$CellType_broad=="GABAergic"),], aes(region,per))+
  geom_line(aes(group=sample), color="grey60")+
  geom_point(aes(fill=CellType), shape=21, size=2)+
  facet_grid(CellType~treatment)+fillscale_cellType+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Percent of Cells")+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=10),
    axis.title.y = element_text(size=15))+ggtitle("GABAergic")

plot_grid(glial_lines + theme(legend.position = "none"),
          GABAergic_lines + theme(legend.position = "none"),
          Glutamatergic_lines + theme(legend.position = "none"),
          plot_grid(glial, gaba, glut, ncol=1, align="l", axis="h"),ncol=4, rel_widths = c(1,1,1,1))


save_plts(plot_grid(glial_lines + theme(legend.position = "none"),
                    GABAergic_lines + theme(legend.position = "none"),
                    Glutamatergic_lines + theme(legend.position = "none"),
                    plot_grid(glial, gaba, glut, ncol=1, align="l", axis="h"),ncol=4, rel_widths = c(1,1,1,1)), 
          "cell_percentlines", w=15, h=12)
