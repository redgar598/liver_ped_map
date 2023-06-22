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


source("scripts/00_pretty_plots.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/00_plot_gene_exp.R")

load(here("data","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

umap_mat<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)
meta<-d10x.combined@meta.data
meta$cell<-rownames(meta)
plt<-merge(meta, umap_mat, by="cell")

plt_mean<-plt %>% group_by(CellType_refined) %>% summarize(mean_umap1=mean(UMAP_1), mean_umap2=mean(UMAP_2))
plt_mean<-as.data.frame(plt_mean)

len_x_bar<-((range(plt$UMAP_1))[2]-(range(plt$UMAP_1))[1])/10
len_y_bar<-((range(plt$UMAP_2))[2]-(range(plt$UMAP_2))[1])/10
arr <- list(x = min(plt$UMAP_1), y = min(plt$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)


ggplot(plt_mean, aes(mean_umap1,mean_umap2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType_refined),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_cellType+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")+
  geom_label(aes(mean_umap1, mean_umap2, label=CellType_refined), size=1, color="black", fill="white")


###############
## Ped Healthy Output
###############
means<-read.table(here("data/cellphonedb/statistical_analysis_significant_means_06_20_2023_17:10:20.txt"), sep="\t", header=T)
means[1:5,1:10]

CCR_sig<-means[grep("CCR",means$interacting_pair), ]
lapply(1:nrow(CCR_sig), function(x) cbind(CCR_sig[x,c(1:12)],CCR_sig[x,13:496][which(!(is.na(CCR_sig[x,13:496])))]) )

CCL_sig<-means[grep("CCL3|CCL4",means$interacting_pair), ]
CCL_sig[,c(1:12)]

CCL_plt<-melt(CCL_sig, id=colnames(CCL_sig)[1:12])
CCL_plt$variable<-as.character(CCL_plt$variable)

CCL_plt$Cell1<-sapply(1:nrow(CCL_plt), function(x) strsplit(CCL_plt$variable[x], "[.]")[[1]][1])
CCL_plt$Cell2<-sapply(1:nrow(CCL_plt), function(x) strsplit(CCL_plt$variable[x], "[.]")[[1]][1])
CCL_plt<-CCL_plt[which(!(is.na(CCL_plt$value))),]

CCL_plt


# Sample data
points <- data.frame(
  x = c(1, 2, 3, 4, 5),    # x-coordinates of points
  y = c(3, 2, 4, 5, 1)     # y-coordinates of points
)

interactions <- data.frame(
  from = CCL_plt$Cell1,  # Index of the point where interaction starts
  to = CCL_plt$Cell2,    # Index of the point where interaction ends
  value = CCL_plt$value  # Interaction values
)

# Create a new data frame to store line coordinates
lines <- data.frame(
  xstart = points$x[interactions$from],
  xend = points$x[interactions$to],
  ystart = points$y[interactions$from],
  yend = points$y[interactions$to])

# Plotting the network
ggplot() +
  geom_curve(
    data = lines,
    aes(x = xstart, y = ystart, xend = xend, yend = yend, color = as.factor(order)),
    curvature = 0.2,
    size = interactions$value,
    alpha = 0.7,
    lineend = "round"
  ) +
  geom_point(
    data = points,
    aes(x, y),
    size = 3,
    color = "black"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  xlim(0, 6) +
  ylim(0, 6) +
  theme_minimal()