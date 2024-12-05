#################
## Shiny plots
#################

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
library(scales)
library(viridis)

library(broom)

source("xenium_scripts/00_pretty_plots.R")
source("xenium_scripts/00_long_functions.R")



load(here("data","zonation_allcells.RData"))

zonation_scores_plt<-do.call(rbind, zonation_scores)


ggplot(zonation_scores_plt, aes(umap_1, umap_2))+geom_point()



## colored by cell type
len_x_bar<-((range(zonation_scores_plt$umap_1))[2]-(range(zonation_scores_plt$umap_1))[1])/10
len_y_bar<-((range(zonation_scores_plt$umap_2))[2]-(range(zonation_scores_plt$umap_2))[1])/10
arr <- list(x = min(zonation_scores_plt$umap_1)-2, y = min(zonation_scores_plt$umap_2)-2, x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(zonation_scores_plt, aes(umap_1,umap_2))+
  geom_point(aes(fill=CellType),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(zonation_scores_plt[sample(nrow(zonation_scores_plt)),], aes(umap_1,umap_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=CellType),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_cellType+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")+
  annotate("text",x = min(zonation_scores_plt$umap_1)+(0.95*len_x_bar)-1, y = min(zonation_scores_plt$umap_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(nrow(zonation_scores_plt))), size=2)
fanciest_UMAP
save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "UMAP_all_samples", w=8, h=6)

save_plts(nice_legend, "UMAP_all_samples_legend", w=2, h=5)



## colored by individual
len_x_bar<-((range(zonation_scores_plt$umap_1))[2]-(range(zonation_scores_plt$umap_1))[1])/10
len_y_bar<-((range(zonation_scores_plt$umap_2))[2]-(range(zonation_scores_plt$umap_2))[1])/10
arr <- list(x = min(zonation_scores_plt$umap_1)-2, y = min(zonation_scores_plt$umap_2)-2, x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(zonation_scores_plt, aes(umap_1,umap_2))+
  geom_point(aes(fill=sample),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_manual(values=c("#b10026","#fc4e2a","#feb24c","#084594","#6baed6","#41ab5d","#e7298a"))+
  theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)

fanciest_UMAP<-ggplot(zonation_scores_plt[sample(nrow(zonation_scores_plt)),], aes(umap_1,umap_2))+ #
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=sample),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_color_manual(values=c("#b10026","#fc4e2a","#feb24c","#084594","#6baed6","#41ab5d","#e7298a"))+
    annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),legend.position = "none")+#
  annotate("text",x = min(zonation_scores_plt$umap_1)+(0.95*len_x_bar)-1, y = min(zonation_scores_plt$umap_2)+(0.5*len_y_bar)-2, label=paste0("n = ",comma(nrow(zonation_scores_plt))), size=2)
fanciest_UMAP

save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "UMAP_all_samples_colorbysample", w=8, h=6)

save_plts(nice_legend, "UMAP_all_samples_colorbysample_legend", w=2, h=5)
