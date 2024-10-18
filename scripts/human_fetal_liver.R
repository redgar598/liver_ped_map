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



myColors_celltype_fetal <- c("#660cc7","#8642d1","#5612a3","#4b911d","#7a4202","#006d2c","#6bbce8","#3469ad","#d17906","#b01e87","#60ba5d","#207537",
                             "#a0c487","#d9a5a5","#87a4c9","#e8c392","#dea4ce","#79639a",
                             "#207537","#fa61ad","#b80783","#994676","#431039","#cb181d","#cb181d","maroon1","#67000d",
                             "grey","#ce1256","#a6d96a","#d9595c","#801f6f","#c98fbf","#360d2f","#360d14","#754166")
color_possibilities_celltype_fetal<-c("pro B cell","B cell","pre pro B cell ","CD3+ T-cells","Cholangiocytes","gd T-cells","Hepatocyte","Fibroblast","Endothelial cell",
                                      "DC2","NK", "NK and T cells",
                                      "ILC precursor","Cholangiocytes\n(Hepatocyte Like)",
                                      "HSC/MPP","LSEC\n(Hepatocyte Like)","Monocyte-DC precursor",
                                      "pre B cell","NKT cells","Mono-Mac","Kupffer Cell","Neutrophil-myeloid progenitor","Mono-NK",
                                      "Mid  Erythroid","Early Erythroid","Mast cell","MEMP","Doublet",
                                      "VCAM1+ Erythroblastic Island Macrophage","Early lymphoid/T lymphocyte",
                                      "Megakaryocyte","Monocyte","pDC precursor","Erythroblastic Island Macrophage","Late Erythroid","DC1")
names(myColors_celltype_fetal) <- color_possibilities_celltype_fetal
fillscale_cellType_fetal <- scale_fill_manual(name="Cell Type",
                                              values = myColors_celltype_fetal, drop = T, limits=force)
colscale_cellType_fetal <- scale_color_manual(name="Cell Type",
                                              values = myColors_celltype_fetal, drop = T, limits=force)


get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
#################


#################
## Load raw QC'ed data
#################
d10x_fetal<-readRDS(here("../../../projects/macparland/RE/PediatricAdult/processed_data","d10x_fetal_raw.rds"))
#d10x_fetal<-readRDS("/media/redgar/Seagate Portable Drive/fetal_liver/d10x_fetal_raw.rds")

Idents(d10x_fetal) <- "Cell.Labels"
DefaultAssay(d10x_fetal)<-"RNA"

all.markers<-FindAllMarkers(d10x_fetal, features = "MTM1",logfc.threshold = 0)


#all.markers<-FindAllMarkers(d10x_fetal,logfc.threshold = 0)

#all.markers[grep("MTM1", all.markers$gene),]

save(all.markers,  file="../../../projects/macparland/RE/random_side/Fetal_differential_MTM1.RData")

Mtm1_sig<-all.markers[grep("MTM1", all.markers$gene),][which(all.markers[grep("MTM1", all.markers$gene),]$p_val_adj<0.005),]
Mtm1_sig

d10x_fetal <- NormalizeData(d10x_fetal)
d10x_fetal <- FindVariableFeatures(d10x_fetal, selection.method = "vst", nfeatures = 2000)
d10x_fetal <- ScaleData(d10x_fetal)
d10x_fetal <- RunPCA(d10x_fetal, ndims.print = 1:10, nfeatures.print = 10)
d10x_fetal <- RunUMAP(d10x_fetal, dims = 1:30)

umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_fetal, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_fetal@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")



DefaultAssay(d10x_fetal)<-"RNA"
gene="MTM1"
gene_exp<-FetchData(d10x_fetal, vars=gene)
gene_exp$cell<-rownames(gene_exp)
plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')

save(plt_myeloid, file="../../../projects/macparland/RE/random_side/Fetal_UMAP_mtm1.RData")





######################
## transfered from server
######################
load(here("data/Fetal_differential_MTM1.RData"))

all.markers[grep("MTM1", all.markers$gene),]


load(here("data/Fetal_UMAP_mtm1.RData"))


len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=Cell.Labels),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType_fetal+theme_bw()+
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position = "bottom")+guides(fill=guide_legend(nrow=3,byrow=TRUE,override.aes = list(size=3)))
nice_legend<-get_leg(forlegned_plot)


fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=Cell.Labels),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_cellType_fetal+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.5,color="black",
           arrow = arrow(type = "closed", length = unit(4, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=10,hjust = 0.05),
                     axis.title.y = element_text(size=10,hjust = 0.05,angle = 90),
                     legend.position = "none")

fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(-0.5*len_y_bar), label=paste0("n = ",comma(nrow(plt_myeloid))), size=3)
plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,1))




gene="MTM1"
percentile=0.5
exp_limit<-quantile(plt_myeloid[, which(colnames(plt_myeloid)==gene)], percentile)
plt_myeloid$gene_exp_limited<-NA
over_limit<-which(plt_myeloid[, which(colnames(plt_myeloid)==gene)]>exp_limit)
plt_myeloid$gene_exp_limited[over_limit]<-plt_myeloid[over_limit, which(colnames(plt_myeloid)==gene)]
plt_myeloid<-plt_myeloid[rev(order(plt_myeloid$gene_exp_limited)),]

plt_myeloid<-rbind(plt_myeloid[which(is.na(plt_myeloid$gene_exp_limited)),],
                   plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),][(order(plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),]$gene_exp_limited)),])

umap_expression<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(color="black",size=1)+xlab("UMAP 1")+ylab("UMAP 2")+
  geom_point(aes(color=(gene_exp_limited)),size=0.5)+xlab("UMAP 1")+ylab("UMAP 2")+
  #facet_wrap(~age_condition,ncol=2)+
  #geom_text(aes(x = 7, y = -15, label=paste0("n = ",comma(CellCount))), cell_num_all)+
  scale_color_continuous_sequential(palette = "Viridis", rev=F,
                                    name=paste(gene, "\nExpression"),na.value = "grey80")+
  annotate("segment",
           x = arr$x, xend = arr$x + c(arr$x_len, 0),
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.5,color="black",
           arrow = arrow(type = "closed", length = unit(4, 'pt'))) +
  theme_void()+theme(legend.text=element_text(size=8),
                     legend.title=element_text(size=12),
                     legend.position = "bottom",
                     plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=10,hjust = 0.05),
                     axis.title.y = element_text(size=10,hjust = 0.05,angle = 90))

violin_plt<-ggplot()+ geom_point(aes(Cell.Labels, gene_exp_limited), plt_myeloid, size=0.5,position = position_jitter(w = 0.3, h = 0), color="grey")+
  geom_violin( aes(Cell.Labels, MTM1, fill=Cell.Labels), plt_myeloid)+xlab("")+ylab("Mtm1 Expression")+
  theme_bw()+fillscale_cellType_fetal+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12),
                                        axis.title.x = element_text(size=14),
                                        legend.position = "none")
#+  geom_text(aes(cluster, y=max(plt_myeloid$Mtm1)*1.1, label="*"),Mtm1_sig, size=6)





Mmt1<-plot_grid(plot_grid(plot_grid(fanciest_UMAP,umap_expression,align = "hv", axis="lr"),nice_legend, ncol=1,  rel_heights  = c(4,1)), violin_plt, ncol=1, rel_heights = c(1,0.8))
ggsave(Mmt1, file=here("figures/MTM1_liver_human_fetal_scRNAseq.jpeg"), w=20, h=18)



