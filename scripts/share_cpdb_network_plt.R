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

load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

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





##############
## Differential expression with age (output of FindMarkers)
##############

## age differential
de_RR<-read.csv(here("data","differential_age_RR.csv"))
de_KC<-read.csv(here("data","differential_age_KC.csv"))
de_MHCII<-read.csv(here("data","differential_age_MHCII.csv"))


plt_median<-plt %>% group_by(CellType_refined) %>% summarize(mean_umap1=median(UMAP_1), mean_umap2=median(UMAP_2))
plt_median<-as.data.frame(plt_median)

###############
## Ped Healthy Output (output of cellphonedb)
###############
means<-read.table(here("data/cellphonedb/statistical_analysis_significant_means_03_05_2024_15:06:27.txt"), sep="\t", header=T)
means[1:5,1:10]

differentialexp_cpdbsig_a<-means[which(means$gene_a%in%de_KC$X), ]
differentialexp_cpdbsig_a$differential<-"adiff"
differentialexp_cpdbsig_b<-means[which(means$gene_b%in%de_KC$X), ]
differentialexp_cpdbsig_b$differential<-"bdiff"
differentialexp_cpdbsig<-rbind(differentialexp_cpdbsig_a, differentialexp_cpdbsig_b)

differentialexp_plt<-melt(differentialexp_cpdbsig, id=colnames(differentialexp_cpdbsig)[c(1:12,which(colnames(differentialexp_cpdbsig)=="differential"))])
differentialexp_plt$variable<-as.character(differentialexp_plt$variable)

differentialexp_plt$Cell1<-sapply(1:nrow(differentialexp_plt), function(x) strsplit(differentialexp_plt$variable[x], "[.]")[[1]][1])
differentialexp_plt$Cell2<-sapply(1:nrow(differentialexp_plt), function(x) strsplit(differentialexp_plt$variable[x], "[.]")[[1]][2])
differentialexp_plt<-differentialexp_plt[which(!(is.na(differentialexp_plt$value))),]



## fix cell labels
differentialexp_plt$Cell1<-as.factor(differentialexp_plt$Cell1)
levels(differentialexp_plt$Cell1)<-c("CD3+ T-cells","CDC1","Cholangiocytes", "Cycling Myeloid", "Cycling Plasma",
                                     "Doublet","Erythrocytes" ,"gd T-cells",  "Hepatocytes","HSC",
                                     "KC Like", "LSEC", 
                                     "Macrophage\n(MHCII high)","Mature B cells","Mono-Mac","Myeloid Erythrocytes\n(phagocytosis)",
                                     "Neutrophil", "NK-like cells","Plasma cells", "Platelets")
differentialexp_plt$Cell2<-as.factor(differentialexp_plt$Cell2)
levels(differentialexp_plt$Cell2)<-c("CD3+ T-cells","CDC1","Cholangiocytes", "Cycling Myeloid", "Cycling Plasma",
                                     "Doublet","Erythrocytes" ,"gd T-cells",  "Hepatocytes","HSC",
                                     "KC Like", "LSEC", 
                                     "Macrophage\n(MHCII high)","Mature B cells","Mono-Mac","Myeloid Erythrocytes\n(phagocytosis)",
                                     "Neutrophil", "NK-like cells","Plasma cells", "Platelets")


differentialexp_plt<-merge(differentialexp_plt,plt_median, by.x="Cell1", by.y="CellType_refined")
colnames(differentialexp_plt)[which(colnames(differentialexp_plt)%in%c("mean_umap1","mean_umap2"))]<-c("Cell1x","Cell1y")
differentialexp_plt<-merge(differentialexp_plt,plt_median, by.x="Cell2", by.y="CellType_refined")
colnames(differentialexp_plt)[which(colnames(differentialexp_plt)%in%c("mean_umap1","mean_umap2"))]<-c("Cell2x","Cell2y")


## Just sig in KC for now
#differentialexp_plt<-differentialexp_plt[grep("KC_Like",differentialexp_plt$variable),]
differentialexp_plta<-differentialexp_plt[which(differentialexp_plt$differential=="adiff" & differentialexp_plt$Cell1=="KC Like"),]
differentialexp_pltb<-differentialexp_plt[which(differentialexp_plt$differential=="bdiff" & differentialexp_plt$Cell2=="KC Like"),]
differentialexp_plt<-rbind(differentialexp_plta, differentialexp_pltb)



## self interactions
differentialexp_plt_self<-do.call(rbind,lapply(1:nrow(differentialexp_plt), function(x) if(differentialexp_plt$Cell1[x]==differentialexp_plt$Cell2[x]){differentialexp_plt[x,]}else{}))
differentialexp_plt_notself<-do.call(rbind,lapply(1:nrow(differentialexp_plt), function(x) if(differentialexp_plt$Cell1[x]==differentialexp_plt$Cell2[x]){}else{differentialexp_plt[x,]}))

unique(differentialexp_plt_notself$interacting_pair)




###############
#### UMAP/network colored by celltype
###############

interacting_UMAP<-function(interaction_pair){
  differentialexp_plt_notself_pair<-differentialexp_plt_notself[which(differentialexp_plt_notself$interacting_pair==interaction_pair),]
  
  umap_network<-ggplot() +   
    theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=8,hjust = 0.05),
                       axis.title.y = element_text(size=8,hjust = 0.05,angle = 90),
                       legend.position = "none")+xlab("UMAP 1")+ylab("UMAP 2")+
    geom_point(aes(UMAP_1,UMAP_2), data=plt, size = 0.6, colour= "black", stroke = 1)+
    geom_point(aes(UMAP_1,UMAP_2, color=CellType_refined), data=plt,size=0.5)+
    geom_rect(data=plt_median, mapping=aes(xmin=min(mean_umap1)*1.1, xmax=max(mean_umap1)*1.21, ymin=min(mean_umap2)*1.25, ymax=max(mean_umap2)*1.5), fill = "white", alpha=0.05)+
    geom_point(aes(mean_umap1,mean_umap2, fill=CellType_refined), data=plt_median,size=2, shape=21)+
    geom_curve(
      data = differentialexp_plt_notself_pair,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y, alpha=as.factor(differential)), 
      arrow = arrow(length = unit(0.01, "npc"),type = "closed"),
      color = "grey10",curvature = -0.3,
      lineend = "round") +  scale_alpha_manual(values=c(1,0.5)) +
    fillscale_cellType+colscale_cellType+
    annotate("segment", 
             x = arr$x, xend = arr$x + c(arr$x_len, 0), 
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(2, 'pt')))
  
  ligand<-strsplit(interaction_pair,"_")[[1]][1]
  receptor<-strsplit(interaction_pair,"_")[[1]][2]
  
  
  if(sum(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==2){
    umap_network+
      annotate("text", label=ligand, x = min(plt$UMAP_1)*1, y = max(plt$UMAP_2),
               color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("text", label=receptor, x = min(plt$UMAP_1)*1, y = max(plt$UMAP_2)*0.8,
               color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("segment", 
               x = min(plt$UMAP_1)*1.01, xend = min(plt$UMAP_1)*1.01 , 
               y = max(plt$UMAP_2)*0.95, yend = max(plt$UMAP_2)*0.85, size=0.25,color="black",
               arrow = arrow(type = "closed", length = unit(2, 'pt')))
  }else{if(which(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==1){
    umap_network+annotate("text", label=ligand, x = min(plt$UMAP_1)*1,y = max(plt$UMAP_2),
                          color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("text", label=receptor, x = min(plt$UMAP_1)*1, y = max(plt$UMAP_2)*0.8,color="black")+
      annotate("segment", 
               x = min(plt$UMAP_1)*1.01, xend = min(plt$UMAP_1)*1.01 , 
               y = max(plt$UMAP_2)*0.95, yend = max(plt$UMAP_2)*0.85, size=0.25,color="black",
               arrow = arrow(type = "closed", length = unit(2, 'pt')))} else {
                 umap_network+annotate("text", label=ligand, x = min(plt$UMAP_1)*1, y = max(plt$UMAP_2), color="black")+
                   annotate("text", label=receptor,  x = min(plt$UMAP_1)*1,  y = max(plt$UMAP_2)*0.8, color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
                   annotate("segment",  
                            x = min(plt$UMAP_1)*1.01, xend = min(plt$UMAP_1)*1.01 , 
                            y = max(plt$UMAP_2)*0.95, yend = max(plt$UMAP_2)*0.85, size=0.25,color="black",
                            arrow = arrow(type = "closed", length = unit(2, 'pt')))
               }}
}

interacting_UMAP("CCL3_CCR1")



###############
## Color by receptor ligand expression
###############
plot_gene_UMAP_2gene_network_notblend<-function(d10x, gene,interaction_pair, percentile){
  differentialexp_plt_notself_pair<-differentialexp_plt_notself[which(differentialexp_plt_notself$interacting_pair==interaction_pair),]
  
  ### plot individual genes
  DefaultAssay(d10x) <- "RNA"
  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
  
  cell_num_all<-as.data.frame(table(plt_myeloid$age_condition))
  colnames(cell_num_all)<-c("age_condition","CellCount")
  
  len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
  len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
  arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)
  
  
  gene_exp<-FetchData(d10x, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')
  
  
  ### non zero expression over the 90th percentile
  exp_limit1<-quantile(plt_myeloid[, which(colnames(plt_myeloid)==gene[1])], percentile)
  exp_limit2<-quantile(plt_myeloid[, which(colnames(plt_myeloid)==gene[2])], percentile)
  
  both_gene_over<-plt_myeloid[which(plt_myeloid[, which(colnames(plt_myeloid)==gene[1])]>=exp_limit1 &
                                      plt_myeloid[, which(colnames(plt_myeloid)==gene[1])]>0 &
                                      plt_myeloid[, which(colnames(plt_myeloid)==gene[2])]>=exp_limit2 &
                                      plt_myeloid[, which(colnames(plt_myeloid)==gene[2])]>0),]
  both_gene_over$color<-"Both Highly Expressed"
  
  gene1_over<-plt_myeloid[which(plt_myeloid[, which(colnames(plt_myeloid)==gene[1])]>=exp_limit1 &
                                  plt_myeloid[, which(colnames(plt_myeloid)==gene[1])]>0),]
  gene1_over<-gene1_over[which(!(gene1_over$cell%in%both_gene_over$cell)),]
  gene1_over$color<-paste(gene[1],"Highly Expressed")
  
  gene2_over<-plt_myeloid[which(plt_myeloid[, which(colnames(plt_myeloid)==gene[2])]>=exp_limit2 &
                                  plt_myeloid[, which(colnames(plt_myeloid)==gene[2])]>0),]
  gene2_over<-gene2_over[which(!(gene2_over$cell%in%both_gene_over$cell)),]
  gene2_over$color<-paste(gene[2],"Highly Expressed")
  
  plt_point_color<-rbind(both_gene_over, gene1_over, gene2_over)
  plt_point_color<-plt_point_color[order(plt_point_color$cell),]
  
  
  ## Network UMAP overlay
  umap_network<-ggplot()+
    geom_point(aes(UMAP_1,UMAP_2),plt_myeloid, size = 1, colour= "black", stroke = 1)+
    geom_point(aes(UMAP_1,UMAP_2),plt_myeloid, color="grey",size=0.75)+
    geom_point(aes(UMAP_1,UMAP_2, color=color),plt_point_color, size=0.5)+
    geom_rect(data=plt_median, mapping=aes(xmin=min(mean_umap1)*1.21, xmax=max(mean_umap1)*1.75, ymin=min(mean_umap2)*1.25, ymax=max(mean_umap2)*1.5), fill = "white", alpha=0.05)+
    xlab("UMAP 1")+ylab("UMAP 2")+
    scale_color_manual(values=c("#88a000","#b80783","#03008e"))+
    geom_point(aes(mean_umap1,mean_umap2, fill=CellType_refined), data=plt_median,size=4, shape=21, color="white")+
    geom_curve(
      data = differentialexp_plt_notself_pair,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y, alpha=as.factor(differential)), 
      arrow = arrow(length = unit(0.01, "npc"),type = "closed"),
      color = "white",curvature = -0.3,
      lineend = "round", size=1.15) + 
    geom_curve(
      data = differentialexp_plt_notself_pair,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y, alpha=as.factor(differential)), 
      arrow = arrow(length = unit(0.01, "npc"),type = "closed"),
      color = "grey10",curvature = -0.3,
      lineend = "round") +  scale_alpha_manual(values=c(1,0.5)) +
    fillscale_cellType+
    annotate("segment", 
             x = arr$x, xend = arr$x + c(arr$x_len, 0), 
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(4, 'pt'))) +
    theme_void()+theme(legend.text=element_text(size=6),
                       legend.title=element_text(size=8), 
                       plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=7,hjust = 0.05, vjust = 10),
                       axis.title.y = element_text(size=7,hjust = 0.05, vjust=-9,angle = 90),
                       legend.position = "none")
  
  
  ## Cell type node legend
  nice_legend<-get_leg(ggplot()+ 
                         geom_point(aes(mean_umap1,mean_umap2, fill=CellType_refined), 
                                    data=plt_median,size=3, shape=21, color="white")+
                         fillscale_cellType+theme_void()+
                         theme(legend.text=element_text(size=8),legend.title=element_text(size=10),
                               plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))+
                         guides(fill=guide_legend(ncol=2)))
  
  ## gene expression legend
  exp_legend<-get_leg(ggplot()+
                        geom_point(aes(UMAP_1,UMAP_2, color=color),plt_point_color, size=2, alpha=0.5)+
                        scale_color_manual(values=c("#88a000","#b80783","#03008e"), 
                                           name=paste("Gene Expression\n(Higher than ", percentile*100, "th\npercentile for gene",sep=""))+
                        theme(legend.text=element_text(size=8),legend.title=element_text(size=10),
                              plot.margin = margin(0.1,0.1,0.1,0.1, "cm")))
  
  
  ## Interaction Legend
  ligand<-strsplit(interaction_pair,"_")[[1]][1]
  receptor<-strsplit(interaction_pair,"_")[[1]][2]
  
  if(sum(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==2){
    interaction_legend<- ggplot()+
      annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
      annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
      annotate("text", label=ligand, x = 1, y = 1, color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("text", label=receptor, x = 2, y = 1,color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("segment", 
               x = 1.3, xend = 1.7 , 
               y = 1, yend = 1, size=0.25,color="black",
               arrow = arrow(type = "closed", length = unit(2, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
    
  }else{if(which(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==1){
    
    interaction_legend<-ggplot()+
      annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
      annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
      annotate("text", label=ligand, x = 1, y = 1, color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("text", label=receptor, x = 2, y = 1,color="black")+
      annotate("segment", 
               x = 1.3, xend = 1.7 , 
               y = 1, yend = 1, size=0.6,color="black",
               arrow = arrow(type = "closed", length = unit(5, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
  } else {
    
    if(sum(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==2){
      interaction_legend<- ggplot()+
        annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
        annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
        annotate("text", label=ligand, x = 1, y = 1, color="black")+
        annotate("text", label=receptor, x = 2, y = 1,color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
        annotate("segment", 
                 x = 1.3, xend = 1.7 , 
                 y = 1, yend = 1, size=0.25,color="black",
                 arrow = arrow(type = "closed", length = unit(2, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
    }}}
  
  
  plot_grid(umap_network,
            plot_grid(interaction_legend,nice_legend, exp_legend,ncol=1,rel_heights = c(0.2,1,0.5)),
            rel_widths = c(2,1))
  
}


#### with network as well as coloured by expression but not a blend
plot_gene_UMAP_2gene_network_notblend(d10x.combined,c("CCR5","CCL4"),"CCL4_CCR5",0.9)