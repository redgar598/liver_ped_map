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
library(ggsignif)



source("scRNA_seq_scripts/00_pretty_plots.R")

myColors_condition <- c("#D64A56","cornflowerblue","#B4EB65",
                        "cornflowerblue","#B4EB65",
                        "#F4A261","#7D5BA6",
                        "#77baab","#d9c4f2","#42107d","#7ecae6")
names(myColors_condition) <- c( "Adult Healthy","Healthy","IFALD",
                                "Ped Healthy","Ped IFALD",
                                "Biliary Atresia", "Choledochal Cyst",
                                "ALGS","BASm","iBA","NC")
fillscale_condition <- scale_fill_manual(name="Age\nGroup",values = myColors_condition, drop = T, limits=force)
colscale_condition <- scale_color_manual(name="Age\nGroup",values = myColors_condition, drop = T, limits=force)

## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

###########################################################################################
#### Cell labelling
###########################################################################################

d10x_BA_taylor<-readRDS(here("BA_Taylor_d10x_adult_ped_integrated.rds"))

d10x_BA_taylor$condition[which(is.na(d10x_BA_taylor$condition))]<-d10x_BA_taylor$age_condition[which(is.na(d10x_BA_taylor$condition))]
d10x_BA_taylor$condition[which(is.na(d10x_BA_taylor$condition))]<-d10x_BA_taylor$Sample_characteristics_ch1[which(is.na(d10x_BA_taylor$condition))]

plot_grid(DimPlot(d10x_BA_taylor, group.by = "orig.ident"),
          DimPlot(d10x_BA_taylor, group.by = "CellType_refined")+colscale_cellType,
          DimPlot(d10x_BA_taylor, group.by = "condition"),
          DimPlot(d10x_BA_taylor, group.by = "seurat_clusters", label=T), ncol=2)


table(d10x_BA_taylor$seurat_clusters, d10x_BA_taylor$CellType_refined)
tbl <- table(d10x_BA_taylor$seurat_clusters, d10x_BA_taylor$CellType_refined)
max_col_per_row <- apply(tbl, 1, function(row) names(which.max(row)))
df<-as.data.frame(max_col_per_row)
df$cluster<-rownames(df)

d10x_BA_taylor@meta.data$CellType<-as.character(d10x_BA_taylor@meta.data$CellType_refined)
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("0"))]<-"NK-like cells"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("1","2","11"))]<-"CD3+ T-cells"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("3","34"))]<-"Mono-Mac"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("4"))]<-"Mature B-cells"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("5","8","26","33","35"))]<-"LSEC"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("6","10"))]<-"KC Like"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("7","16"))]<-"Macrophage\n(MHCII high)"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("9"))]<-"gd T-cells"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("12"))]<-"Cycling Myeloid"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("13","15","27","29","30"))]<-"Hepatocytes"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("14","23"))]<-"HSC"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("21"))]<-"CDC1"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("22"))]<-"Cholangiocytes"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("28","31","18"))]<-"Doublet"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("17"))]<-"Erythrocytes"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("25"))]<-"Myeloid Erythrocytes\n(phagocytosis)"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("20"))]<-"Neutrophil"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("19"))]<-"Plasma cells"
d10x_BA_taylor@meta.data$CellType[which(d10x_BA_taylor@meta.data$seurat_clusters%in%c("24"))]<-"pDC"

DimPlot(d10x_BA_taylor, group.by = "CellType")+colscale_cellType
DimPlot(d10x_BA_taylor, group.by = "CellType", split.by = "orig.ident")+colscale_cellType

#32, 33,35
FeaturePlot(d10x_BA_taylor, features=c("GZMB","IL3RA","PLAC8"))


table(d10x_BA_taylor$seurat_clusters, d10x_BA_taylor$condition)

table_data <- table(d10x_BA_taylor$seurat_clusters, d10x_BA_taylor$condition)
round(prop.table(table_data, margin = 2) * 100,2 ) # Convert to percentages


#Only see pDC in disease so missed in the other Hepatology paper
# TAylor map captured no cholangiocytes and no HSC


d10x_BA_taylor$study_condition<-as.factor(paste(d10x_BA_taylor$orig.ident, d10x_BA_taylor$condition))
levels(d10x_BA_taylor$study_condition)<-c("Xiao\nBiliary Atresia",  "Xiao\nCholedochal Cyst",
                                          "Ped Atlas\nAdult\nNon-Disease", "Ped Atlas\nPediatric\nNon-Disease","Ped Atlas\nIFALD" ,
                                          "Taylor\nALGS", "Taylor\nBASm","Taylor\niBA","Taylor\nNC"  )

DimPlot(d10x_BA_taylor, group.by = "CellType", split.by = "study_condition")+colscale_cellType


table(d10x_BA_taylor$CellType, d10x_BA_taylor$orig.ident)
table_data<-table(d10x_BA_taylor$CellType, d10x_BA_taylor$study_condition)
per_cell<-melt(prop.table(table_data, margin = 2) * 100)


baplor<-ggplot(per_cell[which(per_cell$Var1=="pDC"),], aes(Var2, value, fill=Var2))+geom_bar(stat="identity", color="black")+theme_bw()+
  scale_fill_manual(values=c("#238443", "#fed98e","grey","grey","#f768a1","#41b6c4","#78c679","#c2e699", "grey"))+
  xlab("Study and Condition")+ylab("pDCs as a Percentage of All Cells")+theme(legend.position = "none")
ggsave(baplor, file="bar_plot_pDCs.jpeg", width=10, height=6)



## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}



umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_BA_taylor, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_BA_taylor@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)


fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
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
                     legend.position = "none")


cell_num_all<-as.data.frame(table(d10x_BA_taylor@meta.data$study_condition))
colnames(cell_num_all)<-c("study_condition","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~study_condition, ncol=2)+  geom_text(aes(x = (min(plt_myeloid$UMAP_1)+(0.95*len_x_bar))-2, y = (min(plt_myeloid$UMAP_2)+(0.5*len_y_bar))-2, label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
ggsave(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), file="ped_public_integrated.jpeg", width=10, height=15)




#######################################
## DEG
#######################################

table(d10x_BA_taylor$condition)

DefaultAssay(d10x_BA_taylor)<-"RNA"

## Cholangiocytes
d10x_raw_Cholangiocytes<-subset(d10x_BA_taylor, subset = CellType %in% c("Cholangiocytes"))
d10x_raw_Cholangiocytes$condition<-as.factor(d10x_raw_Cholangiocytes$condition)
levels(d10x_raw_Cholangiocytes$condition)<-c( "Adult Healthy","Biliary Atresia", "Choledochal Cyst" ,"Ped Healthy" ,"Ped IFALD")
Idents(d10x_raw_Cholangiocytes)<-d10x_raw_Cholangiocytes$condition

table(d10x_raw_Cholangiocytes$condition)
table(d10x_raw_Cholangiocytes$condition, d10x_raw_Cholangiocytes$Sex)

saveRDS(d10x_raw_Cholangiocytes, file = here("~/BA_Taylor_d10x_adult_ped_cholangiocytes.rds"))


## HSC
d10x_raw_HSC<-subset(d10x_BA_taylor, subset = CellType %in% c("HSC"))
d10x_raw_HSC$condition<-as.factor(d10x_raw_HSC$condition)
levels(d10x_raw_HSC$condition)#<-c( "Adult Healthy","Biliary Atresia", "Choledochal Cyst" ,"Ped Healthy" ,"Ped IFALD")
Idents(d10x_raw_HSC)<-d10x_raw_HSC$condition

table(d10x_raw_HSC$condition)
table(d10x_raw_HSC$condition, d10x_raw_HSC$Sex)

saveRDS(d10x_raw_HSC, file = here("~/BA_Taylor_d10x_adult_ped_HSC.rds"))



## Myeloid
d10x_raw_myeloid<-subset(d10x_BA_taylor, subset = CellType %in% c("Mono-Mac","KC Like","Macrophage\n(MHCII high)","Cycling Myeloid","CDC1","Myeloid Erythrocytes\n(phagocytosis)"))
d10x_raw_myeloid$condition<-as.factor(d10x_raw_myeloid$condition)
levels(d10x_raw_myeloid$condition)<-c( "Adult Healthy","ALGS","BASm","Biliary Atresia", "Choledochal Cyst","iBA","NC"  ,"Ped Healthy" ,"Ped IFALD")
Idents(d10x_raw_myeloid)<-d10x_raw_myeloid$condition

table(d10x_raw_myeloid$condition)
table(d10x_raw_myeloid$condition, d10x_raw_myeloid$Sex)

saveRDS(d10x_raw_myeloid, file = here("~/BA_Taylor_d10x_adult_ped_myeloid.rds"))


###################
## Cholangiocytes DEG
###################
d10x_raw_Cholangiocytes<-readRDS(here("data/BA_Taylor_d10x_adult_ped_cholangiocytes.rds"))
d10x_raw_Cholangiocytes <- NormalizeData(d10x_raw_Cholangiocytes,scale.factor = 10000, normalization.method = "LogNormalize")


## IFALD differential
de_ifald<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_ifald<-de_ifald[which(de_ifald$p_val_adj < 0.005 & abs(de_ifald$avg_log2FC) > 1),]
sig_de_ifald[which(sig_de_ifald$avg_log2FC>0),]
sig_de_ifald[which(sig_de_ifald$avg_log2FC<0),]

## Choledochal cyst differential
de_cyst<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Choledochal Cyst", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_cyst<-de_cyst[which(de_cyst$p_val_adj < 0.005 & abs(de_cyst$avg_log2FC) > 1),]
sig_de_cyst[which(sig_de_cyst$avg_log2FC>0),]
sig_de_cyst[which(sig_de_cyst$avg_log2FC<0),]

## BA differential
de_BA<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Biliary Atresia", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_BA<-de_BA[which(de_BA$p_val_adj < 0.005 & abs(de_BA$avg_log2FC) > 1),]
sig_de_BA[which(sig_de_BA$avg_log2FC>0),]
sig_de_BA[which(sig_de_BA$avg_log2FC<0),]


## BA vs cyst
length(rownames(sig_de_cyst[which(sig_de_cyst$avg_log2FC>0),]))
length(rownames(sig_de_BA[which(sig_de_BA$avg_log2FC>0),]))
length(intersect(rownames(sig_de_BA[which(sig_de_BA$avg_log2FC>0),]), rownames(sig_de_cyst[which(sig_de_cyst$avg_log2FC>0),])))



### reproduce paper
## Choledochal cyst vs BA differential
de_cystBA<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Biliary Atresia", ident.2 = "Choledochal Cyst", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_cystBA<-de_cystBA[which(de_cyst$p_val_adj < 0.005 & abs(de_cyst$avg_log2FC) > 1),]
sig_de_cystBA[which(sig_de_cystBA$avg_log2FC>0),]
sig_de_cystBA[which(sig_de_cystBA$avg_log2FC<0),]

de_cystBA[grep("IL32|CCL2|IFITM3|STAT1|TNFRSF12A", rownames(de_cystBA)),]
VlnPlot(d10x_raw_Cholangiocytes, features = c("IL32","CCL2","IFITM3","STAT1","TNFRSF12A"))

## GSE176189 BA paper says "Differentially expressed genes for scRNA-seq data were defined as genes with fold change greater than 1.5 fold, p-value <0.05"
sig_de_cystBA_likepaper<-de_cystBA[which(de_cystBA$p_val_adj < 0.05 & abs(de_cystBA$avg_log2FC) > log2(1.5)),]
sig_de_cystBA_likepaper[grep("IL32|CCL2|IFITM3|STAT1|TNFRSF12A", rownames(sig_de_cystBA_likepaper)),]

sig_de_BA_likepaper<-de_BA[which(de_BA$p_val_adj < 0.05 & abs(de_BA$avg_log2FC) > log2(1.5)),]
sig_de_cyst_likepaper<-de_cyst[which(de_cyst$p_val_adj < 0.05 & abs(de_cyst$avg_log2FC) > log2(1.5)),]
sig_de_ifald_likepaper<-de_ifald[which(de_ifald$p_val_adj < 0.05 & abs(de_ifald$avg_log2FC) > log2(1.5)),]

length(rownames(sig_de_cystBA_likepaper[which(sig_de_cystBA_likepaper$avg_log2FC>0),]))
length(rownames(sig_de_BA_likepaper[which(sig_de_BA_likepaper$avg_log2FC>0),]))
length(intersect(rownames(sig_de_BA_likepaper[which(sig_de_BA_likepaper$avg_log2FC>0),]), rownames(sig_de_cystBA_likepaper[which(sig_de_cystBA_likepaper$avg_log2FC>0),])))


sig_de_BA_likepaper[grep("IL32|CCL2|IFITM3|STAT1|TNFRSF12A", rownames(sig_de_BA_likepaper)),]
de_BA[grep("IL32|CCL2|IFITM3|STAT1|TNFRSF12A", rownames(de_BA)),]

# in our comparison but not to cyst
unique_to_this_compare<-sig_de_BA_likepaper[which(!(rownames(sig_de_BA_likepaper)%in%rownames(sig_de_cystBA_likepaper))),]
unique_to_this_compare[rev(order(unique_to_this_compare$avg_log2FC)),]

## BA genes diff in cyst and healthy
sig_de_cyst_likepaper[grep("IL32|CCL2|IFITM3|STAT1|TNFRSF12A", rownames(sig_de_cyst_likepaper)),]
de_cyst[grep("IL32|CCL2|IFITM3|STAT1|TNFRSF12A", rownames(de_cyst)),]
## check many gene do not pass logfc.threshold = 0.25 in cyst vs healthy ped


## all comparisions for plot
de_cystBA<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Biliary Atresia", ident.2 = "Choledochal Cyst", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_cystBA_likepaper<-de_cystBA[which(de_cystBA$p_val_adj < 0.05 & abs(de_cystBA$avg_log2FC) > log2(1.5)),]

de_cystifald<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Ped IFALD", ident.2 = "Choledochal Cyst", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_cystiflad_likepaper<-de_cystifald[which(de_cystifald$p_val_adj < 0.05 & abs(de_cystifald$avg_log2FC) > log2(1.5)),]

de_BAifald<-FindMarkers(d10x_raw_Cholangiocytes, ident.1 = "Ped IFALD", ident.2 = "Biliary Atresia", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_BAiflad_likepaper<-de_BAifald[which(de_BAifald$p_val_adj < 0.05 & abs(de_BAifald$avg_log2FC) > log2(1.5)),]



## combine stats for plotting
sig_de_BA_likepaper$group1<-"Healthy"
sig_de_BA_likepaper$group2<-"Biliary Atresia"
sig_de_BA_likepaper$gene<-rownames(sig_de_BA_likepaper)
sig_de_cyst_likepaper$group1<-"Healthy"
sig_de_cyst_likepaper$group2<-"Choledochal Cyst"
sig_de_cyst_likepaper$gene<-rownames(sig_de_cyst_likepaper)
sig_de_ifald_likepaper$group1<-"Healthy"
sig_de_ifald_likepaper$group2<-"IFALD"
sig_de_ifald_likepaper$gene<-rownames(sig_de_ifald_likepaper)

sig_de_cystBA_likepaper$group1<-"Choledochal Cyst"
sig_de_cystBA_likepaper$group2<-"Biliary Atresia"
sig_de_cystBA_likepaper$gene<-rownames(sig_de_cystBA_likepaper)
sig_de_cystiflad_likepaper$group1<-"Choledochal Cyst"
sig_de_cystiflad_likepaper$group2<-"IFALD"
sig_de_cystiflad_likepaper$gene<-rownames(sig_de_cystiflad_likepaper)
sig_de_BAiflad_likepaper$group1<-"Biliary Atresia"
sig_de_BAiflad_likepaper$group2<-"IFALD"
sig_de_BAiflad_likepaper$gene<-rownames(sig_de_BAiflad_likepaper)


sig<-rbind(sig_de_BA_likepaper, sig_de_cyst_likepaper,sig_de_ifald_likepaper,
           sig_de_cystBA_likepaper,sig_de_cystiflad_likepaper,sig_de_BAiflad_likepaper)



### plot with relaxed statistics
DefaultAssay(d10x_raw_Cholangiocytes) <- "RNA"
meta_myeloid<-d10x_raw_Cholangiocytes@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)




plt_cho<-function(gene, vjust_sig){
  gene_exp<-FetchData(d10x_raw_Cholangiocytes, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  plt_myeloid<-merge(meta_myeloid, gene_exp, by='cell')
  colnames(plt_myeloid)[which(colnames(plt_myeloid)==gene)]<-"expression"
  plt_myeloid$condition<-as.factor(plt_myeloid$condition)
  levels(plt_myeloid$condition)<-c( "Adult Healthy","Biliary Atresia", "Choledochal Cyst","Healthy" ,"IFALD")
  plt_myeloid$condition<-factor(plt_myeloid$condition, levels = c( "Adult Healthy","Healthy" ,"IFALD","Biliary Atresia", "Choledochal Cyst"))
  
  sig_plt<-sig[which(sig$gene==gene),]
  sig_plt$y<-rev(seq_len(nrow(sig_plt))/10)
  
  ggplot(plt_myeloid[which(plt_myeloid$condition!="Adult Healthy"),], aes(condition,log(expression)))+
    geom_violin(fill="grey90",color="white")+geom_boxplot(aes(fill=condition),width=0.2)+fillscale_condition+
    theme_bw()+th_present+xlab("Pediatric Condition")+ylab(paste(gene, "Expression (log)"))+
    theme(legend.position = "none")+ ggtitle(gene)+
    geom_signif(
      xmin = (sig_plt$group1),
      xmax = (sig_plt$group2),
      y_position = max(log(plt_myeloid$expression))+sig_plt$y,
      annotation = "*",vjust=vjust_sig,
      tip_length = 0.01,textsize = 6)
}

plt_cho("IL32",0)
save_plts(plt_cho("IL32",0), "IL32_BA", w=7,h=7)
plt_cho("TNFRSF12A",0)
save_plts(plt_cho("TNFRSF12A",0), "TNFRSF12A_BA", w=7,h=7)


################ chol GSEA
source("scRNA_seq_scripts/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

### Age
de_BA$gene<-rownames(de_BA)
gene_list = de_BA$avg_log2FC
names(gene_list) = de_BA$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nBiliary Atresia"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Down-regulated in \nBiliary Atresia"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

cho_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=16.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))
cho_GSEA



### WHat might have been missed in comparing to cyst

### what was found between BA and ped that is not different between cyst and ped
BA_not_cyst<-sig_de_BA_likepaper[which(!(rownames(sig_de_BA_likepaper)%in%rownames(sig_de_cyst_likepaper))),]
cyst_not_BA<-sig_de_cyst_likepaper[which(!(rownames(sig_de_cyst_likepaper)%in%rownames(sig_de_BA_likepaper))),]
BA_and_cyst<-sig_de_BA_likepaper[which((rownames(sig_de_BA_likepaper)%in%rownames(sig_de_cyst_likepaper))),]

nrow(BA_not_cyst)
nrow(cyst_not_BA)
nrow(BA_and_cyst)





################ chol GSEA
source("scRNA_seq_scripts/00_GSEA_function.R")
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

### BA
de_BA$gene<-rownames(de_BA)
gene_list = de_BA$avg_log2FC
names(gene_list) = de_BA$gene
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Up-regulated in \nBiliary Atresia"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Down-regulated in \nBiliary Atresia"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
                                    
                                              
cho_GSEA<-ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  scale_fill_manual(values=c("#D64A56","cornflowerblue"))
cho_GSEA


BA_bar_GSEA<-ggplot(plt_path[1:10,], aes(NES, reorder(pathway, NES)))+geom_bar(fill="#F4A261",stat="identity", color="black")+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")
save_plts(BA_bar_GSEA, "BA_bar_GSEA", w=10,h=4)



myGO[grep("RESPONSE TO CYTOKINE", names(myGO))] # $`RESPONSE TO CYTOKINE%GOBP%GO:0034097`
myGO[grep("VIRAL PROCESS", names(myGO))] # $`VIRAL PROCESS%GOBP%GO:0016032`
myGO[grep("REGULATION OF I-KAPPAB KINASE/NF-KAPPAB SIGNALING", names(myGO))] # $`VIRAL PROCESS%GOBP%GO:0016032`



sapply(1:nrow(plt_path), function(x){
  print(paste(plt_path$pathway[x], ":", paste0(intersect(plt_path$leadingEdge[x][[1]], rownames(sig_de_cyst_likepaper)), collapse=",")))
})
                                              
#[1] "INTERFERON SIGNALING : HLA-DRA,HLA-DRB1,HLA-DPA1"


plt_cho("HLA-DRA",0.75)
save_plts(plt_cho("HLA-DRA",0.75), "HLA_BA", w=7,h=7)

plt_cho("CSNK2B",0.75)

sig_de_BA_likepaper[grep("HLA", sig_de_BA_likepaper$gene),]
sig_de_BA_likepaper[grep("CSNK2B", sig_de_BA_likepaper$gene),]


# we see and enchrich of infteferon signalling in out comparison of BA and healthy. Missed as HLA is high in cyst tooS

### new UMAP
DefaultAssay(d10x_raw_Cholangiocytes)<-"integrated"
d10x_raw_Cholangiocytes_pedonly_notaylor<-subset(d10x_raw_Cholangiocytes, subset = condition %in% c("Biliary Atresia", "Choledochal Cyst","Ped Healthy" ,"Ped IFALD"))

d10x_raw_Cholangiocytes_pedonly_notaylor <- RunPCA(d10x_raw_Cholangiocytes_pedonly_notaylor, npcs = 30, verbose = FALSE)
d10x_raw_Cholangiocytes_pedonly_notaylor <- RunUMAP(d10x_raw_Cholangiocytes_pedonly_notaylor, reduction = "pca", dims = 1:20)


umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_raw_Cholangiocytes_pedonly_notaylor, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_raw_Cholangiocytes_pedonly_notaylor@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

plt_myeloid$condition<-as.factor(as.character(plt_myeloid$condition))
levels(plt_myeloid$condition)<-c("Biliary Atresia", "Choledochal Cyst","Healthy" ,"IFALD")
plt_myeloid$condition<-factor(plt_myeloid$condition, levels = c( "Healthy" ,"IFALD","Biliary Atresia", "Choledochal Cyst"))


len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=condition),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_condition+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)



fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(size = 0.06, colour= "black", stroke = 1)+
  geom_point(aes(color=condition),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_condition+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")


fanciest_UMAP <- fanciest_UMAP + annotate("text",x = (min(plt_myeloid$UMAP_1)+(0.95*len_x_bar))-2, y = (min(plt_myeloid$UMAP_2)+(0.5*len_y_bar))-2, 
                                          label=paste0("n = ",comma(ncol(d10x_raw_Cholangiocytes_pedonly_notaylor))), size=2)

plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "ped_public_integrated_cholangiocyte", w=5, h=4.5)










###################
## HSC UMAP and DEG
###################
d10x_raw_HSC<-readRDS(here("data/BA_Taylor_d10x_adult_ped_HSC.rds"))
DefaultAssay(d10x_raw_HSC)<-"integrated"

d10x_raw_HSC <- RunPCA(d10x_raw_HSC, npcs = 30, verbose = FALSE)
d10x_raw_HSC <- RunUMAP(d10x_raw_HSC, reduction = "pca", dims = 1:20)
d10x_raw_HSC <- FindNeighbors(d10x_raw_HSC, reduction = "pca", dims = 1:30)
d10x_raw_HSC <- FindClusters(d10x_raw_HSC, resolution = 0.05)

DimPlot(d10x_raw_HSC, group.by = "seurat_clusters")
DimPlot(d10x_raw_HSC, group.by = "condition")+colscale_condition

table(d10x_raw_HSC$seurat_clusters, d10x_raw_HSC$condition)

d10x_raw_HSC@meta.data$CellType_HSC<-NA
d10x_raw_HSC@meta.data$CellType_HSC[which(d10x_raw_HSC@meta.data$seurat_clusters%in%c("0"))]<-"healthy_ped_HSC"
d10x_raw_HSC@meta.data$CellType_HSC[which(d10x_raw_HSC@meta.data$seurat_clusters%in%c("1"))]<-"adult_IFALD_HSC"
d10x_raw_HSC@meta.data$CellType_HSC[which(d10x_raw_HSC@meta.data$seurat_clusters%in%c("2"))]<-"Outlier HSC1"
d10x_raw_HSC@meta.data$CellType_HSC[which(d10x_raw_HSC@meta.data$seurat_clusters%in%c("3"))]<-"Outlier HSC2"

DimPlot(d10x_raw_HSC, group.by = "CellType_HSC")
FeaturePlot(d10x_raw_HSC, features = c("PDGFRA","CXCL12","COL1A1","IGFBP3"), min.cutoff = "q9", pt.size=0.25)



count_plt<-as.data.frame(table(d10x_raw_HSC@meta.data$CellType_HSC, d10x_raw_HSC$condition))
ggplot(count_plt, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity")

count_plt_percentage <- count_plt %>%
  group_by(Var1, Var2) %>%
  summarise(count = sum(Freq), .groups = 'drop') %>%
  group_by(Var1) %>%
  mutate(total_count = sum(count),
         percentage = (count / total_count) * 100) %>%
  ungroup()

count_percent<-ggplot(count_plt_percentage[which(count_plt_percentage$Var1!="Outlier HSC"),], aes(Var1,percentage))+
  geom_bar(aes(fill=Var2),stat="identity", color="black")+
  fillscale_condition+theme_bw()+
  xlab("")+ylab("Percent of Cells")+guides(fill=guide_legend(ncol=3))+coord_flip()+
  theme(
    axis.text = element_text(size=12),
    strip.text = element_text(size=15),
    axis.title.y = element_text(size=15))
count_percent

### DEG
DefaultAssay(d10x_raw_HSC)<-"RNA"

## IFALD differential
de_ifald<-FindMarkers(d10x_raw_HSC, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_ifald<-de_ifald[which(de_ifald$p_val_adj < 0.005 & abs(de_ifald$avg_log2FC) > 1),]
sig_de_ifald[which(sig_de_ifald$avg_log2FC>0),]
sig_de_ifald[which(sig_de_ifald$avg_log2FC<0),]

## Choledochal cyst differential
de_cyst<-FindMarkers(d10x_raw_HSC, ident.1 = "Choledochal Cyst", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_cyst<-de_cyst[which(de_cyst$p_val_adj < 0.005 & abs(de_cyst$avg_log2FC) > 1),]
sig_de_cyst[which(sig_de_cyst$avg_log2FC>0),]
sig_de_cyst[which(sig_de_cyst$avg_log2FC<0),]

## BA differential
de_BA<-FindMarkers(d10x_raw_HSC, ident.1 = "Biliary Atresia", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_BA<-de_BA[which(de_BA$p_val_adj < 0.005 & abs(de_BA$avg_log2FC) > 1),]
sig_de_BA[which(sig_de_BA$avg_log2FC>0),]
sig_de_BA[which(sig_de_BA$avg_log2FC<0),]

## BA vs cyst
length(rownames(sig_de_cyst[which(sig_de_cyst$avg_log2FC>0),]))
length(rownames(sig_de_BA[which(sig_de_BA$avg_log2FC>0),]))
length(intersect(rownames(sig_de_BA[which(sig_de_BA$avg_log2FC>0),]), rownames(sig_de_cyst[which(sig_de_cyst$avg_log2FC>0),])))





###################
## Myeloid UMAP and DEG
###################
d10x_raw_myeloid<-readRDS(here("data/BA_Taylor_d10x_adult_ped_myeloid.rds"))


### new UMAP
DefaultAssay(d10x_raw_myeloid)<-"integrated"

d10x_raw_myeloid <- RunPCA(d10x_raw_myeloid, npcs = 30, verbose = FALSE)
d10x_raw_myeloid <- RunUMAP(d10x_raw_myeloid, reduction = "pca", dims = 1:20)
d10x_raw_myeloid <- FindNeighbors(d10x_raw_myeloid, reduction = "pca", dims = 1:30)
d10x_raw_myeloid <- FindClusters(d10x_raw_myeloid, resolution = 0.05)

DimPlot(d10x_raw_myeloid, group.by = "seurat_clusters")
DimPlot(d10x_raw_myeloid, group.by = "CellType")+colscale_cellType
DimPlot(d10x_raw_myeloid, group.by = "CellType", split.by = "condition")+colscale_cellType

DimPlot(d10x_raw_myeloid, group.by = "condition")+colscale_condition

table(d10x_raw_myeloid$CellType, d10x_raw_myeloid$condition)

count_plt<-as.data.frame(table(d10x_raw_myeloid@meta.data$CellType, d10x_raw_myeloid$condition))
ggplot(count_plt, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="stack", stat="identity")+fillscale_condition



### UMAP just ped no taylor
d10x_raw_myeloid<-readRDS(here("data/BA_Taylor_d10x_adult_ped_myeloid.rds"))
DefaultAssay(d10x_raw_myeloid)<-"integrated"
d10x_raw_myeloid_pedonly_notaylor<-subset(d10x_raw_myeloid, subset = condition %in% c("Biliary Atresia", "Choledochal Cyst","Ped Healthy" ,"Ped IFALD"))
rm(d10x_raw_myeloid)
gc()

d10x_raw_myeloid_pedonly_notaylor <- RunPCA(d10x_raw_myeloid_pedonly_notaylor, npcs = 30, verbose = FALSE)
d10x_raw_myeloid_pedonly_notaylor <- RunUMAP(d10x_raw_myeloid_pedonly_notaylor, reduction = "pca", dims = 1:20)

umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x_raw_myeloid_pedonly_notaylor, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x_raw_myeloid_pedonly_notaylor@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

plt_myeloid$condition<-as.factor(as.character(plt_myeloid$condition))
levels(plt_myeloid$condition)<-c("Biliary Atresia", "Choledochal Cyst","Healthy" ,"IFALD")
plt_myeloid$condition<-factor(plt_myeloid$condition, levels = c( "Healthy" ,"IFALD","Biliary Atresia", "Choledochal Cyst"))


len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1)-2, y = min(plt_myeloid$UMAP_2)-2, x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)


fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
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
                     legend.position = "none")


cell_num_all<-as.data.frame(table(d10x_raw_myeloid_pedonly_notaylor@meta.data$condition))
colnames(cell_num_all)<-c("condition","CellCount")
cell_num_all<-cell_num_all[which(cell_num_all$condition%in%c("Biliary Atresia", "Choledochal Cyst","Ped Healthy" ,"Ped IFALD")),]

cell_num_all$condition<-as.factor(as.character(cell_num_all$condition))
levels(cell_num_all$condition)<-c("Biliary Atresia", "Choledochal Cyst","Healthy" ,"IFALD")
cell_num_all$condition<-factor(cell_num_all$condition, levels = c( "Healthy" ,"IFALD","Biliary Atresia", "Choledochal Cyst"))

fanciest_UMAP <- fanciest_UMAP + facet_wrap(~condition, ncol=2)+  geom_text(aes(x = (min(plt_myeloid$UMAP_1)+(0.95*len_x_bar))-2, y = (min(plt_myeloid$UMAP_2)+(0.5*len_y_bar))-2, label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)

plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
save_plts(plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2)), "ped_public_integrated_KC", w=5, h=4.5)







### DEG
DefaultAssay(d10x_raw_myeloid)<-"RNA"
Idents(d10x_raw_myeloid)<-"condition"

d10x_raw_myeloid <- NormalizeData(d10x_raw_myeloid,scale.factor = 10000, normalization.method = "LogNormalize")
d10x_raw_KC<-subset(d10x_raw_myeloid, subset = CellType %in% c("KC Like"))

## IFALD differential
de_ifald<-FindMarkers(d10x_raw_KC, ident.1 = "Ped IFALD", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_ifald<-de_ifald[which(de_ifald$p_val_adj < 0.005 & abs(de_ifald$avg_log2FC) > 1),]
sig_de_ifald[which(sig_de_ifald$avg_log2FC>0),]
sig_de_ifald[which(sig_de_ifald$avg_log2FC<0),]

## Choledochal cyst differential
de_cyst<-FindMarkers(d10x_raw_KC, ident.1 = "Choledochal Cyst", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_cyst<-de_cyst[which(de_cyst$p_val_adj < 0.005 & abs(de_cyst$avg_log2FC) > 1),]
sig_de_cyst[which(sig_de_cyst$avg_log2FC>0),]
sig_de_cyst[which(sig_de_cyst$avg_log2FC<0),]

## BA differential
de_BA<-FindMarkers(d10x_raw_KC, ident.1 = "Biliary Atresia", ident.2 = "Ped Healthy", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_BA<-de_BA[which(de_BA$p_val_adj < 0.005 & abs(de_BA$avg_log2FC) > 1),]
sig_de_BA[which(sig_de_BA$avg_log2FC>0),]
sig_de_BA[which(sig_de_BA$avg_log2FC<0),]

## BA vs cyst
# up reg
Cyst_up<-rownames(sig_de_cyst[which(sig_de_cyst$avg_log2FC>0),])
BA_up<-rownames(sig_de_BA[which(sig_de_BA$avg_log2FC>0),])
IFALD_up<-rownames(sig_de_ifald[which(sig_de_ifald$avg_log2FC>0),])

length(IFALD_up)
length(BA_up)
length(Cyst_up)

intersect(BA_up, Cyst_up)
intersect(IFALD_up, Cyst_up)
intersect(BA_up, IFALD_up)

intersect(IFALD_up,intersect(BA_up, Cyst_up))

#unique to
IFALD_up[which(!(IFALD_up%in%c(Cyst_up, BA_up)))]
BA_up[which(!(BA_up%in%c(Cyst_up, IFALD_up)))]
Cyst_up[which(!(Cyst_up%in%c(IFALD_up, BA_up)))]



## venn diagram
write.csv(IFALD_up, file="data/iflad_kc_genes.csv", quote=F)
write.csv(BA_up, file="data/BA_kc_genes.csv", quote=F)
write.csv(Cyst_up, file="data/cyst_kc_genes.csv", quote=F)


## disease comparisions
de_cystBA<-FindMarkers(d10x_raw_KC, ident.1 = "Biliary Atresia", ident.2 = "Choledochal Cyst", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_cystBA<-de_cystBA[which(de_cystBA$p_val_adj < 0.005 & abs(de_cystBA$avg_log2FC) > 1),]
de_cystiflad<-FindMarkers(d10x_raw_KC, ident.1 = "Choledochal Cyst", ident.2 = "Ped IFALD", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_cystiflad<-de_cystiflad[which(de_cystiflad$p_val_adj < 0.005 & abs(de_cystiflad$avg_log2FC) > 1),]
de_BAiflad<-FindMarkers(d10x_raw_KC, ident.1 = "Ped IFALD", ident.2 = "Biliary Atresia", test.use = "MAST",latent.vars=c("nFeature_RNA"), verbose=F)
sig_de_BAiflad<-de_BAiflad[which(de_BAiflad$p_val_adj < 0.005 & abs(de_BAiflad$avg_log2FC) > 1),]


## combine stats for plotting
sig_de_BA$group1<-"Healthy"
sig_de_BA$group2<-"Biliary Atresia"
sig_de_BA$gene<-rownames(sig_de_BA)
sig_de_cyst$group1<-"Healthy"
sig_de_cyst$group2<-"Choledochal Cyst"
sig_de_cyst$gene<-rownames(sig_de_cyst)
sig_de_ifald$group1<-"Healthy"
sig_de_ifald$group2<-"IFALD"
sig_de_ifald$gene<-rownames(sig_de_ifald)

sig_de_cystBA$group1<-"Choledochal Cyst"
sig_de_cystBA$group2<-"Biliary Atresia"
sig_de_cystBA$gene<-rownames(sig_de_cystBA)
sig_de_cystiflad$group1<-"Choledochal Cyst"
sig_de_cystiflad$group2<-"IFALD"
sig_de_cystiflad$gene<-rownames(sig_de_cystiflad)
sig_de_BAiflad$group1<-"Biliary Atresia"
sig_de_BAiflad$group2<-"IFALD"
sig_de_BAiflad$gene<-rownames(sig_de_BAiflad)

sig<-rbind(sig_de_BA, sig_de_cyst,sig_de_ifald,sig_de_cystBA,sig_de_cystiflad,sig_de_BAiflad)

### plot 
DefaultAssay(d10x_raw_KC) <- "RNA"
meta_myeloid<-d10x_raw_KC@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)


plt_kc<-function(gene, vjust_sig){
  gene_exp<-FetchData(d10x_raw_KC, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  plt_myeloid<-merge(meta_myeloid, gene_exp, by='cell')
  colnames(plt_myeloid)[which(colnames(plt_myeloid)==gene)]<-"expression"
  
  plt_myeloid<-plt_myeloid[which(plt_myeloid$condition%in%c("Biliary Atresia", "Choledochal Cyst","Ped Healthy" ,"Ped IFALD")),]
  plt_myeloid$condition<-as.factor(as.character(plt_myeloid$condition))
  levels(plt_myeloid$condition)<-c("Biliary Atresia", "Choledochal Cyst","Healthy" ,"IFALD")
  plt_myeloid$condition<-factor(plt_myeloid$condition, levels = c( "Healthy" ,"IFALD","Biliary Atresia", "Choledochal Cyst"))
  
  sig_plt<-sig[which(sig$gene==gene),]
  sig_plt$y<-rev(seq_len(nrow(sig_plt))/10)
  
  ggplot(plt_myeloid, aes(condition,log(expression)))+
    geom_violin(fill="grey90",color="white")+geom_boxplot(aes(fill=condition),width=0.2)+fillscale_condition+
    theme_bw()+th_present+xlab("Pediatric Condition")+ylab(paste(gene, "Expression (log)"))+
    theme(legend.position = "none")+ ggtitle(gene)+
    geom_signif(
      xmin = (sig_plt$group1),
      xmax = (sig_plt$group2),
      y_position = max(log(plt_myeloid$expression))+sig_plt$y,
      annotation = "*",vjust=vjust_sig,
      tip_length = 0.01,textsize = 6)
}

plt_kc("LY96",0)
save_plts(plt_kc("LY96",0), "LY96_KC", w=7,h=7)


plt_kc("SAA1",0)
