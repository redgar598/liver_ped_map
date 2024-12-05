### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(SoupX)
library(colorspace)
library(cowplot)
library(DropletQC)




source("scRNA_seq_scripts/00_pretty_plots.R")
source("scRNA_seq_scripts/00_entropy_d10x.R")
source("scRNA_seq_scripts/00_plot_gene_exp.R")
source("scRNA_seq_scripts/00_fanciest_UMAP.R")



####################################
## Fancy dot plot
####################################
load("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/IFALD_adult_ped_integrated_healthy_only.rds")
#d10x.combined_healthy <- NormalizeData(d10x.combined_healthy,scale.factor = 10000, normalization.method = "LogNormalize")

# ## add cell type labels
# load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
# cell_label$index<-rownames(cell_label)
# cell_label<-cell_label[match(colnames(d10x), cell_label$index),]
# identical(colnames(d10x), cell_label$index)
# 
# d10x <- AddMetaData(d10x, metadata = cell_label)

DefaultAssay(d10x.combined_healthy)<-"RNA"

## markers
Macrophage_genes<-c( "PTPRC", "CD74")
LEC_genes<-c("CALCRL","RAMP2")
Hepatocyte_genes<-c("ALB", "CYP3A4")
Cholangiocytes_genes<-c( "EPCAM", "KRT7")
HSCs_genes<-c( "IGFBP7",  "SPARC")
T_genes<-c("CD3D","CD8A")
NK_genes<-c("NKG7","CD7")
CDC1_genes<-c("CLEC9A","XCR1")
gd_genes<-c("GNLY")
RBC<-c("HBA1","FCGR3A")
MAST<-c("TPSAB1", "AREG")
recent_recruit_myeloid<-c("S100A8","S100A9","CD68","LYZ")
kuffer_signature<-c("MARCO","VSIG4","CD5L")
neutro_gene<-c("CSF3R","FCGR3B")
MHCII<-c("HLA-DRA","HLA-DPB1")
b_genes_noIG<-c("MS4A1", "CD79B")
immunoglobins<-c("IGKC","IGHG1")
platelet_genes<-c("PPBP","NRGN")
cycle_genes<-c("MKI67","TOP2A")



gene_exp<-FetchData(d10x.combined_healthy, vars=c(Macrophage_genes,LEC_genes,Hepatocyte_genes,Cholangiocytes_genes,HSCs_genes,T_genes,NK_genes,gd_genes,RBC,
                                 recent_recruit_myeloid, kuffer_signature, neutro_gene, MHCII,CDC1_genes, b_genes_noIG, immunoglobins, 
                                 platelet_genes,cycle_genes))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt<-merge(gene_exp, d10x.combined_healthy@meta.data, by="cell")


## summarize
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_summary<-plt %>% group_by(CellType_refined, variable) %>% 
  summarise(mn=mean(value), count=length(which(value>0)), percent_exp=(length(which(value>0))/length(value))*100)
plt_summary <- plt_summary %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(mn))
plt_summary<-as.data.frame(plt_summary)

# remove dots where 0 cell expressing marker
plt_summary<-plt_summary[(which(plt_summary$count>0)),]

plt_summary$variable<-factor(plt_summary$variable, levels=rev(c(b_genes_noIG, immunoglobins, T_genes,gd_genes,NK_genes, 
                                                            Macrophage_genes,recent_recruit_myeloid,kuffer_signature,MHCII, CDC1_genes,
                                                            RBC, neutro_gene, platelet_genes,
                                                            Cholangiocytes_genes,LEC_genes,HSCs_genes,Hepatocyte_genes,cycle_genes)))

plt_summary$CellType_refined<-factor(plt_summary$CellType_refined, levels=(c("Mature B-cells","Plasma cells","Cycling Plasma","CD3+ T-cells",
                                                                                "gd T-cells","NK-like cells",
                                                                                "Mono-Mac","KC Like","Macrophage\n(MHCII high)","CDC1","Cycling Myeloid",
                                                                                "Myeloid Erythrocytes\n(phagocytosis)","Erythrocytes","Neutrophil","Platelets",
                                                                                "Cholangiocytes","LSEC","HSC","Hepatocytes","Doublet" )))

fancy_dotplot<-plot_grid(
  ggplot(plt_summary, aes(CellType_refined, variable, color=scaled, size=percent_exp))+geom_point()+
    th+theme_classic()+
    scale_color_continuous_sequential(palette = "Oslo", rev=F, name="Scaled\nMean\nExpression")+
    scale_size(name="Percent\nCells\nExpressing")+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),axis.ticks.x = element_blank(),
          plot.margin = margin(0.25,0.25,0,0.25,"cm"))+
    geom_hline(yintercept = 34.5, color="grey70")+ 
    geom_hline(yintercept = 32.5, color="grey70")+ 
    geom_hline(yintercept = 29.5, color="grey70")+ 
    geom_hline(yintercept = 27.5, color="grey70")+ 
    geom_hline(yintercept = 23.5, color="grey70")+ 
    geom_hline(yintercept = 20.5, color="grey70")+ 
    geom_hline(yintercept = 18.5, color="grey70")+ 
    geom_hline(yintercept = 16.5, color="grey70")+
    geom_hline(yintercept = 14.5, color="grey70")+
    geom_hline(yintercept = 12.5, color="grey70")+
    geom_hline(yintercept = 10.5, color="grey70")+  
    geom_hline(yintercept = 8.5, color="grey70")+
    geom_hline(yintercept = 6.5, color="grey70")+
    geom_hline(yintercept = 4.5, color="grey70")+ 
    geom_hline(yintercept = 2.5, color="grey70")  ,
  ggplot(plt_summary, aes(CellType_refined, y=1, fill=CellType_refined))+geom_tile(color="black")+
    th+theme_classic()+fillscale_cellType+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),
          legend.position = "none",axis.line  = element_blank(),
          plot.margin = margin(0,0,1,1,"cm")),
  ncol=1, rel_heights = c(6,1.3), align = "v", axis="lr")
fancy_dotplot
save_plts(fancy_dotplot, "Healthy_only_dot_plot_celltype", w=6,h=10)



####################################
## individual gene expression
####################################
load("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/IFALD_adult_ped_integrated_healthy_only.rds")

allcell_keymakers<-plot_grid(plot_gene_UMAP(d10x.combined_healthy,"PTPRC", 0.9),
                             plot_gene_UMAP(d10x.combined_healthy,"ALB", 0.9))
allcell_keymakers
save_plts(allcell_keymakers, "Healthyonly_all_cell_key_markers", w=11,h=5)


allcell_keymakers<-plot_grid(plot_gene_UMAP(d10x.combined_healthy[,1:100],"PTPRC", 0.9),
                             plot_gene_UMAP(d10x.combined_healthy[,1:100],"ALB", 0.9))
allcell_keymakers
save_plts(allcell_keymakers, "Healthyonly_all_cell_key_markers_legend_save", w=11,h=5)


#####################################
## fancy UMAP individual
#####################################
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_integrated_refinedlabels_withDropletQC.rds"))

## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x.combined, reduction = "umap"))#
umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
meta_myeloid<-d10x.combined@meta.data
meta_myeloid$cell<-rownames(meta_myeloid)
plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")

len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)


forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
  geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_cellType+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)
#guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))

fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
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
                     legend.position = "none")
## cell count
cell_num_all<-as.data.frame(table(d10x.combined@meta.data$age_id))
colnames(cell_num_all)<-c("age_id","CellCount")
fanciest_UMAP <- fanciest_UMAP + facet_wrap(~age_id, ncol=4)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)


fancy_individual<-plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,1))
save_plts(fancy_individual, "IFALD_individual_fancy", w=20,h=18)



####################################
## Cell type counts in map
####################################
load(here("/media/redgar/Seagate Portable Drive/ped_map_update_feb2024/","IFALD_adult_ped_cellRefined_withDropletQC.rds"))
table(cell_label$CellType_refined, cell_label$age_condition)
table(cell_label$age_condition)


cell_label_immune<-cell_label[which(!(cell_label$CellType_refined %in% c("Cholangiocytes","LSEC","HSC","Hepatocytes","Doublet","Low Quality"))),]
table(cell_label_immune$CellType_refined, cell_label_immune$age_condition)

cell_label_immune_ped<-cell_label_immune[which(!(cell_label_immune$age_condition %in% c("Adult Healthy"))),]
table(cell_label_immune_ped$CellType_refined)
nrow(cell_label_immune_ped)

# mean in individual
cell_label_qc<-cell_label[which(!(cell_label$CellType_refined %in% c("Doublet","Low Quality"))),]
cell_label_all_ped<-cell_label_qc[which(!(cell_label_qc$age_condition %in% c("Adult Healthy"))),]
table(cell_label_all_ped$CellType_refined)


means<-table(cell_label_immune_ped$individual)/table(cell_label_all_ped$individual)


means<-table(cell_label_immune$individual)/table(cell_label_qc$individual)

### Bar plot
max_count<-as.data.frame(cell_label %>% 
                           group_by(individual, CellType_refined,CellType_rough, Age) %>% 
                           summarise(n = n()))

max_count$individual_age<-sapply(1:nrow(max_count), function(x) paste(strsplit(max_count$individual[x], "_")[[1]][1],"\n(", max_count$Age[x],")", sep=""))

bar_individual<-ggplot() + 
  geom_bar(aes(fill=CellType_refined, y=n, x=reorder(individual_age, Age)),max_count, position="stack", stat="identity", color="black")+
  theme_bw()+th+ylab("Cell Count")+xlab("Individual")+ fillscale_cellType +facet_grid(CellType_rough~., scales="free_y")
bar_individual



### Bar plot immune
max_count<-as.data.frame(cell_label_immune %>% 
                           group_by(individual, CellType_refined,CellType_rough, Age,age_condition) %>% 
                           summarise(n = n()))

max_count$individual_age<-sapply(1:nrow(max_count), function(x) paste(strsplit(max_count$individual[x], "_")[[1]][1],"\n(", max_count$Age[x],")", sep=""))

bar_individual<-ggplot() + 
  geom_bar(aes(fill=CellType_refined, y=n, x=reorder(individual_age, Age)),max_count, position="stack", stat="identity", color="black")+
  theme_bw()+th+ylab("Cell Count")+xlab("Individual")+ fillscale_cellType 
bar_individual


max_count_percent<-as.data.frame(cell_label_immune %>%  group_by(individual,CellType_refined) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)))

max_count_percent<-merge(max_count_percent, max_count[,c("individual","Age","individual_age","CellType_refined","CellType_rough","age_condition")], by=c("individual","CellType_refined"))

bar_individual<-ggplot() + 
  geom_bar(aes(fill=CellType_refined, y=freq, x=reorder(individual_age, Age)),max_count_percent, position="stack", stat="identity", color="black")+
  theme_bw()+th+ylab("Cell Count")+xlab("Individual")+ fillscale_cellType +facet_grid(CellType_rough~., scales="free_y")
bar_individual

ggplot(max_count_percent, aes(age_condition, freq))+geom_boxplot()+geom_point()+facet_wrap(~CellType_refined, scales="free_y")




####################################
## Stimulation Data
####################################
stim<-read.csv("data/stimulation_data.csv")

t.test(IL1B.Positive..~Age, data=stim[which(stim$Gate=="CD11B_CD68"),])
t.test(IL1B.Positive..~Age, data=stim[which(stim$Gate=="CD14_CD68"),])

stderror <- function(x) sd(x)/sqrt(length(x))


mean_data<- stim %>% group_by(Age, Gate) %>% summarise(mn=mean(IL1B.Positive..), SEM = stderror(IL1B.Positive..))

ggplot()+geom_bar(aes(Age,IL1B.Positive.., fill=Age),stim, stat='summary', color="black")+
  geom_point(aes(Age,IL1B.Positive.., fill=Age),stim, shape=21)+
  facet_wrap(~Gate)+theme_bw()+scale_fill_manual(values = c("#D64A56","cornflowerblue"))+
  geom_errorbar(aes(Age, ymin = mn-SEM, ymax = mn+SEM), mean_data, width = 0.2)+
  ylab("% of CD11B+ CD68+ ")



CD14_CD68<-ggplot()+geom_bar(aes(Age,mn, fill=Age),mean_data[which(mean_data$Gate=="CD14_CD68"),], stat='identity', color="black")+
  geom_point(aes(Age,IL1B.Positive.., fill=Age),stim[which(stim$Gate=="CD14_CD68"),], shape=21)+
  theme_bw()+scale_fill_manual(values = c("#D64A56","cornflowerblue"))+
  geom_errorbar(aes(Age, ymin = mn-SEM, ymax = mn+SEM), mean_data[which(mean_data$Gate=="CD14_CD68"),], width = 0.2)+
  ylab("Percentage Myeloid Cells (CD68+ CD14+)\nSecreating IL1B")+theme(legend.position = "none")+th_present+
  annotate("text", x=1.5, y=45.25, label="*", size=8)+
  annotate("segment", x = 1, xend = 2, y = 45, yend = 45)+
  annotate("segment", x = 1, xend = 1, y = 45, yend = 44.5)+
  annotate("segment", x = 2, xend = 2, y = 45, yend = 44.5)
CD14_CD68
save_plts(CD14_CD68, "CD14_CD68_ICS", h=5,w=3)

CD11B_CD68<-ggplot()+geom_bar(aes(Age,mn, fill=Age),mean_data[which(mean_data$Gate=="CD11B_CD68"),], stat='identity', color="black")+
  geom_point(aes(Age,IL1B.Positive.., fill=Age),stim[which(stim$Gate=="CD11B_CD68"),], shape=21)+
  theme_bw()+scale_fill_manual(values = c("#D64A56","cornflowerblue"))+
  geom_errorbar(aes(Age, ymin = mn-SEM, ymax = mn+SEM), mean_data[which(mean_data$Gate=="CD11B_CD68"),], width = 0.2)+
  ylab("Percentage Myeloid Cells (CD68+ CD11B+)\nSecreating IL1B")+theme(legend.position = "none")+th_present+
  annotate("text", x=1.5, y=46, label="ns", size=4)+
  annotate("segment", x = 1, xend = 2, y = 45, yend = 45)+
  annotate("segment", x = 1, xend = 1, y = 45, yend = 44.5)+
  annotate("segment", x = 2, xend = 2, y = 45, yend = 44.5)
CD11B_CD68
save_plts(CD11B_CD68, "CD11B_CD68_ICS", h=5,w=3)





####################################
## IF Data
####################################
IF<-read.csv("data/IF_percent_positive.csv")

wilcox.test(MD2_CD86~Condition, data=IF)

stderror <- function(x) sd(x)/sqrt(length(x))


mean_data<- IF %>% group_by(Condition) %>% summarise(mn=mean(MD2_CD86), SEM = stderror(MD2_CD86))

ggplot()+geom_bar(aes(Condition,MD2_CD86, fill=Condition),IF, stat='summary', color="black")+
  geom_point(aes(Condition,MD2_CD86, fill=Condition),IF, shape=21)+
  theme_bw()+scale_fill_manual(values = c("cornflowerblue","#B4EB65"))+
  geom_errorbar(aes(Condition, ymin = mn-SEM, ymax = mn+SEM), mean_data, width = 0.2)+
  ylab("% pos")

ggplot()+geom_bar(aes(Condition,MD2_noCD68, fill=Condition),IF, stat='summary', color="black")+
  geom_point(aes(Condition,MD2_noCD68, fill=Condition),IF, shape=21)+
  theme_bw()+scale_fill_manual(values = c("#B4EB65","cornflowerblue"))+
  geom_errorbar(aes(Condition, ymin = mn-SEM, ymax = mn+SEM), mean_data, width = 0.2)+
  ylab("% pos")

ggplot()+geom_bar(aes(Condition,noMD2_CD68, fill=Condition),IF, stat='summary', color="black")+
  geom_point(aes(Condition,noMD2_CD68, fill=Condition),IF, shape=21)+
  theme_bw()+scale_fill_manual(values = c("#B4EB65","cornflowerblue"))+
  geom_errorbar(aes(Condition, ymin = mn-SEM, ymax = mn+SEM), mean_data, width = 0.2)+
  ylab("% pos")



MD2_CD68<-ggplot()+geom_bar(aes(Condition,MD2_CD86, fill=Condition),IF, stat='summary', color="black")+
  geom_point(aes(Condition,MD2_CD86, fill=Condition),IF, shape=21, size=2)+
  theme_bw()+scale_fill_manual(values = c("cornflowerblue","#B4EB65"))+
  geom_errorbar(aes(Condition, ymin = mn-SEM, ymax = mn+SEM), mean_data, width = 0.2)+
  ylab("MD2+ Myeloid Cells (CD68+)\n (% of All Cells)")+theme(legend.position = "none")+th_present+
  annotate("text", x=1.5, y=26.25, label="*", size=8)+
  annotate("segment", x = 1, xend = 2, y = 26, yend = 26)+
  annotate("segment", x = 1, xend = 1, y = 26, yend = 25.5)+
  annotate("segment", x = 2, xend = 2, y = 26, yend = 25.5)
MD2_CD68
save_plts(MD2_CD68, "MD2_CD68_IF", h=5,w=3)
save_plts(MD2_CD68, "MD2_CD68_IF_narrow", h=4,w=2)
