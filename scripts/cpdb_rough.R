

plt_median<-plt %>% group_by(CellType_refined) %>% summarize(mean_umap1=median(UMAP_1), mean_umap2=median(UMAP_2))
plt_median<-as.data.frame(plt_median)

###############
## Ped Healthy Output
###############
means<-read.table(here("data/cellphonedb/statistical_analysis_significant_means_06_22_2023_09:47:59.txt"), sep="\t", header=T)
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
levels(differentialexp_plt$Cell1)<-c("CD3+ T-cells","Cholangiocytes" , "CLNK T-cells", "Cycling Myeloid",
                                     "Cycling T-cells","Doublet","Erythrocytes" ,"gd T-cells",  "Hepatocytes","HSC",
                                     "KC Like","Low_Quality", "LSEC", "Macrophage\n(CLEC9A high)",
                                     "Macrophage\n(MHCII high)","Mature B cells","Myeloid Erythrocytes\n(phagocytosis)",
                                     "Neutrophil", "NK-like cells","Plasma cells", "Platelets","RR Myeloid")
differentialexp_plt$Cell2<-as.factor(differentialexp_plt$Cell2)
levels(differentialexp_plt$Cell2)<-c("CD3+ T-cells","Cholangiocytes" , "CLNK T-cells", "Cycling Myeloid",
                                     "Cycling T-cells","Doublet","Erythrocytes" ,"gd T-cells",  "Hepatocytes","HSC",
                                     "KC Like","Low_Quality", "LSEC", "Macrophage\n(CLEC9A high)",
                                     "Macrophage\n(MHCII high)","Mature B cells","Myeloid Erythrocytes\n(phagocytosis)",
                                     "Neutrophil", "NK-like cells","Plasma cells", "Platelets","RR Myeloid")


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
  
interacting_UMAP("APP_CD74")
interacting_UMAP("CCL3_CCR5")




unique(differentialexp_plt_notself$interacting_pair)

pdf(file = here("figures/ped_cpdb_significant.pdf"), w=5, h=5,)
interacting_UMAP("CCL3_CCR5")
interacting_UMAP("CCL4_CCR5")
interacting_UMAP("CXCL2_CXCR2")
interacting_UMAP("CXCL8_CXCR2")
interacting_UMAP("CXCL3_CXCR2")
interacting_UMAP("CCL8_CCR1")
interacting_UMAP("CCL3L1_CCR1")
interacting_UMAP("CCL3_CCR1")
interacting_UMAP("CCL5_CCR1")
interacting_UMAP("PTPRC_MRC1")
interacting_UMAP("IL10_IL10_receptor")
interacting_UMAP("IL1B_IL1_receptor")
dev.off()


interacting_UMAP("CCL3_CCR5")
save_plts(interacting_UMAP("CCL3_CCR5"), "CCL3_CCR5_cpdb_pedhealthy", w=5,h=5)
interacting_UMAP("IL10_IL10_receptor")
save_plts(interacting_UMAP("IL1B_IL1_receptor"), "IL1B_IL1_receptor_cpdb_pedhealthy", w=5,h=5)
interacting_UMAP("CCL3_CCR1")
save_plts(interacting_UMAP("CCL3_CCR1"), "CCL3_CCR1_cpdb_pedhealthy", w=5,h=5)

