#library(scales)
save_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/png/"),name,".png", sep=""), w=w, h=h)}


myColors_celltype <- c("#660cc7","#8642d1","#5612a3","#4b911d","#7a4202","#006d2c","#6bbce8","#3469ad","#d17906","#b01e87","#60ba5d","#207537",
                       "#a0c487","#d9a5a5","#87a4c9","#e8c392","#dea4ce","#79639a",
                       "#207537","#fa61ad","#b80783","#994676","#431039","#cb181d","maroon1","#67000d",
                       "grey","#ce1256","#a6d96a","#750c32","#d9667f")
color_possibilities_celltype<-c("B-cells","Mature B-cells","Plasma cells","CD3+ T-cells","Cholangiocytes","gd T-cells","Hepatocytes","HSC","LSEC","Myeloid cells","NK-like cells", "NK and T cells",
                                "NKT cells\n(Hepatocyte Like)","Cholangiocytes\n(Hepatocyte Like)","HSC\n(Hepatocyte Like)","LSEC\n(Hepatocyte Like)","Myeloid cells\n(Hepatocyte Like)","B-cells\n(Hepatocyte Like)",
                                "NKT cells","RR Myeloid","KC Like","Neutrophil","Neutrophil\n(DEFA+)","Erythrocytes","Mast cell","Myeloid Erythrocytes\n(phagocytosis)",
                                "Doublet","Macrophage\n(MHCII high)","Cycling T-cells","Macrophage\n(CLEC9A high)","Platelets")
names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)
colscale_cellType <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)


myColors_age <- c("#348595","#d6604d")
color_possibilities_age<-c( "Adult","Ped")
names(myColors_age) <- color_possibilities_age
fillscale_age <- scale_fill_manual(name="Age\nGroup",
                                         values = myColors_age, drop = T, limits=force)
colscale_age <- scale_color_manual(name="Age\nGroup",
                                         values = myColors_age, drop = T, limits=force)




th <-   theme(axis.text=element_text(size=10),
              axis.title=element_text(size=12),
              strip.text = element_text(size = 12),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14))

th_present <- theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    strip.text = element_text(size = 12),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14))



#################
## Fetal
#################
save_fetal_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/fetal/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/fetal/jpeg/"),name,".jpeg", sep=""), w=w, h=h)}


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
