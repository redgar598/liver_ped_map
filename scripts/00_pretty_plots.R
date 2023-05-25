#library(scales)
save_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/png/"),name,".png", sep=""), w=w, h=h)}


myColors_celltype <- c("#660cc7","#8642d1","#5612a3","#4b911d","#7a4202",
                       "#2c6e02","#6bbce8","#3469ad","#d17906","#b01e87",
                       "#60ba5d","#207537","#a0c487","#d9a5a5","#87a4c9",
                       "#e8c392","#dea4ce","#79639a","#207537","#fa61ad",
                       "#b80783","#994676","#431039","#cb181d","maroon1",
                       "#b01629","grey","#ce1256","#a6d96a","#750c32",
                       "#d9667f","#1b4003","#e0a8ce","#8a68b0","#3d1b63",
                       "#c9a8ed","#c48db4","#a3588d","#6ca647","#3a7d31",
                       "#60ba5d","#3d1b63","#2dc918","#F4355B")
color_possibilities_celltype<-c("B-cells","Mature B-cells","Plasma cells","CD3+ T-cells","Cholangiocytes",
                                "gd T-cells","Hepatocytes","HSC","LSEC","Myeloid cells",
                                "NK-like cells", "NK and T cells","NKT cells\n(Hepatocyte Like)","Cholangiocytes\n(Hepatocyte Like)","HSC\n(Hepatocyte Like)",
                                "LSEC\n(Hepatocyte Like)","Myeloid cells\n(Hepatocyte Like)","B-cells\n(Hepatocyte Like)","NKT cells","RR Myeloid",
                                "KC Like","Neutrophil","Neutrophil\n(DEFA+)","Erythrocytes","Mast cell",
                                "Myeloid Erythrocytes\n(phagocytosis)","Doublet","Macrophage\n(MHCII high)","Cycling T-cells","Macrophage\n(CLEC9A high)",
                                "Platelets","CLNK T-cells","Cycling Myeloid","Mature B-cells (104)","pDC",
                                "pre B-cell","KC Like\n(Hepatocyte Like)","KC Like (C97)","Naive CD4 T-cells","Memory CD4 T-cells",
                                "NK cells","DC","CD8 T-cells","CD14+ Mono")
names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)
colscale_cellType <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)



myColors_age <- c("#D64A56","cornflowerblue")
color_possibilities_age<-c( "Adult","Ped")
names(myColors_age) <- color_possibilities_age
fillscale_age <- scale_fill_manual(name="Age\nGroup",
                                         values = myColors_age, drop = T, limits=force)
colscale_age <- scale_color_manual(name="Age\nGroup",
                                         values = myColors_age, drop = T, limits=force)

myColors_agecondition <- c("#D64A56","cornflowerblue","#374eb8")
color_possibilities_agecondition<-c( "Adult Healthy","Ped Healthy","Ped IFALD")
names(myColors_agecondition) <- color_possibilities_agecondition
fillscale_agecondition <- scale_fill_manual(name="Age\nGroup",
                                   values = myColors_agecondition, drop = T, limits=force)
colscale_agecondition <- scale_color_manual(name="Age\nGroup",
                                   values = myColors_agecondition, drop = T, limits=force)



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



myColors_celltype_fetalonly <- c("#a48cbf","#ba9797","#95a38b", 
                                 "#798694","#b09d9d","#c1a1c2",
                                 "#a3869a","#ad939b","#d4a5a6",
                                 "#917a8d")
color_possibilities_celltype_fetalonly<-c( "Fetal B cell","Fetal Erythrocytes","Fetal ILC precursor",
                                           "Fetal HSC/MPP","Fetal Early Erythrocytes","Fetal Mast cell",
                                           "Fetal Kupffer Cell","Fetal VCAM1+ Erythroblastic Island Macrophage", "Fetal Megakaryocyte",
                                           "Fetal Mono-NK"  )  

names(myColors_celltype_fetalonly) <- color_possibilities_celltype_fetalonly

#### Combined
combo_colors<-c(myColors_celltype, myColors_celltype_fetal,myColors_celltype_fetalonly)
combo_colors<-combo_colors[which(names(combo_colors)%in%c("B-cells","Mature B-cells","Plasma cells","pre B cell","pro B cell","B cell","pre pro B cell ","Fetal B cell",
                                                          
                                                          "CD3+ T-cells","gd T-cells","NK-like cells", "NK and T cells","Early lymphoid/T lymphocyte","CLNK T-cells",
                                                          "Cycling T-cells","NK","ILC precursor","NKT cells","Fetal ILC precursor",
                                                          
                                                          "Cholangiocytes", "LSEC","Endothelial cell",
                                                          
                                                          "Myeloid cells","RR Myeloid","Mono-Mac","Kupffer Cell","KC Like","Macrophage\n(MHCII high)","Macrophage\n(CLEC9A high)",
                                                          "Monocyte","pDC precursor","Megakaryocyte","Fetal Megakaryocyte",
                                                          "Cycling Myeloid","DC1","DC2","Monocyte-DC precursor","Neutrophil-myeloid progenitor","Mono-NK",
                                                          "Fetal Kupffer Cell","Fetal Mono-NK", 
                                                          
                                                          "Myeloid Erythrocytes\n(phagocytosis)","VCAM1+ Erythroblastic Island Macrophage",
                                                          "Fetal VCAM1+ Erythroblastic Island Macrophage",
                                                          "Early Erythrocytes","Mid  Erythroid","Late Erythroid","Erythrocytes","MEMP","Fetal Erythrocytes",
                                                          "Fetal Early Erythrocytes", 
                                                          "Neutrophil","Neutrophil\n(DEFA+)",
                                                          "Mast cell","Fetal Mast cell", 
                                                          "HSC","HSC/MPP","Fetal HSC/MPP",
                                                          "Hepatocytes","Hepatocyte","Fetal Hepatocytes",
                                                          "Doublet","Fetal NA"))]

combo_colors<-combo_colors[match(c("B-cells","Mature B-cells","Plasma cells","pre B cell","pro B cell","B cell","pre pro B cell ","Fetal B cell",
                                   
                                   "CD3+ T-cells","gd T-cells","NK-like cells", "NK and T cells","Early lymphoid/T lymphocyte","CLNK T-cells",
                                   "Cycling T-cells","NK","ILC precursor","NKT cells","Fetal ILC precursor",
                                   
                                   "Cholangiocytes", "LSEC","Endothelial cell",
                                   
                                   "Myeloid cells","RR Myeloid","Mono-Mac","Kupffer Cell","KC Like","Macrophage\n(MHCII high)","Macrophage\n(CLEC9A high)",
                                   "Monocyte","pDC precursor","Megakaryocyte","Fetal Megakaryocyte",
                                   "Cycling Myeloid","DC1","DC2","Monocyte-DC precursor","Neutrophil-myeloid progenitor","Mono-NK",
                                   "Fetal Kupffer Cell","Fetal Mono-NK", 
                                   
                                   "Myeloid Erythrocytes\n(phagocytosis)","VCAM1+ Erythroblastic Island Macrophage",
                                   "Fetal VCAM1+ Erythroblastic Island Macrophage",
                                   "Early Erythrocytes","Mid  Erythroid","Late Erythroid","Erythrocytes","MEMP","Fetal Erythrocytes",
                                   "Fetal Early Erythrocytes", 
                                   "Neutrophil","Neutrophil\n(DEFA+)",
                                   "Mast cell","Fetal Mast cell", 
                                   "HSC","HSC/MPP","Fetal HSC/MPP",
                                   "Hepatocytes","Hepatocyte","Fetal Hepatocytes",
                                   "Doublet","Fetal NA"),names(combo_colors))]

