#library(scales)
save_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h)}

save_fetal_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/fetal/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/fetal/jpeg/"),name,".jpeg", sep=""), w=w, h=h)}




myColors_celltype <- c("#660cc7","#8642d1","#5612a3","#4b911d","#7a4202","#006d2c","#6bbce8","#3469ad","#d17906","#b01e87","#60ba5d","#207537",
                       "#a0c487","#d9a5a5","#87a4c9","#e8c392","#dea4ce","#79639a",
                       "#207537","#fa61ad","#b80783","#994676","#431039","#cb181d","maroon1","#67000d","grey","#ce1256","#a6d96a")
color_possibilities_celltype<-c("B-cells","Mature B-cells","Plasma cells","CD3+ T-cells","Cholangiocytes","gd T-cells","Hepatocytes","HSC","LSEC","Myeloid cells","NK-like cells", "NK and T cells",
                                "NKT cells\n(Hepatocyte Like)","Cholangiocytes\n(Hepatocyte Like)","HSC\n(Hepatocyte Like)","LSEC\n(Hepatocyte Like)","Myeloid cells\n(Hepatocyte Like)","B-cells\n(Hepatocyte Like)",
                                "NKT cells","RR Myeloid","KC Like","Neutrophil","Neutrophil\n(DEFA+)","Erythrocytes","Mast cell","Myeloid Erythrocytes\n(phagocytosis)",
                                "Doublet","Macrophage\n(MHCII high)","Cycling T-cells")
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




