#library(scales)
save_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h)}


myColors_celltype <- c("#4b911d","#b83b3b","#006d2c","#6bbce8","#3469ad","#d17906","#b01e87","#60ba5d","#a0c487","#d9a5a5","#87a4c9","#e8c392","#dea4ce")
color_possibilities_celltype<-c("CD3+ T-cells","Cholangiocytes","gd T-cells","Hepatocytes","HSC","LSEC","Myeloid cells","NK-like cells", "NKT cells\n(Hepatocyte Like)","Cholangiocytes\n(Hepatocyte Like)","HSC\n(Hepatocyte Like)","LSEC\n(Hepatocyte Like)","Myeloid cells\n(Hepatocyte Like)")
names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T)
colscale_cellType <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T)



myColors_age <- c("#348595","#d6604d")
color_possibilities_age<-c( "Adult","Ped")
names(myColors_age) <- color_possibilities_age
fillscale_age <- scale_fill_manual(name="Age\nGroup",
                                         values = myColors_age, drop = T)
colscale_age <- scale_color_manual(name="Age\nGroup",
                                         values = myColors_age, drop = T)






myColors_sampsite <- c("#1a9850","#a6d96a","cornflowerblue","cornflowerblue","#a6d96a","grey")
color_possibilities_sampsite<-c( "AC","SC","TI","proximal","distal","other")
names(myColors_sampsite) <- color_possibilities_sampsite
fillscale_sampsite <- scale_fill_manual(name="Sample Site",
                                         values = myColors_sampsite, drop = T)
colscale_sampsite <- scale_color_manual(name="Sample Site",
                                        values = myColors_sampsite, drop = T)





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



pass_col<-c("grey30","#9E0142", "#D53E4F", "#F46D43", "#FDAE61","#FEC776","#FEE08B", "#FFFFBF","#E6F598", "#C9E99E","#ABDDA4","#66C2A5","#4CA5B1","#3288BD", "#5E4FA2","#762a83","#3f007d")
names(pass_col)<-c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,21,23)



myColors_diagnosis_sample_samplesize <- c("lightgrey","grey80","darkgoldenrod1","darkgoldenrod1","dodgerblue3","#74abe1")
names(myColors_diagnosis_sample_samplesize) <- c( "Control full", "Control sub", "UC full", "UC sub","CD full", "CD sub")
fillscale_diagnosis_sample <- scale_fill_manual(name="Diagnosis",
                                                values = myColors_diagnosis_sample_samplesize, drop = T)
colscale_diagnosis_sample <- scale_color_manual(name="Diagnosis",
                                                values = myColors_diagnosis_sample_samplesize, drop = T)

