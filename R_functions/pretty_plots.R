#library(scales)
save_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h)}


myColors_diagnosis <- c("lightgrey","darkgoldenrod1","dodgerblue3","lightskyblue","khaki1","#819377","#993399","darkgrey")
color_possibilities_diagnosis<-c( "Control","UC","CD","IBD-U (CD-like)","IBD-U (UC-like)", "IBD-U", "Other.GI","Neonatal")
names(myColors_diagnosis) <- color_possibilities_diagnosis
fillscale_diagnosis <- scale_fill_manual(name="Diagnosis",
                               values = myColors_diagnosis, drop = T)
colscale_diagnosis <- scale_color_manual(name="Diagnosis",
                                         values = myColors_diagnosis, drop = T)


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

