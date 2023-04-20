load(here("data","adult_ped_cellRefined_withDropletQC.rds"))

meta_min<-cell_label[,c("individual","Age")]
meta_min<-meta_min[!duplicated(meta_min),]

cell_label$index<-rownames(cell_label)

cell_counts<-cell_label %>% 
  group_by(individual, CellType_refined) %>% 
  summarise(count=length(unique(cell))) %>% 
  group_by(individual) %>%
  mutate(countT= sum(count)) %>%
  group_by(CellType_refined, add=TRUE) %>%
  mutate(per=100*count/countT)

cell_counts<-merge(cell_counts, meta_min, by="individual")

cell_counts$label<-sapply(1:nrow(cell_counts), function(x){
  if(length(grep("NPC", cell_counts$individual[x]))==1){
    paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1]," NPC)", sep="")
  }else{if(length(grep("TLH", cell_counts$individual[x])==1)){
    paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1]," TLH)", sep="")
  }else{
      paste(cell_counts$Age[x], "\n(", strsplit(cell_counts$individual[x],"_")[[1]][1],")", sep="")}}
})

cell_counts$label<-factor(cell_counts$label, c("2\n(C104)","11\n(C85)","12\n(C93)", "17\n(C64)", "17\n(C96)", 
                                               "26\n(C82)","48\n(C70)", "57\n(C97)","61\n(C68)",
                                               "65\n(C39 NPC)", "65\n(C39 TLH)", "67\n(C54)","69\n(C88)"))


all_composistion<-ggplot(cell_counts, aes(label, per, fill=CellType_refined))+geom_bar(stat = "identity", color="black")+
  theme_bw()+th+fillscale_cellType+xlab("")+ylab("Percent of Cells in Sample")
save_plts(all_composistion, "all_composistion", w=12,h=8)

#########
## Just myeloid
#########

cell_counts_myeloid<-cell_label[which(cell_label$CellType_refined%in%c("Macrophage\n(MHCII high)","KC Like","RR Myeloid" )),] %>% 
  group_by(individual, CellType_refined) %>% 
  summarise(count=length(unique(cell))) %>% 
  group_by(individual) %>%
  mutate(countT= sum(count)) %>%
  group_by(CellType_refined, add=TRUE) %>%
  mutate(per=100*count/countT)

cell_counts_myeloid<-merge(cell_counts_myeloid, meta_min, by="individual")

cell_counts_myeloid$label<-sapply(1:nrow(cell_counts_myeloid), function(x){
  if(length(grep("NPC", cell_counts_myeloid$individual[x]))==1){
    paste(cell_counts_myeloid$Age[x], "\n(", strsplit(cell_counts_myeloid$individual[x],"_")[[1]][1]," NPC)", sep="")
  }else{if(length(grep("TLH", cell_counts_myeloid$individual[x])==1)){
    paste(cell_counts_myeloid$Age[x], "\n(", strsplit(cell_counts_myeloid$individual[x],"_")[[1]][1]," TLH)", sep="")
  }else{
    paste(cell_counts_myeloid$Age[x], "\n(", strsplit(cell_counts_myeloid$individual[x],"_")[[1]][1],")", sep="")}}
})

cell_counts_myeloid$label<-factor(cell_counts_myeloid$label, c("2\n(C104)","11\n(C85)","12\n(C93)", "17\n(C64)", "17\n(C96)", 
                                               "26\n(C82)","48\n(C70)", "57\n(C97)","61\n(C68)",
                                               "65\n(C39 NPC)", "65\n(C39 TLH)", "67\n(C54)","69\n(C88)"))

myeloid_composistion<-ggplot(cell_counts_myeloid, aes(label, per, fill=CellType_refined))+geom_bar(stat = "identity", color="black")+
  theme_bw()+th+fillscale_cellType+xlab("")+ylab("Percent of Cells in Sample")
save_plts(myeloid_composistion, "myeloid_composistion", w=12,h=8)
