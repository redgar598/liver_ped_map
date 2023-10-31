
#compute Shannon entropy
entropy <- function(target) {
  freq <- table(target)/length(target)
  # vectorize
  vec <- as.data.frame(freq)[,2]
  #drop 0 to avoid NaN resulting from log2
  vec<-vec[vec>0]
  #compute entropy
  -sum(vec * log2(vec))
}

broad<-d10x@meta.data[,c("CellType_refined", "CellType_rough")]
broad<-broad[!duplicated(broad),]
broad$CellType_rough[grep("Cycling",broad$CellType_refined)]<-"Cycling"
broad$CellType_rough<-as.factor(broad$CellType_rough)
levels(broad$CellType_rough)<-c( "B-cells",  "Non-Immune", "Cycling","Erythrocytes" ,"Non-Immune","Non-Immune","Non-Immune","Myeloid cells","NK and T cells",  "Doublet")


## 95% from one age group
level_num<-length(unique(d10x@meta.data[,"age_condition"]))
print(paste("Entrophy theshold if 95% of samples in a cluster from 1 covariate level (and the other 5% a random mix of the",level_num-1,"other factor levels):",
            round(entropy(c(sample(as.character(unique(d10x@meta.data[,"age_condition"])[1:(level_num-1)]), 100, replace=T), 
                            rep(as.character(unique(d10x@meta.data[,"age_condition"])[level_num]),1900))),2)))

entrophy_cluster_df<-do.call(rbind, lapply(as.character(unique(d10x$CellType_refined)), function(cluster){
  data.frame(CellType_refined=cluster, entropy=entropy(d10x@meta.data[,"age_condition"][which(d10x$CellType_refined==cluster)]))
}))



cell_cluster_count<-d10x@meta.data %>%  group_by(CellType_refined,age_condition) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

cell_cluster_count<-as.data.frame(cell_cluster_count)
cell_cluster_count<-merge(cell_cluster_count, entrophy_cluster_df, by="CellType_refined")
cell_cluster_count


max_count<-as.data.frame(cell_cluster_count %>% 
                           group_by(CellType_refined) %>% 
                           summarise(total = sum(n)))

entrophy_label<-cell_cluster_count[!duplicated(cell_cluster_count[,c("CellType_refined",  "entropy")]), c("CellType_refined",  "entropy")]
entrophy_label<-merge(entrophy_label, max_count, by="CellType_refined")

entrophy_label<-merge(entrophy_label, broad,  by="CellType_refined")
cell_cluster_count<-merge(cell_cluster_count, broad,  by="CellType_refined")


cell_cluster_count$CellType_rough<-factor(cell_cluster_count$CellType_rough, levels=c( 
  "Non-Immune",
  "Myeloid cells","NK and T cells",  "B-cells",
  "Cycling","Erythrocytes" ,  "Doublet"))
entrophy_label$CellType_rough<-factor(entrophy_label$CellType_rough, levels=c( 
  "Non-Immune",
  "Myeloid cells","NK and T cells",  "B-cells",
  "Cycling","Erythrocytes" , "Doublet"))


bar_condition<-ggplot() + 
  geom_bar(aes(fill=age_condition, y=n, x=reorder(CellType_refined, -n)),cell_cluster_count, position="stack", stat="identity", color="black")+
  theme_bw()+th+ylab("Cell Count")+xlab("Seurat Cluster")+
  geom_text(aes(label=round(entropy,2), y=(total+(0.05*max(total))), x=CellType_refined), entrophy_label, size=3)+fillscale_agecondition+
  facet_grid(.~CellType_rough, scales="free_x", space = "free_x")+theme(strip.background =element_rect(fill="white"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
bar_condition
save_plts(bar_condition, "entrophy_celltype_age_condition", w=12,h=4)

