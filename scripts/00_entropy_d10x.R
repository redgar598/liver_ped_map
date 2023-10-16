
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



## entropy on a seurat object
entropy_d10<-function(d10x, covariate){
  
  ## 95% from one age group
  level_num<-length(unique(d10x@meta.data[,covariate]))
  print(paste("Entrophy theshold if 95% of samples in a cluster from 1 covariate level (and the other 5% a random mix of the",level_num-1,"other factor levels):",
              round(entropy(c(sample(as.character(unique(d10x@meta.data[,covariate])[1:(level_num-1)]), 100, replace=T), 
                              rep(as.character(unique(d10x@meta.data[,covariate])[level_num]),1900))),2)))
  
  entrophy_cluster_df<-do.call(rbind, lapply(as.character(unique(d10x$seurat_clusters)), function(cluster){
    data.frame(seurat_clusters=cluster, entropy=entropy(d10x@meta.data[,covariate][which(d10x$seurat_clusters==cluster)]))
  }))
  

  individual_UMAP<-DimPlot(d10x, group.by = covariate) + theme(legend.position = "none") + ggtitle("")
  
  cell_cluster_count<-d10x@meta.data %>%  group_by(seurat_clusters,rlang::parse_expr(covariate) %>% eval) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  cell_cluster_count<-as.data.frame(cell_cluster_count)
  cell_cluster_count<-merge(cell_cluster_count, entrophy_cluster_df, by="seurat_clusters")
  colnames(cell_cluster_count)[2]<-covariate
  cell_cluster_count}


##############
### plots for entropy
##############

entropy_plt<-function(plt_entropy, covariate, d10x){
  max_count<-as.data.frame(plt_entropy %>% 
                             group_by(seurat_clusters) %>% 
                             summarise(total = sum(n)))
  
  entrophy_label<-plt_entropy[!duplicated(plt_entropy[,c("seurat_clusters",  "entropy")]), c("seurat_clusters",  "entropy")]
  entrophy_label<-merge(entrophy_label, max_count, by="seurat_clusters")
  
  bar_individual<-ggplot() + 
    geom_bar(aes(fill=rlang::parse_expr(covariate) %>% eval, y=n, x=seurat_clusters),plt_entropy, position="stack", stat="identity", color="black")+
    theme_bw()+th+ylab("Cell Count")+xlab("Seurat Cluster")+
    geom_text(aes(label=round(entropy,2), y=(total+(0.1*max(total))), x=seurat_clusters), entrophy_label)+
    scale_fill_manual(values=sequential_hcl(length(unique(d10x@meta.data[,covariate])), "Batlow"), name=covariate)

  individal_UMAP<-DimPlot(d10x, group.by = covariate,pt.size=1) + theme(legend.position = "none") +ggtitle("")+
    scale_color_manual(values=sequential_hcl(length(unique(d10x@meta.data[,covariate])), "Batlow"))
  cluster_UMAP<-DimPlot(d10x, pt.size=1, label=T) + theme(legend.position = "none") +ggtitle("")
  
  plot_grid(plot_grid(cluster_UMAP, individal_UMAP), bar_individual, ncol=1, rel_heights = c(1.5,1))}
