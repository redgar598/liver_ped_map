
plot_gene_UMAP<-function(d10x, gene, percentile){
  ### plot individual genes
  
  DefaultAssay(d10x) <- "RNA"
  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
  
  cell_num_all<-as.data.frame(table(plt_myeloid$age_condition))
  colnames(cell_num_all)<-c("age_condition","CellCount")
  
  len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
  len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
  arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)
  
  
  gene_exp<-FetchData(d10x, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')
  
  exp_limit<-quantile(plt_myeloid[, which(colnames(plt_myeloid)==gene)], percentile)
  plt_myeloid$gene_exp_limited<-NA
  over_limit<-which(plt_myeloid[, which(colnames(plt_myeloid)==gene)]>exp_limit)
  plt_myeloid$gene_exp_limited[over_limit]<-plt_myeloid[over_limit, which(colnames(plt_myeloid)==gene)]
  plt_myeloid<-plt_myeloid[rev(order(plt_myeloid$gene_exp_limited)),]
  
  plt_myeloid<-rbind(plt_myeloid[which(is.na(plt_myeloid$gene_exp_limited)),],
                     plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),][(order(plt_myeloid[which(!(is.na(plt_myeloid$gene_exp_limited))),]$gene_exp_limited)),])
  
  ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(aes(color=log(gene_exp_limited)),size=0.75)+xlab("UMAP 1")+ylab("UMAP 2")+
    #facet_wrap(~age_condition,ncol=2)+
    #geom_text(aes(x = 7, y = -15, label=paste0("n = ",comma(CellCount))), cell_num_all)+
    scale_color_continuous_sequential(palette = "Blues 3", rev=T, 
                                      name=paste(gene, "\nExpression\n(log)"),na.value = "grey80")+
    annotate("segment", 
             x = arr$x, xend = arr$x + c(arr$x_len, 0), 
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(4, 'pt'))) +
    theme_void()+theme(legend.text=element_text(size=6),
                       legend.title=element_text(size=8), 
                       plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=5,hjust = 0.05),
                       axis.title.y = element_text(size=5,hjust = 0.05,angle = 90))
  }




plot_gene_violin<-function(d10x, gene){
  ### plot individual genes
  
  DefaultAssay(d10x) <- "RNA"
  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
  
  cell_num_all<-as.data.frame(table(plt_myeloid$age_condition))
  colnames(cell_num_all)<-c("age_condition","CellCount")
  
  gene_exp<-FetchData(d10x, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')
  
  colnames(plt_myeloid)[which(colnames(plt_myeloid)==gene)]<-"expression"
  
  ggplot(plt_myeloid, aes(age_condition,log(expression)))+
    geom_violin(fill="grey80",color="white")+geom_boxplot(aes(fill=age_condition),width=0.1)+fillscale_agecondition+
    theme_bw()+th_present+xlab("Age Group")+ylab(paste(gene, "Expression (log)"))+
    theme(legend.position = "none")
    
}




plot_heat_map<-function(d10x, gene, cellsubset){ 
  DefaultAssay(d10x) <- "RNA"
  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
  
  if(is.character(cellsubset)){plt_myeloid<-plt_myeloid[which(plt_myeloid$CellType_refined%in%cellsubset),]}
  
  cell_num_all<-as.data.frame(table(plt_myeloid$age_condition, plt_myeloid$CellType_refined))
  colnames(cell_num_all)<-c("age_condition","CellType_refined","CellCount")
  
  gene_exp<-FetchData(d10x, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  
  plt_myeloid<-plt_myeloid[,c("age_condition","CellType_refined","cell")]
  plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')
  
  melt_exp<-melt(plt_myeloid)
  
  
  ### scale each genen across cells then take mean for gene in each cell type
  scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}
  # 
  # plt_exp_scaled <- melt_exp %>% group_by(variable) %>%
  #   dplyr::mutate(scaled = scale_this(value))
  
  plt_exp_summary <- melt_exp %>% 
    group_by(variable, CellType_refined, age_condition) %>%
    dplyr::summarize(Mean = mean(value, na.rm=TRUE))
  
  plt_exp_scaled <- plt_exp_summary %>% group_by(variable) %>%
    dplyr::mutate(scaled = scale_this(Mean))
  plt_exp_summary<-as.data.frame(plt_exp_scaled)
  
  plt_exp_summary$variable<-factor(plt_exp_scaled$variable, levels=rev(gene))
  
  ggplot(plt_exp_summary, aes( age_condition,variable, fill=scaled))+
    geom_tile()+facet_grid(.~CellType_refined, scales = "free_y", space = "free_y")+
    th+theme_classic()+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, 'RdBu')), space='Lab')(100), name="Scaled\nMean\nExpression")+
    ylab("")+xlab("")
}
