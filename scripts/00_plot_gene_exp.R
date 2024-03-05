
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
    geom_point(color="black",size=1)+xlab("UMAP 1")+ylab("UMAP 2")+
    geom_point(aes(color=log(gene_exp_limited)),size=0.5)+xlab("UMAP 1")+ylab("UMAP 2")+
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




plot_gene_violin<-function(d10x, gene, log_exp=F){
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
  
  if(log_exp==T){
  ggplot(plt_myeloid, aes(age_condition,log(expression)))+
    geom_violin(fill="grey80",color="white")+geom_boxplot(aes(fill=age_condition),width=0.1)+fillscale_agecondition+
    theme_bw()+th_present+xlab("Age Group")+ylab(paste(gene, "Expression (log)"))+
    theme(legend.position = "none")}else{
  ggplot(plt_myeloid, aes(age_condition,expression))+
    geom_violin(fill="grey80",color="white")+geom_boxplot(aes(fill=age_condition),width=0.1)+fillscale_agecondition+
    theme_bw()+th_present+xlab("Age Group")+ylab(paste(gene, "Expression (log)"))+
    theme(legend.position = "none")}
    
}


plot_gene_violin_fetal<-function(d10x, gene){
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
  
  plt_myeloid$age_condition<-factor(plt_myeloid$age_condition, levels=c("Fetal Healthy","Ped Healthy","Ped IFALD","Adult Healthy"))
  
  ggplot(plt_myeloid, aes(age_condition,log(expression)))+
    geom_violin(fill="grey80",color="white")+geom_boxplot(aes(fill=age_condition),width=0.1)+fillscale_agecondition_fetal+
    theme_bw()+th_present+xlab("Age Group")+ylab(paste(gene, "Expression (log)"))+
    theme(legend.position = "none")
  
}

plot_heat_map<-function(d10x, gene, cellsubset, sig_label){ 
  DefaultAssay(d10x) <- "RNA"
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  
  if(is.character(cellsubset)){meta_myeloid<-meta_myeloid[which(meta_myeloid$CellType_refined%in%cellsubset),]}
  
  cell_num_all<-as.data.frame(table(meta_myeloid$age_condition, meta_myeloid$CellType_refined))
  colnames(cell_num_all)<-c("age_condition","CellType_refined","CellCount")
  
  gene_exp<-FetchData(d10x, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  
  meta_myeloid<-meta_myeloid[,c("age_condition","CellType_refined","cell")]
  plt_myeloid<-merge(meta_myeloid, gene_exp, by='cell')
  
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
  
  plt_exp_summary$age_condition<-as.factor(plt_exp_summary$age_condition)
  levels(plt_exp_summary$age_condition)<-c("Ped\nHealthy","Ped\nIFALD","Adult\nHealthy")
  
  plt_exp_summary$CellType_refined<-factor(plt_exp_scaled$CellType_refined, levels=cellsubset)
  
  if(sig_label==T){
    sig<-de[which(de$variable%in%gene),]
    ggplot()+
      geom_tile(aes( age_condition,variable, fill=scaled), plt_exp_summary)+facet_grid(.~CellType_refined, scales = "free_y", space = "free_y")+
      th+theme_classic()+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, 'RdBu')), space='Lab')(100), name="Scaled\nMean\nExpression")+
      ylab("")+xlab("")+geom_text(data=sig, aes(age_condition,variable, label=label))
  }else{  
    ggplot(plt_exp_summary, aes( age_condition,variable, fill=scaled))+
      geom_tile()+facet_grid(.~CellType_refined, scales = "free_y", space = "free_y")+
      th+theme_classic()+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, 'RdBu')), space='Lab')(100), name="Scaled\nMean\nExpression")+
      ylab("")+xlab("")
  }
}


plot_heat_map_fetal<-function(d10x, gene, cellsubset){ 
  DefaultAssay(d10x) <- "RNA"
  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
  
  if(is.character(cellsubset)){plt_myeloid<-plt_myeloid[which(plt_myeloid$CellType_harmonized%in%cellsubset),]}
  
  cell_num_all<-as.data.frame(table(plt_myeloid$age_condition, plt_myeloid$CellType_harmonized))
  colnames(cell_num_all)<-c("age_condition","CellType_harmonized","CellCount")
  
  gene_exp<-FetchData(d10x, vars=gene)
  gene_exp$cell<-rownames(gene_exp)
  
  plt_myeloid<-plt_myeloid[,c("age_condition","CellType_harmonized","cell")]
  plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')
  
  melt_exp<-melt(plt_myeloid)
  
  
  ### scale each genen across cells then take mean for gene in each cell type
  scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}
  # 
  # plt_exp_scaled <- melt_exp %>% group_by(variable) %>%
  #   dplyr::mutate(scaled = scale_this(value))
  
  plt_exp_summary <- melt_exp %>% 
    group_by(variable, CellType_harmonized, age_condition) %>%
    dplyr::summarize(Mean = mean(value, na.rm=TRUE))
  
  plt_exp_scaled <- plt_exp_summary %>% group_by(variable) %>%
    dplyr::mutate(scaled = scale_this(Mean))
  plt_exp_summary<-as.data.frame(plt_exp_scaled)
  
  plt_exp_summary$variable<-factor(plt_exp_scaled$variable, levels=rev(gene))
  
  plt_exp_summary$age_condition<-factor(plt_exp_summary$age_condition, levels=c("Fetal Healthy","Ped Healthy","Ped IFALD","Adult Healthy"))
  
  ggplot(plt_exp_summary, aes( age_condition,variable, fill=scaled))+
    geom_tile()+facet_grid(.~CellType_harmonized, scales = "free_y", space = "free_y")+
    th+theme_classic()+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, 'RdBu')), space='Lab')(100), name="Scaled\nMean\nExpression")+
    ylab("")+xlab("")
}

## expression PCA highlight cell type
plot_gene_PCA<-function(d10x, gene, percentile, split, highlight){
  ### plot individual genes
  if(length(grep("PC_1", colnames(d10x)))>0){
    plt_myeloid<-d10x@meta.data
  }else{
    DefaultAssay(d10x) <- "RNA"
    umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "pca"))#
    umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
    meta_myeloid<-d10x@meta.data
    meta_myeloid$cell<-rownames(meta_myeloid)
    plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")}
  
  cell_num_all<-as.data.frame(table(plt_myeloid$age_condition))
  colnames(cell_num_all)<-c("age_condition","CellCount")
  
  len_x_bar<-((range(plt_myeloid$PC_1))[2]-(range(plt_myeloid$PC_1))[1])/10
  len_y_bar<-((range(plt_myeloid$PC_2))[2]-(range(plt_myeloid$PC_2))[1])/10
  arr <- list(x = min(plt_myeloid$PC_1), y = min(plt_myeloid$PC_2), x_len = len_x_bar, y_len = len_y_bar)
  
  
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
  
  
  if(is.na(highlight)){
    exp_umap<-ggplot(plt_myeloid, aes(PC_1,PC_2))+
      geom_point(aes(color=log(gene_exp_limited)),size=0.75)+xlab("PC1")+ylab("PC2")+
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
  }else{
    plt_myeloid$highlight<-"0"
    plt_myeloid$highlight[which(plt_myeloid$CellType_refined==highlight)]<-"1"
    
    exp_umap<-ggplot(plt_myeloid, aes(PC_1,PC_2))+
      geom_point(aes(color=log(gene_exp_limited)),size=0.75)+
      geom_point(data=plt_myeloid[which(plt_myeloid$highlight==1),], size = 0.9, colour= "black", stroke = 1)+
      geom_point(aes(color=log(gene_exp_limited)),data=plt_myeloid[which(plt_myeloid$highlight==1),], size=0.75)+
      xlab("PC1")+ylab("PC2")+
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
  
  
  if(split==T){ exp_umap + facet_wrap(~age_condition, ncol=3)}else{exp_umap}
}







################
## 2 gene UMAP
################
scale_values <- function(x){round(((x-min(x))/(max(x)-min(x)))*100)}


plot_gene_UMAP_2gene<-function(d10x, gene){
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
  
  gene_exp$gene1_scaled<-scale_values(gene_exp[,1])
  gene_exp$gene2_scaled<-scale_values(gene_exp[,2])
  gene_exp$gene1_round<-round(gene_exp[,1])
  gene_exp$gene2_round<-round(gene_exp[,2])
  
  gene1_levels<-unique(gene_exp$gene1_scaled)[order(unique(gene_exp$gene1_scaled))]
  gene2_levels<-unique(gene_exp$gene2_scaled)[order(unique(gene_exp$gene2_scaled))]
  
  dat <- expand.grid(gene1_scaled=gene1_levels, gene2_scaled=gene2_levels)
  dat <- within(dat, mix <- rgb(green=gene1_scaled, red=gene2_scaled, blue=40, maxColorValue=100)) #maybe maxColorValue=150
  
  gene_exp<-merge(gene_exp, dat, by=c("gene1_scaled", "gene2_scaled"))
  
  plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')
  plt_myeloid$point_order<-sapply(1:nrow(plt_myeloid), function(x) sum(plt_myeloid$gene1_scaled[x], plt_myeloid$gene2_scaled[x]))
  plt_myeloid<-plt_myeloid[order(plt_myeloid$point_order),]
  plt_myeloid$mix[which(plt_myeloid$point_order==0)]<-"grey"
  
  UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(size = 0.06, colour= "black", stroke = 1)+
    geom_point(color=plt_myeloid$mix,size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
    annotate("segment", 
             x = arr$x, xend = arr$x + c(arr$x_len, 0), 
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(4, 'pt'))) +
    theme_void()+theme(legend.text=element_text(size=6),
                       legend.title=element_text(size=8), 
                       plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=5,hjust = 0.05),
                       axis.title.y = element_text(size=5,hjust = 0.05,angle = 90))
  
  dat_mini<-dat[which(dat$gene1_scaled %in% seq(0,100, by=10) & dat$gene2_scaled %in% seq(0,100, by=10)),]
  
  convert_scale_gene1<-plt_myeloid[!duplicated(plt_myeloid[,c("gene1_scaled",  "gene1_round")]),c("gene1_scaled", "gene1_round")]
  convert_scale_gene2<-plt_myeloid[!duplicated(plt_myeloid[,c("gene2_scaled", "gene2_round")]),c("gene2_scaled", "gene2_round")]
  
  dat_mini<-merge(dat_mini, convert_scale_gene1, by="gene1_scaled")
  dat_mini<-merge(dat_mini, convert_scale_gene2, by="gene2_scaled")
  dat_mini$mix[1]<-"grey"
  
  dat_mini$label1<-as.factor(dat_mini$gene1_round+(dat_mini$gene1_scaled/100))
  dat_mini$label2<-as.factor(dat_mini$gene2_round+(dat_mini$gene2_scaled/100))
  
  legend<-ggplot(dat_mini, aes(x=label1, y=label2)) + 
    geom_tile(aes(fill=mix), color="white") + theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    scale_fill_identity()+
    xlab(paste(gene[1], "Expression"))+
    ylab(paste(gene[2], "Expression"))+theme_minimal()+
    scale_x_discrete(breaks = c(0,  round(max(dat_mini$gene1_round+(dat_mini$gene1_scaled/100)))))+
    scale_y_discrete(breaks = c(0,  round(max(dat_mini$gene2_round+(dat_mini$gene2_scaled/100)))))
  
  
  plot_grid(UMAP,
            plot_grid(NULL,legend,NULL,ncol=1,rel_heights = c(1,2,1)),
            rel_widths = c(2,1))
}




## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}




plot_gene_UMAP_2gene_network<-function(d10x, gene,interaction_pair){
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
  
  gene_exp$gene1_scaled<-scale_values(gene_exp[,1])
  gene_exp$gene2_scaled<-scale_values(gene_exp[,2])
  gene_exp$gene1_round<-round(gene_exp[,1])
  gene_exp$gene2_round<-round(gene_exp[,2])
  
  gene1_levels<-unique(gene_exp$gene1_scaled)[order(unique(gene_exp$gene1_scaled))]
  gene2_levels<-unique(gene_exp$gene2_scaled)[order(unique(gene_exp$gene2_scaled))]
  
  dat <- expand.grid(gene1_scaled=gene1_levels, gene2_scaled=gene2_levels)
  dat <- within(dat, mix <- rgb(green=gene1_scaled, red=gene2_scaled, blue=40, maxColorValue=100)) #maybe maxColorValue=150
  
  gene_exp<-merge(gene_exp, dat, by=c("gene1_scaled", "gene2_scaled"))
  
  plt_myeloid<-merge(plt_myeloid, gene_exp, by='cell')
  plt_myeloid$point_order<-sapply(1:nrow(plt_myeloid), function(x) sum(plt_myeloid$gene1_scaled[x], plt_myeloid$gene2_scaled[x]))
  plt_myeloid<-plt_myeloid[order(plt_myeloid$point_order),]
  plt_myeloid$mix[which(plt_myeloid$point_order==0)]<-"grey"
  
  differentialexp_plt_notself_pair<-differentialexp_plt_notself[which(differentialexp_plt_notself$interacting_pair==interaction_pair),]
  
  
  ## Network UMAP overlay
  umap_network<-ggplot()+
    geom_point(aes(UMAP_1,UMAP_2),plt_myeloid, size = 1, colour= "black", stroke = 1)+
    geom_point(aes(UMAP_1,UMAP_2),plt_myeloid, color=plt_myeloid$mix,size=0.75)+xlab("UMAP 1")+ylab("UMAP 2")+
    #geom_rect(data=plt_median, mapping=aes(xmin=min(mean_umap1)*1.1, xmax=max(mean_umap1)*1.21, ymin=min(mean_umap2)*1.25, ymax=max(mean_umap2)*1.5), fill = "white", alpha=0.05)+
    geom_point(aes(mean_umap1,mean_umap2, fill=CellType_refined), data=plt_median,size=4, shape=21, color="white")+
    geom_curve(
      data = differentialexp_plt_notself_pair,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y, alpha=as.factor(differential)), 
      arrow = arrow(length = unit(0.01, "npc"),type = "closed"),
      color = "white",curvature = -0.3,
      lineend = "round", size=1.15) + 
    geom_curve(
      data = differentialexp_plt_notself_pair,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y, alpha=as.factor(differential)), 
      arrow = arrow(length = unit(0.01, "npc"),type = "closed"),
      color = "grey10",curvature = -0.3,
      lineend = "round") +  scale_alpha_manual(values=c(1,0.5)) +
    fillscale_cellType+
    annotate("segment", 
             x = arr$x, xend = arr$x + c(arr$x_len, 0), 
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(4, 'pt'))) +
    theme_void()+theme(legend.text=element_text(size=6),
                       legend.title=element_text(size=8), 
                       plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=7,hjust = 0.05, vjust = 10),
                       axis.title.y = element_text(size=7,hjust = 0.05, vjust=-9,angle = 90),
                       legend.position = "none")
  
  
  ## Cell type node legend
  nice_legend<-get_leg(ggplot()+ 
                         geom_point(aes(mean_umap1,mean_umap2, fill=CellType_refined), 
                                    data=plt_median,size=3, shape=21, color="white")+
                         fillscale_cellType+theme_void()+
                         theme(legend.text=element_text(size=6),legend.title=element_text(size=8)))
  
  ## Expression level grid legend
  dat_mini<-dat[which(dat$gene1_scaled %in% seq(0,100, by=10) & dat$gene2_scaled %in% seq(0,100, by=10)),]
  
  convert_scale_gene1<-plt_myeloid[!duplicated(plt_myeloid[,c("gene1_scaled",  "gene1_round")]),c("gene1_scaled", "gene1_round")]
  convert_scale_gene2<-plt_myeloid[!duplicated(plt_myeloid[,c("gene2_scaled", "gene2_round")]),c("gene2_scaled", "gene2_round")]
  
  dat_mini<-merge(dat_mini, convert_scale_gene1, by="gene1_scaled")
  dat_mini<-merge(dat_mini, convert_scale_gene2, by="gene2_scaled")
  dat_mini$mix[1]<-"grey"
  
  dat_mini$label1<-as.factor(dat_mini$gene1_round+(dat_mini$gene1_scaled/100))
  dat_mini$label2<-as.factor(dat_mini$gene2_round+(dat_mini$gene2_scaled/100))
  
  exp_legend<-ggplot(dat_mini, aes(x=label1, y=label2)) + 
    geom_tile(aes(fill=mix), color="white") + 
    scale_fill_identity()+
    xlab(paste(gene[1], "Expression"))+
    ylab(paste(gene[2], "Expression"))+theme_minimal()+
    scale_x_discrete(breaks = c(0,  round(max(dat_mini$gene1_round+(dat_mini$gene1_scaled/100)))))+
    scale_y_discrete(breaks = c(0,  round(max(dat_mini$gene2_round+(dat_mini$gene2_scaled/100)))))+
    theme(axis.title = element_text(size=9),
          axis.text = element_text(size=7))

  
  
  ## Interaction Legend
  ligand<-strsplit(interaction_pair,"_")[[1]][1]
  receptor<-strsplit(interaction_pair,"_")[[1]][2]
  
  if(sum(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==2){
    interaction_legend<- ggplot()+
      annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
      annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
      annotate("text", label=ligand, x = 1, y = 1, color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("text", label=receptor, x = 2, y = 1,color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("segment", 
               x = 1.3, xend = 1.7 , 
               y = 1, yend = 1, size=0.25,color="black",
               arrow = arrow(type = "closed", length = unit(2, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
    
  }else{if(which(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==1){
    
    interaction_legend<-ggplot()+
      annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
      annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
      annotate("text", label=ligand, x = 1, y = 1, color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("text", label=receptor, x = 2, y = 1,color="black")+
      annotate("segment", 
               x = 1.3, xend = 1.7 , 
               y = 1, yend = 1, size=0.6,color="black",
               arrow = arrow(type = "closed", length = unit(5, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
  } else {
    
    if(sum(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==2){
      interaction_legend<- ggplot()+
        annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
        annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
        annotate("text", label=ligand, x = 1, y = 1, color="black")+
        annotate("text", label=receptor, x = 2, y = 1,color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
        annotate("segment", 
                 x = 1.3, xend = 1.7 , 
                 y = 1, yend = 1, size=0.25,color="black",
                 arrow = arrow(type = "closed", length = unit(2, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
    }}}
  
  
  # plot_grid(umap_network,
  #           plot_grid(interaction_legend,nice_legend, exp_legend,ncol=1,rel_heights = c(0.4,1.5,1.6)),
  #           rel_widths = c(2,1))
  
  plot_grid(umap_network,
            plot_grid(interaction_legend,nice_legend, plot_grid(NULL,exp_legend,NULL, ncol=3, rel_widths = c(0.3,1,0.3)) ,ncol=1,rel_heights = c(0.25,1,0.6)),
            rel_widths = c(2,1))
  
}





plot_gene_UMAP_2gene_network_notblend<-function(d10x, gene,interaction_pair, percentile){
  differentialexp_plt_notself_pair<-differentialexp_plt_notself[which(differentialexp_plt_notself$interacting_pair==interaction_pair),]
  
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
  
  
  ### non zero expression over the 90th percentile
  exp_limit1<-quantile(plt_myeloid[, which(colnames(plt_myeloid)==gene[1])], percentile)
  exp_limit2<-quantile(plt_myeloid[, which(colnames(plt_myeloid)==gene[2])], percentile)
  
  both_gene_over<-plt_myeloid[which(plt_myeloid[, which(colnames(plt_myeloid)==gene[1])]>=exp_limit1 &
                                      plt_myeloid[, which(colnames(plt_myeloid)==gene[1])]>0 &
                                      plt_myeloid[, which(colnames(plt_myeloid)==gene[2])]>=exp_limit2 &
                                      plt_myeloid[, which(colnames(plt_myeloid)==gene[2])]>0),]
  both_gene_over$color<-"Both Highly Expressed"
  
  gene1_over<-plt_myeloid[which(plt_myeloid[, which(colnames(plt_myeloid)==gene[1])]>=exp_limit1 &
                                  plt_myeloid[, which(colnames(plt_myeloid)==gene[1])]>0),]
  gene1_over<-gene1_over[which(!(gene1_over$cell%in%both_gene_over$cell)),]
  gene1_over$color<-paste(gene[1],"Highly Expressed")
  
  gene2_over<-plt_myeloid[which(plt_myeloid[, which(colnames(plt_myeloid)==gene[2])]>=exp_limit2 &
                                  plt_myeloid[, which(colnames(plt_myeloid)==gene[2])]>0),]
  gene2_over<-gene2_over[which(!(gene2_over$cell%in%both_gene_over$cell)),]
  gene2_over$color<-paste(gene[2],"Highly Expressed")
  
  plt_point_color<-rbind(both_gene_over, gene1_over, gene2_over)
  plt_point_color<-plt_point_color[order(plt_point_color$cell),]
  
  
  ## Network UMAP overlay
  umap_network<-ggplot()+
    geom_point(aes(UMAP_1,UMAP_2),plt_myeloid, size = 1, colour= "black", stroke = 1)+
    geom_point(aes(UMAP_1,UMAP_2),plt_myeloid, color="grey",size=0.75)+
    geom_point(aes(UMAP_1,UMAP_2, color=color),plt_point_color, size=0.5)+
    geom_rect(data=plt_median, mapping=aes(xmin=min(mean_umap1)*1.1, xmax=max(mean_umap1)*1.21, ymin=min(mean_umap2)*1.25, ymax=max(mean_umap2)*1.5), fill = "white", alpha=0.05)+
    xlab("UMAP 1")+ylab("UMAP 2")+
    scale_color_manual(values=c("#88a000","#b80783","#03008e"))+
    geom_point(aes(mean_umap1,mean_umap2, fill=CellType_refined), data=plt_median,size=4, shape=21, color="white")+
    geom_curve(
      data = differentialexp_plt_notself_pair,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y, alpha=as.factor(differential)), 
      arrow = arrow(length = unit(0.01, "npc"),type = "closed"),
      color = "white",curvature = -0.3,
      lineend = "round", size=1.15) + 
    geom_curve(
      data = differentialexp_plt_notself_pair,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y, alpha=as.factor(differential)), 
      arrow = arrow(length = unit(0.01, "npc"),type = "closed"),
      color = "grey10",curvature = -0.3,
      lineend = "round") +  scale_alpha_manual(values=c(1,0.5)) +
    fillscale_cellType+
    annotate("segment", 
             x = arr$x, xend = arr$x + c(arr$x_len, 0), 
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(4, 'pt'))) +
    theme_void()+theme(legend.text=element_text(size=6),
                       legend.title=element_text(size=8), 
                       plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=7,hjust = 0.05, vjust = 10),
                       axis.title.y = element_text(size=7,hjust = 0.05, vjust=-9,angle = 90),
                       legend.position = "none")
  
  
  ## Cell type node legend
  nice_legend<-get_leg(ggplot()+ 
                         geom_point(aes(mean_umap1,mean_umap2, fill=CellType_refined), 
                                    data=plt_median,size=3, shape=21, color="white")+
                         fillscale_cellType+theme_void()+
                         theme(legend.text=element_text(size=8),legend.title=element_text(size=10),
                               plot.margin = margin(0.1,0.1,0.1,0.1, "cm")))
  
  ## gene expression legend
  exp_legend<-get_leg(ggplot()+
                        geom_point(aes(UMAP_1,UMAP_2, color=color),plt_point_color, size=2, alpha=0.5)+
                        scale_color_manual(values=c("#88a000","#b80783","#03008e"), 
                                           name=paste("Gene Expression\n(Higher than ", percentile*100, "th\npercentile for gene",sep=""))+
                        theme(legend.text=element_text(size=8),legend.title=element_text(size=10),
                              plot.margin = margin(0.1,0.1,0.1,0.1, "cm")))
  
  
  ## Interaction Legend
  ligand<-strsplit(interaction_pair,"_")[[1]][1]
  receptor<-strsplit(interaction_pair,"_")[[1]][2]
  
  if(sum(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==2){
    interaction_legend<- ggplot()+
      annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
      annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
      annotate("text", label=ligand, x = 1, y = 1, color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("text", label=receptor, x = 2, y = 1,color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("segment", 
               x = 1.3, xend = 1.7 , 
               y = 1, yend = 1, size=0.25,color="black",
               arrow = arrow(type = "closed", length = unit(2, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
    
  }else{if(which(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==1){
    
    interaction_legend<-ggplot()+
      annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
      annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
      annotate("text", label=ligand, x = 1, y = 1, color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
      annotate("text", label=receptor, x = 2, y = 1,color="black")+
      annotate("segment", 
               x = 1.3, xend = 1.7 , 
               y = 1, yend = 1, size=0.6,color="black",
               arrow = arrow(type = "closed", length = unit(5, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
  } else {
    
    if(sum(strsplit(interaction_pair,"_")[[1]]%in%de_KC$X)==2){
      interaction_legend<- ggplot()+
        annotate("text", label="Ligand", x = 1, y = 1.01, color="black",  fontface =2)+
        annotate("text", label="Receptor", x = 2, y = 1.01,color="black", fontface =2)+
        annotate("text", label=ligand, x = 1, y = 1, color="black")+
        annotate("text", label=receptor, x = 2, y = 1,color=myColors_celltype[which(names(myColors_celltype)=="KC Like")])+
        annotate("segment", 
                 x = 1.3, xend = 1.7 , 
                 y = 1, yend = 1, size=0.25,color="black",
                 arrow = arrow(type = "closed", length = unit(2, 'pt')))+theme_void()+ylim(0.99,1.02)+xlim(0.3,3)
    }}}
  
  
  plot_grid(umap_network,
            plot_grid(interaction_legend,nice_legend, exp_legend,ncol=1,rel_heights = c(0.2,1,0.5)),
            rel_widths = c(2,1))
  
}


