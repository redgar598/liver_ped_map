d10x<-d10x.combined
gene<-c("CCR1","CCL3")
interaction_pair<-("CCL3_CCR1")


################
## 2 gene UMAP
################
scale_values <- function(x){round(((x-min(x))/(max(x)-min(x)))*100)}

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
    geom_point(aes(UMAP_1,UMAP_2),plt_myeloid, color=plt_myeloid$mix,size=0.9)+xlab("UMAP 1")+ylab("UMAP 2")+
    #geom_rect(data=plt_median, mapping=aes(xmin=min(mean_umap1)*1.1, xmax=max(mean_umap1)*1.21, ymin=min(mean_umap2)*1.25, ymax=max(mean_umap2)*1.5), fill = "white", alpha=0.05)+
    geom_point(aes(mean_umap1,mean_umap2, fill=CellType_refined), data=plt_median,size=4, shape=21, color="white")+
    geom_curve(
      data = differentialexp_plt_notself_pair,
      aes(x = Cell1x, y = Cell1y, xend = Cell2x, yend = Cell2y, alpha=as.factor(differential)), 
      arrow = arrow(length = unit(0.01, "npc"),type = "closed"),
      color = "white",curvature = -0.3,
      lineend = "round", size=1.5) + 
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
