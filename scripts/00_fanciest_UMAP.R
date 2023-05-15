## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


fanciest_UMAP<-function(d10x, highlight, split){
  umap_mat_myeloid<-as.data.frame(Embeddings(object = d10x, reduction = "umap"))#
  umap_mat_myeloid$cell<-rownames(umap_mat_myeloid)
  meta_myeloid<-d10x@meta.data
  meta_myeloid$cell<-rownames(meta_myeloid)
  plt_myeloid<-merge(meta_myeloid, umap_mat_myeloid, by="cell")
  
  len_x_bar<-((range(plt_myeloid$UMAP_1))[2]-(range(plt_myeloid$UMAP_1))[1])/10
  len_y_bar<-((range(plt_myeloid$UMAP_2))[2]-(range(plt_myeloid$UMAP_2))[1])/10
  arr <- list(x = min(plt_myeloid$UMAP_1), y = min(plt_myeloid$UMAP_2), x_len = len_x_bar, y_len = len_y_bar)
  
  
  forlegned_plot<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(aes(fill=CellType_refined),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
    fillscale_cellType+theme_bw()+
    theme(legend.text = element_text(size=5),
          legend.title = element_text(size=6))
  nice_legend<-get_leg(forlegned_plot)
  #guides(colour = guide_legend(override.aes = list(size=0.5),byrow = TRUE))
  
  if(is.na(highlight)){
  fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
    geom_point(size = 0.06, colour= "black", stroke = 1)+
    geom_point(aes(color=CellType_refined),size=0.05)+xlab("UMAP 1")+ylab("UMAP 2")+
    colscale_cellType+
    annotate("segment", 
             x = arr$x, xend = arr$x + c(arr$x_len, 0), 
             y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
             arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
    theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                       axis.title.x = element_text(size=5,hjust = 0.05),
                       axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                       legend.position = "none")
  }else{
    plt_myeloid$highlight<-"0"
    plt_myeloid$highlight[which(plt_myeloid$CellType_refined==highlight)]<-"1"
    
    fanciest_UMAP<-ggplot(plt_myeloid, aes(UMAP_1,UMAP_2))+
      geom_point(aes(color=CellType_refined),size=0.05)+
      geom_point(data=plt_myeloid[which(plt_myeloid$highlight==1),], size = 0.06, colour= "black", stroke = 1)+
      geom_point(aes(color=CellType_refined),data=plt_myeloid[which(plt_myeloid$highlight==1),], size=0.05)+
      xlab("UMAP 1")+ylab("UMAP 2")+
      colscale_cellType+
      annotate("segment", 
               x = arr$x, xend = arr$x + c(arr$x_len, 0), 
               y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
               arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
      theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                         axis.title.x = element_text(size=5,hjust = 0.05),
                         axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                         legend.position = "none")
    }
  
  ## cell count
  
  if(split==T & !is.na(highlight)){  
    cell_num_all<-as.data.frame(table(plt_myeloid[which(plt_myeloid$highlight==1),]$age_condition))
    colnames(cell_num_all)<-c("age_condition","CellCount")
    fanciest_UMAP <- fanciest_UMAP + facet_wrap(~age_condition, ncol=2)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)
  }else{if(split==T){
    cell_num_all<-as.data.frame(table(d10x@meta.data$age_condition))
    colnames(cell_num_all)<-c("age_condition","CellCount")
    fanciest_UMAP <- fanciest_UMAP + facet_wrap(~age_condition, ncol=2)+  geom_text(aes(x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(CellCount))), cell_num_all, size=2)
  }else{if(!is.na(highlight)){
    fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(nrow(plt_myeloid[which(plt_myeloid$highlight==1),]))), size=2)
  }else{
    fanciest_UMAP <- fanciest_UMAP + annotate("text",x = min(plt_myeloid$UMAP_1)+(0.95*len_x_bar), y = min(plt_myeloid$UMAP_2)+(0.5*len_y_bar), label=paste0("n = ",comma(ncol(d10x))), size=2)
  }}}
  
  plot_grid(fanciest_UMAP,nice_legend, rel_widths = c(5,2))
  
  }




