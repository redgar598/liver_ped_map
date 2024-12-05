save_plts<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h)
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h, bg = "white")
  ggsave(plt, file=paste(here("figures/png/"),name,".png", sep=""), w=w, h=h, bg = "white")}

save_plts_black<-function(plt, name, w,h){
  ggsave(plt, file=paste(here("figures/"),name,".pdf", sep=""), w=w, h=h, bg = "black")
  ggsave(plt, file=paste(here("figures/jpeg/"),name,".jpeg", sep=""), w=w, h=h, bg = "black")
  ggsave(plt, file=paste(here("figures/png/"),name,".png", sep=""), w=w, h=h, bg = "black")}

theme_presentation<- function(base_size = 15, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text = element_blank(),
      axis.ticks =  element_blank(), 
      axis.title= element_blank(),
      panel.background = element_rect(fill="black"), 
      panel.border =element_blank(),  
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.margin = unit(1.0, "lines"), 
      plot.background = element_rect(fill="black"), 
      plot.title =element_text(size=15,colour="white"), 
      plot.margin = unit(c(1,  1, 1, 1), "lines"),
      legend.background=element_rect(fill='black'),
      legend.title=element_text(size=15,colour="white"),
      legend.text=element_text(size=15,colour="white"),
      legend.key = element_rect( fill = 'black'),
      legend.key.size = unit(c(1, 1), "lines"),
      axis.line = element_blank()
    )
}

th_present <- theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    strip.text = element_text(size = 12),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14))



## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


######### 
## Cell Type Color
######### 
myColors_celltype <- c("#2b8cbe", "#7a0177","#7bccc4",  
                       "#cccc16","#f0efa8","#f0efa8",
                       "#a63603","#f03b20","#fcbba1",
                       "#0a15f2","#3e44b8","#5d64f5","#929bda","#5d64f5",
                       "#058205","#8e11bf","#67038f",           
                       "#d1100d","#7a4870", "#27751e",  
                       "#8ede85", "#509418",
                       "#ad235c","#f781b0","#ed0ce2","#ad235c","grey",
                       "hotpink","green","#e0a8ce")  


color_possibilities_celltype<-c("Hepatocyte (Periportal)", "Hepatocyte (Pericentral)", "Hepatocyte (Cycling)",  
                                "Cholangiocytes","Cholangiocyte (Biliary)","Cholangiocytes (Biliary tree)",
                                "LSEC","LSEC II","VEC",
                                "HSC (Periportal)","HSC (Quiescent)","HSC","HSC (Activated)","Mesenchyme Low Gene Count (Biliary tree)",
                                "NK-like cells","Mature B-cells","Plasma Cells",           
                                "Erythrocytes","Neutrophil", "gd T-cells",  
                                "Cycling (T Cell)", "CD3+ T-cells",
                                "Macrophage MHCII High","KC Like","Mono-Mac","Myeloid","Doublet",
                                "Check Posistion Hepatocyte (Periportal)","Check Posistion Hepatocyte (middle)","Cycling Myeloid")  
names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)
colscale_cellType <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)


 
#############
## ped map cell labels
############

myColors_celltype_ped <- c("#660cc7","#8642d1","#5612a3","#4b911d","#7a4202",
                       "#2c6e02","#6bbce8","#3469ad","#d17906","#b01e87",
                       "#60ba5d","#207537","#a0c487","#d9a5a5","#87a4c9",
                       "#e8c392","#dea4ce","#79639a","#207537","#fa61ad",
                       "#b80783","#994676","#431039","#cb181d","maroon1",
                       "#b01629","grey","#ce1256","#a6d96a","#750c32",
                       "#d9667f","#1b4003","#e0a8ce","#8a68b0","#3d1b63",
                       "#c9a8ed","#c48db4","#a3588d","#6ca647","#3a7d31",
                       "#60ba5d","#3d1b63","#2dc918","#F4355B","#e4d5f2")
color_possibilities_celltype_ped<-c("B-cells","Mature B-cells","Plasma cells","CD3+ T-cells","Cholangiocytes",
                                "gd T-cells","Hepatocytes","HSC","LSEC","Myeloid cells",
                                "NK-like cells", "NK and T cells","NKT cells\n(Hepatocyte Like)","Cholangiocytes\n(Hepatocyte Like)","HSC\n(Hepatocyte Like)",
                                "LSEC\n(Hepatocyte Like)","Myeloid cells\n(Hepatocyte Like)","B-cells\n(Hepatocyte Like)","NKT cells","Mono-Mac",
                                "KC Like","Neutrophil","Neutrophil\n(DEFA+)","Erythrocytes","Mast cell",
                                "Myeloid Erythrocytes\n(phagocytosis)","Doublet","Macrophage\n(MHCII high)","Cycling T-cells","CDC1",
                                "Platelets","CLNK T-cells","Cycling Myeloid","Mature B-cells (High MT)","pDC",
                                "pre B-cell","KC Like\n(Hepatocyte Like)","KC Like (C97)","Naive CD4 T-cells","Memory CD4 T-cells",
                                "NK cells","DC","CD8 T-cells","CD14+ Mono","Cycling Plasma")
names(myColors_celltype_ped) <- color_possibilities_celltype_ped
fillscale_cellType_ped <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype_ped, drop = T, limits=force)
colscale_cellType_ped <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype_ped, drop = T, limits=force)




myColors_age <- c("#D64A56","cornflowerblue")
color_possibilities_age<-c( "Adult","Pediatric")
names(myColors_age) <- color_possibilities_age
fillscale_age <- scale_fill_manual(name="Age\nGroup",
                                   values = myColors_age, drop = T, limits=force)
colscale_age <- scale_color_manual(name="Age\nGroup",
                                   values = myColors_age, drop = T, limits=force)

       

######### 
## highlight cell type
######### 
highlight_cell<-function(fav_celltype){
  cell_types<-unique(plot_freely$CellType)
  myColors_celltype <- c(rep("grey20", length(cell_types)))
  myColors_celltype[which(cell_types==fav_celltype)]<-"green"

  names(myColors_celltype) <- cell_types
  scale_color_manual(name="Cell Type",values = myColors_celltype, drop = T, limits=force)}


#################
## Density of a cell type
#################
cell_density<-function(bin_number,min_cell_in_bin ,subset){
  
  proportion_overall<-length(plot_freely$x[which(plot_freely$CellType==subset)])/nrow(plot_freely)
  print(paste(round(proportion_overall,3), subset, "across whole section"))
  
  # Bin width
  bin_width_x <- max(plot_freely$x)/bin_number
  bin_width_y <- max(plot_freely$y)/bin_number
  
  # Create binning grid
  x_bins <- seq(0, max(plot_freely$x), by = bin_width_x)
  y_bins <- seq(0, max(plot_freely$y), by = bin_width_y)
  bins <- list(x = x_bins, y = y_bins)
  
  # Cut points into bins
  x_cut <- cut(plot_freely$x, breaks = x_bins, include.lowest = TRUE)
  y_cut <- cut(plot_freely$y, breaks = y_bins, include.lowest = TRUE)
  
  # Create a table of bin counts
  point_counts_all <- table(x_cut, y_cut)
  point_counts_all[point_counts_all < min_cell_in_bin] <- NA
  
  # Cut points into bins  subset
  x_cut_subset <- cut(plot_freely$x[which(plot_freely$CellType==subset)], breaks = x_bins, include.lowest = TRUE)
  y_cut_subset <- cut(plot_freely$y[which(plot_freely$CellType==subset)], breaks = y_bins, include.lowest = TRUE)
  
  point_counts_subset <- table(x_cut_subset, y_cut_subset)
  
  # Display the point counts in each bin
  density_difference<-as.matrix(point_counts_subset) / as.matrix(point_counts_all)
  colnames(density_difference)<-1:bin_number
  rownames(density_difference)<-1:bin_number
  density_plt<-as.data.frame(density_difference)
  
  
  ggplot(density_plt, aes(x=x_cut_subset, y=rev(y_cut_subset), fill=Freq))+geom_tile( color="black")+
    theme(
      axis.text = element_blank(),
      axis.ticks =  element_blank(), 
      axis.title= element_blank(),
      panel.border =element_blank(),  
      legend.title=element_text(size=15,colour="black"),
      legend.text=element_text(size=15,colour="black")
    ) + 
    scale_fill_gradientn(
      colours=c("#0C6291", "#ebe8e8", "#A63446"), 
      na.value = "grey60", name = "Freq",
      values = c(0,  proportion_overall, max(density_plt$Freq, na.rm=T))/max(density_plt$Freq, na.rm=T)
    )
  
}


