
# Function to find optimal parameters
find_optimal_parameters <- function(x, y, z, points_per_maxima = 100) {
  n_points <- length(x)
  target_maxima <- n_points / points_per_maxima
  print(paste("target maxima: ",target_maxima, sep=""))
  best_params <- list(threshold = NA, k = NA, maxima_count = NA)
  best_diff <- Inf
  
  # Grid search over a range of thresholds and k values
  for (threshold in seq(5, 10, by = 1)) {
    for (k in seq(5, 30, by = 5)) {
      print(paste(threshold, k))
      local_maxima <- find_local_maxima(x, y, z, threshold, k)
      maxima_count <- nrow(local_maxima)
      print(maxima_count)
      diff <- abs(maxima_count - target_maxima)
      if (diff < best_diff) {
        best_diff <- diff
        best_params <- list(threshold = threshold, k = k, maxima_count = maxima_count)
      }
    }
  }
  
  return(best_params)
}



find_local_maxima <- function(x, y, z, threshold, k) {
  # Combine x and y into a matrix
  xy_matrix <- cbind(x, y)
  
  # Find the nearest neighbors for each point
  nearest_neighbors <- knn(xy_matrix, xy_matrix, z, k = k)
  
  # Initialize a vector to store indices of local maxima
  maxima_indices <- c()
  
  # Loop through each point
  for (i in 1:length(x)) {
    # Find the z value of the nearest neighbor
    nearest_z <- z[nearest_neighbors[i]]
    
    # Check if the z value is greater than the threshold
    if (z[i] > nearest_z + threshold) {
      maxima_indices <- c(maxima_indices, i)
    }
  }
  
  # Return the coordinates of the local maxima
  return(data.frame(x = x[maxima_indices], y = y[maxima_indices], z = z[maxima_indices]))
}


calculate_distance_to_maxima <- function(new_x, new_y, maxima_x, maxima_y) {
  distances <- c()
  
  for (i in 1:length(new_x)) {
    # Calculate distances to each local maximum
    dist_to_maxima <- sqrt((new_x[i] - maxima_x)^2 + (new_y[i] - maxima_y)^2)
    
    # Find the minimum distance
    min_dist <- min(dist_to_maxima)
    
    distances <- c(distances, min_dist)
  }
  
  return(distances)
}





























############################
### plot a gene expression in all samples
############################
plot_gene_expression<-function(gene, celltype){
  
  count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv")
  samples<-c("C94_2", "C94_3","C94_4","C105","C85")
  models<-c("2024_03_21_16_22_10","2024_03_21_16_38_48","2024_03_21_17_03_05","2024_05_17_13_32_03","2024_05_17_13_29_40")
  
  lapply(1:length(samples), function(x){
    data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
    tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
    smpl<-samples[x]
    print(smpl)
    
    ##### centroids
    tiff_res <- readTIFF(tiff_path, as.is = TRUE)
    tiff_res <- reshape2::melt(tiff_res)
    tiff_res <- tiff_res[tiff_res$value != 0, ]
    
    colnames(tiff_res) <- c("coord_y", "coord_x", "cell_id")
    tiff_res <- tiff_res[, c("coord_x", "coord_y", "cell_id")]
    tiff_res$cell_id <- as.numeric(tiff_res$cell_id)
    tiff_res$coord_x <- as.numeric(tiff_res$coord_x)
    tiff_res$coord_y <- as.numeric(tiff_res$coord_y)
    tiff_res <- tiff_res[order(tiff_res$cell_id), ]
    rownames(tiff_res) <- NULL
    
    centroids<-as.data.frame(
      tiff_res %>%
        group_by(cell_id) %>%
        summarise(centroid_x = sum(coord_x) / length(coord_x), centroid_y = sum(coord_y) / length(coord_y)))
    
    plt_umap_xenium_sample<-plt_umap_xenium[which(plt_umap_xenium$sample == smpl),]
    plt_umap_xenium_sample$cell<-sapply(1:nrow(plt_umap_xenium_sample), function(x) strsplit(plt_umap_xenium_sample$cell[x],"_")[[1]][1])
    
    cell_centroid<-merge(plt_umap_xenium_sample, centroids, by.x="cell",by.y="cell_id")
    
    ## gene expression
    counts<-read.csv(count_files[x])
    print(samples[x])
    counts$X<-NULL
    rownames(counts) <- counts$cell_id
    counts$cell_id <- NULL
    
    # Transpose the data and convert to sparse matrix.
    mat <- as(t(as.matrix(counts)), "sparseMatrix")
    
    seu <- CreateSeuratObject(counts=mat)
    seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
    seu
    
    
    cell_centroid<-cell_centroid[match(colnames(seu), cell_centroid$cell),]
    identical(colnames(seu), cell_centroid$cell)
    rownames(cell_centroid)<-cell_centroid$cell
    
    seu<-AddMetaData(seu, cell_centroid)
    seu <- NormalizeData(seu)
    
    gene_exp<-as.matrix(GetAssayData(object = seu, assay = "RNA", slot = "data"))
    
    meta<-seu@meta.data
    meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
    
    if(is.na(celltype)){
      gene_plt<-ggplot()+
        geom_point(aes(centroid_x,-centroid_y, color=exp),meta,size=0.75, shape=19)+
        theme_presentation()+ggtitle(smpl)+scale_color_gradientn(colours = plasma(10), name=gene,  na.value = "grey10")
      save_plts_black(gene_plt,paste(smpl,"_", gene,"_allcells" ,sep=""),  w=18, h=12)
    }else{
      gene_plt<-ggplot()+
        geom_point(aes(centroid_x,-centroid_y),meta[which(meta$CellType!="celltype"),],color="grey20",size=0.75, shape=19)+
        geom_point(aes(centroid_x,-centroid_y, color=exp),meta[which(meta$CellType==celltype),],size=0.75, shape=19)+
        theme_presentation()+ggtitle(smpl)+scale_color_gradientn(colours = plasma(10), name=gene,  na.value = "grey10")
      save_plts_black(gene_plt,paste(smpl,"_", gene, "_",celltype ,sep=""),  w=18, h=12)
    }
    
    })
  }



distance_gene_expression<-function(gene, distance ,celltype){
  count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")
  samples<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")
  models<-c("2024_03_21_16_22_10","2024_03_21_16_38_48","2024_03_21_17_03_05","2024_05_17_13_32_03","2024_05_17_13_29_40","2024_07_04_12_35_29","2024_06_28_09_23_40")
  
  load(here("data","zonation_allcells.RData"))
  zonation_scores<-do.call(rbind, zonation_scores)
  
  meta_exp<-lapply(2:length(samples), function(x){
    data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
    tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
    smpl<-samples[x]
    print(smpl)
    
    ## gene expression
    counts<-read.csv(count_files[x])
    print(samples[x])
    counts$X<-NULL
    rownames(counts) <- counts$cell_id
    counts$cell_id <- NULL
    
    # Transpose the data and convert to sparse matrix.
    mat <- as(t(as.matrix(counts)), "sparseMatrix")
    
    seu <- CreateSeuratObject(counts=mat)
    seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
    seu
    
    zonation_scores_sample<-zonation_scores[which(zonation_scores$sample==smpl),]
    zonation_scores_sample<-zonation_scores_sample[match(colnames(seu), zonation_scores_sample$cell),]
    identical(colnames(seu), zonation_scores_sample$cell)
    rownames(zonation_scores_sample)<-zonation_scores_sample$cell
    
    seu<-AddMetaData(seu, zonation_scores_sample)
    seu <- NormalizeData(seu)
    
    gene_exp<-as.matrix(GetAssayData(object = seu, assay = "RNA", slot = "data"))
    
    meta<-seu@meta.data
    meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
    
    meta})
  
  meta_exp<-do.call(rbind, meta_exp)
    
  if(distance=="Portal"){
    if(is.na(celltype)){
      distance_plt<- ggplot()+
        stat_smooth(aes(distance_to_periportal, exp, color=sample),meta_exp[which(!(meta_exp$CellType%in%c("Hepatocyte (Periportal)","Hepatocyte (Cycling)","Hepatocyte (Pericentral)"))),] ,method="lm", se=T,fill = "grey80")+
        facet_grid(CellType~., scale="free_y")+theme_bw()+
        theme(strip.background = element_rect(fill="white"))+
        xlab("Periportal Distance")+ylab(paste(gene, "Expression"))+
        scale_color_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+theme(strip.text.y = element_text(angle = 0))
      save_plts(distance_plt,paste(gene,"_periportal_allcells" ,sep=""),  w=5, h=12)
    }else{
      distance_plt<- ggplot()+
        stat_smooth(aes(distance_to_periportal, exp, color=sample),meta_exp[which((meta_exp$CellType==celltype)),] ,method="lm", se=T,fill = "grey80")+
        theme_bw()+
        theme(strip.background = element_rect(fill="white"))+
        xlab("Periportal Distance")+ylab(paste(gene, "Expression"))+
        scale_color_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+theme(strip.text.y = element_text(angle = 0))
      save_plts_black(distance_plt,paste(gene, "_periportal_",celltype ,sep=""),  w=6, h=6)}
  }else{
    if(is.na(celltype)){
      distance_plt<- ggplot()+
        stat_smooth(aes(distance_to_pericentral, exp, color=sample),meta_exp[which(!(meta_exp$CellType%in%c("Hepatocyte (Periportal)","Hepatocyte (Cycling)","Hepatocyte (Pericentral)"))),] ,method="lm", se=T,fill = "grey80")+
        facet_grid(CellType~., scale="free_y")+theme_bw()+
        theme(strip.background = element_rect(fill="white"))+
        xlab("Pericentral Distance")+ylab(paste(gene, "Expression"))+
        scale_color_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+theme(strip.text.y = element_text(angle = 0))
      save_plts(distance_plt,paste(gene,"_pericentral_allcells" ,sep=""),  w=5, h=12)
    }else{
      distance_plt<- ggplot()+
        stat_smooth(aes(distance_to_pericentral, exp, color=sample),meta_exp[which((meta_exp$CellType==celltype)),] ,method="lm", se=T,fill = "grey80")+
        theme_bw()+
        theme(strip.background = element_rect(fill="white"))+
        xlab("Pericentral Distance")+ylab(paste(gene, "Expression"))+
        scale_color_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+theme(strip.text.y = element_text(angle = 0))
      save_plts_black(distance_plt,paste(gene, "_pericentral_",celltype ,sep=""),  w=6, h=6)}
  }}



distance_gene_expression_combined<-function(gene, distance ,celltype){
  count_files<-c("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_2/cell_gene_matrices/2024_03_21_16_22_10/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_3/cell_gene_matrices/2024_03_21_16_38_48/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C94_4/cell_gene_matrices/2024_03_21_17_03_05/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C105/cell_gene_matrices/2024_05_17_13_32_03/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C85/cell_gene_matrices/2024_05_17_13_29_40/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C95/cell_gene_matrices/2024_07_04_12_35_29/expr_mat.csv",
                 "/media/redgar/Seagate Portable Drive/liver_BIDCell_output/C101/cell_gene_matrices/2024_06_28_09_23_40/expr_mat.csv")
  samples<-c("C94_2", "C94_3","C94_4","C105","C85","C95","C101")
  models<-c("2024_03_21_16_22_10","2024_03_21_16_38_48","2024_03_21_17_03_05","2024_05_17_13_32_03","2024_05_17_13_29_40","2024_07_04_12_35_29","2024_06_28_09_23_40")
  
  load(here("data","zonation_allcells.RData"))
  zonation_scores<-do.call(rbind, zonation_scores)
  
  meta_exp<-lapply(2:length(samples), function(x){
    data_dir<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x],"/cell_gene_matrices/", models[x], "/expr_mat.csv", sep="")
    tiff_path<-paste("/media/redgar/Seagate Portable Drive/liver_BIDCell_output/",samples[x], "/model_outputs/",models[x],"/test_output/epoch_1_step_4000_connected.tif", sep="")
    smpl<-samples[x]
    print(smpl)
    
    ## gene expression
    counts<-read.csv(count_files[x])
    print(samples[x])
    counts$X<-NULL
    rownames(counts) <- counts$cell_id
    counts$cell_id <- NULL
    
    # Transpose the data and convert to sparse matrix.
    mat <- as(t(as.matrix(counts)), "sparseMatrix")
    
    seu <- CreateSeuratObject(counts=mat)
    seu$sample<-strsplit(strsplit(count_files[x],"/")[[1]][6],"__")[[1]][3]
    seu
    
    zonation_scores_sample<-zonation_scores[which(zonation_scores$sample==smpl),]
    zonation_scores_sample<-zonation_scores_sample[match(colnames(seu), zonation_scores_sample$cell),]
    identical(colnames(seu), zonation_scores_sample$cell)
    rownames(zonation_scores_sample)<-zonation_scores_sample$cell
    
    seu<-AddMetaData(seu, zonation_scores_sample)
    seu <- NormalizeData(seu)
    
    gene_exp<-as.matrix(GetAssayData(object = seu, assay = "RNA", slot = "data"))
    
    meta<-seu@meta.data
    meta$exp<-gene_exp[which(rownames(gene_exp)==gene),]
    
    meta})
  
  meta_exp<-do.call(rbind, meta_exp)
  
  if(distance=="Portal"){
    if(is.na(celltype)){
      distance_plt<- ggplot()+
        stat_smooth(aes(distance_to_periportal, exp),meta_exp[which(!(meta_exp$CellType%in%c("Hepatocyte (Periportal)","Hepatocyte (Cycling)","Hepatocyte (Pericentral)"))),] ,method="lm", se=T,fill = "grey80")+
        facet_grid(CellType~., scale="free_y")+theme_bw()+
        theme(strip.background = element_rect(fill="white"))+
        xlab("Periportal Distance")+ylab(paste(gene, "Expression"))+
        scale_color_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+theme(strip.text.y = element_text(angle = 0))
      save_plts(distance_plt,paste(gene,"_periportal_allcells" ,sep=""),  w=5, h=12)
    }else{
      distance_plt<- ggplot()+
        stat_smooth(aes(distance_to_periportal, exp),meta_exp[which((meta_exp$CellType==celltype)),] ,method="lm", se=T,fill = "grey80")+
        theme_bw()+
        theme(strip.background = element_rect(fill="white"))+
        xlab("Periportal Distance")+ylab(paste(gene, "Expression"))+
        scale_color_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+theme(strip.text.y = element_text(angle = 0))
      save_plts_black(distance_plt,paste(gene, "_periportal_",celltype ,sep=""),  w=6, h=6)}
  }else{
    if(is.na(celltype)){
      distance_plt<- ggplot()+
        stat_smooth(aes(distance_to_pericentral, exp),meta_exp[which(!(meta_exp$CellType%in%c("Hepatocyte (Periportal)","Hepatocyte (Cycling)","Hepatocyte (Pericentral)"))),] ,method="lm", se=T,fill = "grey80")+
        facet_grid(CellType~., scale="free_y")+theme_bw()+
        theme(strip.background = element_rect(fill="white"))+
        xlab("Pericentral Distance")+ylab(paste(gene, "Expression"))+
        scale_color_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+theme(strip.text.y = element_text(angle = 0))
      save_plts(distance_plt,paste(gene,"_pericentral_allcells" ,sep=""),  w=5, h=12)
    }else{
      distance_plt<- ggplot()+
        stat_smooth(aes(distance_to_pericentral, exp),meta_exp[which((meta_exp$CellType==celltype)),] ,method="lm", se=T,fill = "grey80")+
        theme_bw()+
        theme(strip.background = element_rect(fill="white"))+
        xlab("Pericentral Distance")+ylab(paste(gene, "Expression"))+
        scale_color_manual(values=c("#D64A56","#ac3c27","cornflowerblue","#6db0dd"))+theme(strip.text.y = element_text(angle = 0))
      save_plts_black(distance_plt,paste(gene, "_pericentral_",celltype ,sep=""),  w=6, h=6)}
  }}


### add scale bar
xenium_scale_bar<-function(dataframe, size=4){
  list(annotate("rect", xmin=max(dataframe$centroid_x)-500, xmax=max(dataframe$centroid_x),
                ymin=min(-dataframe$centroid_y), ymax=(min(-dataframe$centroid_y))*0.995),
       annotate("text", x=max(dataframe$centroid_x)-250, y=(min(-dataframe$centroid_y))*1.01, label="0.5mm", size=size))
}


