
cell_typeMeans_allindex<-read.csv(file="/media/redgar/Seagate Portable Drive/xenium_liver/ped_map_sc_means.csv")

colnames(cell_typeMeans_allindex)[which(colnames(cell_typeMeans_allindex)%in% c("HLA.DQA1","HLA.DQB1","HLA.DQB2"))]<-c('HLA-DQA1', 'HLA-DQB1', 'HLA-DQB2')

cell_typeMeans_allindex$UnassignedCodeword_0006<-0
cell_typeMeans_allindex$UnassignedCodeword_0037<-0
cell_typeMeans_allindex$UnassignedCodeword_0069<-0
cell_typeMeans_allindex$X<-NULL

cell_typeMeans_allindex<-cell_typeMeans_allindex[,c(1:477, 481,482,483, 478:480)]

write.csv(cell_typeMeans_allindex, file="/media/redgar/Seagate Portable Drive/xenium_liver/ped_map_sc_means_unassigned.csv", quote = F)
