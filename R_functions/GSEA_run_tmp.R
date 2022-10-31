library(dplyr)
#http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/Human_GOBP_AllPathways_no_GO_iea_October_26_2022_symbol.gmt
GO_file = here("data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")


top10_pathways<-do.call(rbind,lapply(0:6, function(cluster){
  diff_exp_all_cluster<-diff_exp_all[which(diff_exp_all$cluster==cluster),]
  print(cluster)

  gene_list = diff_exp_all_cluster$avg_log2FC
  names(gene_list) = diff_exp_all_cluster$gene
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]

  res = GSEA(gene_list, GO_file, pval = 0.05)
  
  data.frame(cluster=cluster, pathway=sapply(1:10, function(x) strsplit(res$Results$pathway[x], "%")[[1]][1]), direction=sapply(1:10, function(x) res$Results$Enrichment[x]))}))
  

res$Plot
View(top10_pathways)
