library(caTools)
library(randomForest)

colnames(plt_varimax) <- paste0('varPC_', 1:ncol(plt_varimax))
identical(rownames(d10x.combined_myeloid@meta.data), rownames(plt_varimax))
plt_varimax$age_group <-as.factor(d10x.combined_myeloid@meta.data$AgeGroup)
table(plt_varimax$age_group)

### splitting the train and test data
sample = sample.split(plt_varimax$age_group, SplitRatio = .75)
train = subset(plt_varimax, sample == TRUE) 
test = subset(plt_varimax, sample == FALSE)

### training the RF model
age_group_pred_RF <- randomForest(age_group ~ . , data = train, importance = TRUE)


pred = predict(age_group_pred_RF, newdata=test[,-ncol(test)])
### generating a confusion matrix
cm = table(label=test[,ncol(test)], prediction=pred)
cm
gridExtra::grid.table(cm)
#### evaluating feature importance 
imp.df = data.frame(importance(age_group_pred_RF))        
imp.df[order(imp.df$MeanDecreaseAccuracy, decreasing = T),]
dev.off()
varImpPlot(age_group_pred_RF,main = 'age_group Prediction Based on Varimax-PCs')   

ggplot(plt_varimax_meta, aes(X1, X4, color=age_condition))+geom_point()+colscale_agecondition
ggplot(plt_varimax_meta, aes(age_condition, X4))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition+theme_bw()

ggplot(plt_varimax_meta, aes(X1, X4, color=age_condition))+geom_point()+colscale_agecondition
ggplot(plt_varimax_meta, aes(age_condition, X4))+geom_violin(fill="grey", color="grey")+geom_boxplot(aes(fill=age_condition))+fillscale_agecondition+theme_bw()



### PC4 loadings
gene_list = as.data.frame(unclass(df$rotLoadings))$PC_4
names(gene_list) = rownames(as.data.frame(unclass(df$rotLoadings)))
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]

res = GSEA(gene_list, GO_file, pval = 0.05)

plt_path<-res$Results
plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
plt_path$Enrichment_Cell<-"Down in PC5"
plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up in PC5"

plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))

plt_path$direction_label<-as.factor(plt_path$Enrichment)
levels(plt_path$direction_label)<-c(0.1,-0.1)
plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))

# top and bottom 15
plt_path<-rbind(plt_path[1:15,], plt_path[(nrow(plt_path)-15):(nrow(plt_path)),])

ggplot(plt_path, aes(NES, reorder(pathway, NES)))+geom_point(aes(size=size, fill=Enrichment_Cell), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_hline(yintercept=30.5, color="grey")+scale_fill_manual(values=c("#D64A56","cornflowerblue"))
