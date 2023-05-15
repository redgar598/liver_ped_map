gene_df<-data.frame(celltype=c(rep("KC Like", length(kc_ifald)),rep("Macrophage\n(MHCII high)", length(mhcII_ifald)),rep("RR Myeloid", length(rr_ifald))),
                    gene=c(kc_ifald,mhcII_ifald,rr_ifald ))

gene<-gene_df$gene


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



plt_exp_summary<-merge(plt_exp_summary, gene_df, by.x="variable", by.y="gene")

plt_exp_summary$variable<-factor(plt_exp_scaled$variable, levels=rev(gene))

plot_grid(
  ggplot(gene_df, aes(1,gene, fill=celltype))+
    geom_tile(color="black")+facet_grid(celltype~., scales = "free_y", space = "free_y",switch="both")+
    th+theme_classic()+fillscale_cellType+scale_x_continuous(position="top")+
    ylab("")+xlab("Significant\nCell Type")+theme(legend.position = "none",
                            strip.background = element_blank(),
                            strip.text.y.left = element_text(angle = 0),
                            axis.text=element_blank(),axis.ticks = element_blank(),axis.line =element_blank()),
  ggplot(plt_exp_summary, aes( age_condition,variable, fill=scaled))+
    geom_tile()+facet_grid(celltype~CellType_refined, scales = "free_y", space = "free_y")+
    th+theme_classic()+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(8, 'RdBu')), space='Lab')(100), name="Scaled\nMean\nExpression")+
    ylab("")+xlab("")+theme(panel.spacing = unit(0,"line"),
                            strip.text.y = element_blank()),
  align = "v", axis = "tb", rel_widths = c(1,10))
