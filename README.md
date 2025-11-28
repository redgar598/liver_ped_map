# A Single-Cell Atlas Of Human Pediatric Liver Reveals Age-Related Hepatic Gene Signatures

Analysis to build the map of the pediatric liver and compared healthy pediatric liver to IFALD

# Workflow
Analysis is seperated into dissociated single cell analysis and the xenium spatial analysis. The scripts were run in the order indicated by the numbers in each pipeline with 00 being scripts sourced throughout the workflow (common plots and colour schemes)

The data folder referenced in the analysis is not on github. However, all data is deposited pubically.

# Data Links
**Dissociated single cell map**

[CELLxGENE](https://cellxgene.cziscience.com/collections/ff69f0ee-fef6-4895-9f48-6c64a68c8289)

[Seurat objects](https://drive.google.com/drive/folders/1lk76_1P8Jo0TeaXy4hpoxEqTt1nxL9r_?usp=sharing)

[Raw data](https://explore.data.humancellatlas.org/projects/febdaddd-ad3c-4f4a-820f-ade15c48545a)




**Xenium Spatial Transcriptomics**

[Shiny application](https://macparlandlab.shinyapps.io/pediatric_liver_spatial/)

[Raw data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286254)

** Xenium Annotation**
For the cell type annotation of the Xenium presented in the paper we used the output of BIDCell. The output of the BIDCell segmentation and the annotation of cell type and zonation (with segmentation centroids) are available [here](https://drive.google.com/file/d/1rARs4UAhtwmsFD2uGewnyLTbLtHI0ZmD/view?usp=sharing)

To load this output as a Seurat object we have provided some sample [code](https://github.com/redgar598/liver_ped_map/blob/main/xenium_scripts/BIDCell_xenium_Seurat.R)
