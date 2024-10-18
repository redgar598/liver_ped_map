#################
## ASGARD
#################
#devtools::install_github("lanagarmire/Asgard")

library('Asgard')
library(Hmisc)


################
## download drug libraries
################


# on compute canada
# cd scratch
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz


PrepareReference(cell.info="/home/redgar25/scratch/L1000/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt",
                 gene.info="/home/redgar25/scratch/L1000/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt",
                 GSE70138.sig.info = "/home/redgar25/scratch/L1000/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt",
                 GSE92742.sig.info = "/home/redgar25/scratch/L1000/GSE92742_Broad_LINCS_sig_info.txt",
                 GSE70138.gctx = "/home/redgar25/scratch/L1000/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx",
                 GSE92742.gctx = "/home/redgar25/scratch/L1000/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx",
                 Output.Dir = "/home/redgar25/scratch/L1000/DrugReference/"
)
