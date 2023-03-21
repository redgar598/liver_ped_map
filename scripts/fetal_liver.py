import anndata
adata = anndata.read("/media/redgar/Seagate Portable Drive/fetal_liver/download.h5ad")
adata.write("/media/redgar/Seagate Portable Drive/fetal_liver/fetal_liver.h5ad", compression="gzip")
