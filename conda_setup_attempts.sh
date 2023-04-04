
#Pandas 0.24.0 requires anndata version 0.6.18 and a scanpy version > 1.37.0

conda create -n fetal_liver python=2.7 pandas==0.23 
conda activate fetal_liver
conda install -c bioconda anndata
conda install -c conda-forge jupyter_core
conda install jupyter

conda install python=3.6.6

conda activate fetal_liver
pip install notebook
jupyter notebook

 


#cd /Documents/liver_ped_map/scripts
#python fetal_liver.py

import anndata
adata = anndata.read("/cluster/projects/macparland/RE/PediatricAdult/fetal_liver/download.h5ad")
adata.write("/cluster/projects/macparland/RE/PediatricAdult/fetal_liver/fetal_liver.h5ad", compression="gzip")
