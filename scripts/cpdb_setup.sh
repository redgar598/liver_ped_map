conda create -n cpdb python=3.8
conda activate cpdb 
pip install cellphonedb





scp /home/redgar/Documents/liver_ped_map/data/cellphonedb/v4.1.0/cellphonedb.zip t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/home/t117652uhn/liver_ped_map/data



cd liver_ped_map
git pull
sbatch cpdb_liver.sh


nano cpdb_ped_healthy.out