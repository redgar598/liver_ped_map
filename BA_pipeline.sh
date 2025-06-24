
scp -r /home/redgar/Documents/liver_ped_map/data/BA_GSE176189 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/
scp -r /home/redgar/Documents/liver_ped_map/data/Taylor_GSE163650 t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/


salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/PediatricAdult


cd liver_ped_map
## on h4h
salloc -c 1 -t 3:0:0 --mem 25G

cd /cluster/home/t117652uhn/liver_ped_map
module load R/4.2.1



sbatch BA_integration.sh



../../../projects/macparland/RE/PediatricAdult/BA_GSE176189","BA_Taylor_d10x_adult_ped_integrated.rds

scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/PediatricAdult/BA_GSE176189/BA_Taylor_d10x_adult_ped_integrated.rds /media/redgar/Seagate\ Portable\ Drive

scp /media/redgar/Seagate\ Portable\ Drive/BA_Taylor_d10x_adult_ped_integrated.rds redgar@rc01.ccbr.utoronto.ca:/home/baderlab/redgar

scp redgar@rc01.ccbr.utoronto.ca:/home/baderlab/redgar/BA_Taylor_d10x_adult_ped_HSC.rds /home/redgar/Documents/liver_ped_map/data
scp redgar@rc01.ccbr.utoronto.ca:/home/baderlab/redgar/BA_Taylor_d10x_adult_ped_cholangiocytes.rds /home/redgar/Documents/liver_ped_map/data
