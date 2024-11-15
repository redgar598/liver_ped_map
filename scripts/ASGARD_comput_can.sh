ssh redgar25@cedar.alliancecan.ca

cd scratch/liver_ped_map

salloc -c 8 -t 4:0:0 --mem 100G
module load r/4.4.0

R





sbatch scripts/ASGARD_ref.sh



scp redgar25@cedar.alliancecan.ca:/home/redgar25/scratch/L1000/DrugReference/*txt /media/redgar/Seagate\ Portable\ Drive/L1000/DrugReference
