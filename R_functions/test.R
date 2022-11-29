# library(here)
# dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")
# print(list.files(dataset_loc))

# 
# library(targets)
# library(tarchetypes)
# tar_make()
# 

# 
# tar_load(d10x.list.mt)
# print(d10x.list.mt)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MAST")

