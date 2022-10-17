# library(here)
# dataset_loc <- here("../../../projects/macparland/RE/PediatricAdult")
# print(list.files(dataset_loc))


library(targets)
library(tarchetypes)
tar_load(d10x.list.mt)
print(head(d10x.list.mt@meta.data))
