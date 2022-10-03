library(targets)
library(tarchetypes)


source("R_functions/functions.R")
tar_option_set(packages = c("here","dplyr", "ggplot2","scales","Seurat"))


list(
  tar_target(file, "data/GSE136831_mini.RData", format = "file"),
  tar_target(d10x, load_d10x(file)),
  tar_target(d10x_MT, MT(d10x)),
  tar_target(d10x_QC, QC(d10x_MT)),
  tar_render(report, "R_functions/visualize_functions.Rmd")
  )

#' library(dplyr)
#' library(Seurat)
#' library(patchwork)
#' library(here)
#' library(ggplot2)
#' library(reshape2)
#' library(gridExtra)
#' library(limma)
#' library(cowplot)
#' library(gtools)
#' library(ggsignif)