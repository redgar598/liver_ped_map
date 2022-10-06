library(targets)
library(tarchetypes)


source("R_functions/functions.R")
tar_option_set(packages = c("here","dplyr", "ggplot2","scales","Seurat"))


list(
  tar_target(load_d10x_raw, load_d10x_raw(here("../../../projects/macparland/RE/PediatricAdult"))),
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