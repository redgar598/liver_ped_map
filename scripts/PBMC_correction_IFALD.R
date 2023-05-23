### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(SoupX)
library(colorspace)
library(cowplot)
library(DropletQC)
library(SCINA)




source("scripts/00_pretty_plots.R")
source("scripts/00_entropy_d10x.R")
source("scripts/00_fanciest_UMAP.R")


load("/media/redgar/Seagate Portable Drive/processed_data/IFALD_adult_ped_PBMC_integrated.rds")


load(here("data","IFALD_adult_ped_PBMC_integrated.rds"))

readRDS("/home/redgar/Documents/liver_ped_map/data/IFALD_adult_ped_PBMC_integrated.rds")
