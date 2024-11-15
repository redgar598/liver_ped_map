library(ggplot2)
library(here)
library(dplyr)
library(scales)
library(reshape2)
library(cowplot)

source("scripts/00_pretty_plots.R")


## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

### meta data
metas_tracking<-read.csv(here("data/data_transfer_updated_feb12_2024_IFALD_PBMC.csv"))

## ancestry file location
dataset_loc <- here("../ancestry_calling/scadmix/pediatric")

samples<-c(list.files(dataset_loc))
print(samples)

samples<-samples[grep("_ancestry.qopt", samples)]
samples<-gsub("_ancestry.qopt","",samples)
print(samples)

scadmix_ancestry<-do.call(rbind, lapply(1:length(samples), function(x){
  ancestry<-read.table(paste(dataset_loc, "/",samples[x], "_ancestry.qopt", sep=""), header=T)
  ancestry$indivID<-samples[x]
  ancestry}))

scadmix_ancestry<-merge(metas_tracking[,c("Sample_ID","file","Age")], scadmix_ancestry,by.x="file", by.y="indivID")

scadmix_ancestry$donor<-sapply(1:nrow(scadmix_ancestry),function(x) strsplit(scadmix_ancestry$Sample_ID[x], "_")[[1]][1])
scadmix_ancestry$Age<-as.character(scadmix_ancestry$Age)


############
## Plot
############
scadmix_proportions<-melt(scadmix_ancestry)
scadmix_proportions$label<-sapply(1:nrow(scadmix_proportions), function(x) paste(scadmix_proportions$donor[x],"-",scadmix_proportions$Age[x],  sep=" "))



nice_legend<-get_leg( ggplot(scadmix_proportions, aes(fill=variable, x=value, y=label)) + 
                        geom_bar(position="fill", stat="identity", color="black")+
                        scale_fill_manual(values=c("#fdb462","#b2df8a","#8dd3c7","#80b1d3","#cc689b","#fb9a99","#ac8ed1"), name="HGDP Superpopulation") +
                        theme_bw()+ylab("Donor ID - Age")+xlab("Proportion of genome attributed to each ancestral region"))

  
ancestry_proportions <- plot_grid(
  ggplot(scadmix_proportions, aes(fill=variable, x=value, y=label)) + 
    geom_bar(position="fill", stat="identity", color="black")+
    scale_fill_manual(values=c("#fdb462","#b2df8a","#8dd3c7","#80b1d3","#cc689b","#fb9a99","#ac8ed1"), name="HGDP Superpopulation") +
    theme_bw()+ylab("Donor ID - Age")+xlab("Proportion of genome attributed to each ancestral region")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.text.y = element_text(size = 12),
          legend.position = "none"),
  nice_legend, ncol=2, rel_widths = c(4,1))

ancestry_proportions

save_plts(ancestry_proportions, "scadmix_pediatric", w=10, h=4)

