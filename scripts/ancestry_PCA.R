library(ggplot2)
library(here)
library(dplyr)
library(scales)
library(reshape2)

source("scripts/00_pretty_plots.R")

args = commandArgs(trailingOnly=TRUE)
vega_20 = muted(c(
  '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
  '#8C564B', '#E377C2', '#BCBD22', '#17BECF', '#AEC7E8',
  '#FFBB78', '#98DF8A', '#FF9896', '#C5B0D5', '#C49C94',
  '#F7B6D2', '#DBDB8D', '#9EDAE5', '#AD494A', '#8C6D31'),l = 70, c = 40)

vega_20 <- c(vega_20[4], vega_20[1], vega_20[2], vega_20[3], vega_20[seq(5, length(vega_20),1)])

ref <- read.csv(file="../ancestry_calling/Monopogen/resource/HGDP.PC.csv",head=T)
colnames(ref)[1]<-"PopID"
ref <- ref[,seq(1,5,1)]
ref$PC3<-NULL

ref_pop <- read.csv(file="../ancestry_calling/igsr-human genome diversity project.tsv",head=T, sep="\t")
print(ref_pop %>% group_by(Superpopulation.name) %>%  mutate(value_string = paste(unique(Population.description), collapse = ",")) %>% select(Superpopulation.name, value_string) %>% unique(), n=50)
subpops<-as.data.frame(ref_pop %>% group_by(Superpopulation.name) %>%  mutate(value_string = paste(unique(Population.description), collapse = ",")) %>% select(Superpopulation.name, value_string) %>% unique(), n=50)
subpops$value_string<-sapply(1:nrow(subpops), function(x) gsub("\\(|\\)|HGDP","", subpops$value_string[x]))
subpops$Superpopulation.name<-sapply(1:nrow(subpops), function(x) gsub("\\(|\\)|HGDP","", subpops$Superpopulation.name[x]))
write.table(subpops, file=here("../ancestry_calling/HGDP_subpopulations.txt"), quote=F, row.names = F, sep="\t")




samples=c(
  "C93",
  "C64",
  "C104",
  "C105",
  "C39_NPC",
  "C39_TLH",
  "C54",
  "C70",
  "C82",
  "C85",
  "C88",
  "C96",
  "IFALD006",
  "IFALD030",
  "IFALD073_PBMC",
  "IFALD073",
  "C102", "C97",
  "P1TLH",
  "P2TLH",
  "P3TLH",
  "P5TLH",
  "C61"#"C86",
  )




PC_samples<-lapply(samples, function(smple){
  print(smple)
  df <- read.table(file=paste("../ancestry_calling/Monopogen/",smple,"_trace.ProPC.coord", sep=""),head=T,sep="\t")
  df$PopID<-"scRNAseq_liver"
  df$indivID<-smple
  df
})

study_coord<-do.call(rbind, PC_samples)

study_coord<-study_coord[,c("PopID","indivID","PC1" ,"PC2")]

dt_plt <- rbind(ref,study_coord)

ancestry_PCA<-ggplot() + geom_point(aes(x=PC1, y=PC2, color=PopID),dt_plt, size=2) +  
  scale_color_manual(values=c("#b2df8a","#fb9a99","#fdb462","#80b1d3","#8dd3c7","#cc689b","#ac8ed1","grey30")) + theme_bw() + 
  guides(colour=guide_legend(override.aes = list(alpha=0.9)))+
  geom_text(aes(x=PC1, y=PC2, label=indivID),dt_plt[which(dt_plt$PopID=="scRNAseq_liver"),], hjust=-0.25, vjust=0, size=2)
save_plts(ancestry_PCA, "ancestry_ped_PCA", w=7,h=5)


###########
## 5 liver map too


samples=c(
  "C93",
  "C64",
  "C104",
  "C105",
  "C39_NPC",
  "C39_TLH",
  "C54",
  "C70",
  "C82",
  "C85",
  "C88",
  "C96",
  "IFALD006",
  "IFALD030",
  "IFALD073_PBMC",
  "IFALD073",
  "P1TLH",
  "P2TLH",
  "P3TLH",
  "P5TLH",
  "C61"
) #"C102", "C86","C97", "P4TLH


PC_samples<-lapply(samples, function(smple){
  print(smple)
  df <- read.table(file=paste("../ancestry_calling/Monopogen/",smple,"_trace.ProPC.coord", sep=""),head=T,sep="\t")
  df$PopID<-"scRNAseq_liver"
  df$indivID<-smple
  df
})

study_coord<-do.call(rbind, PC_samples)

study_coord<-study_coord[,c("PopID","indivID","PC1" ,"PC2")]

dt_plt <- rbind(ref,study_coord)

ancestry_PCA<-ggplot() + geom_point(aes(x=PC1, y=PC2, color=PopID),dt_plt, size=2) +  
  scale_color_manual(values=c(vega_20[1:7],"grey60")) + theme_bw() + 
  guides(colour=guide_legend(override.aes = list(alpha=0.9)))+
  geom_text(aes(x=PC1, y=PC2, label=indivID),dt_plt[which(dt_plt$indivID%in%c(  "P1TLH",
                                                                                "P2TLH",
                                                                                "P3TLH",
                                                                                "P5TLH",
                                                                                "C61")),], hjust=-0.25, vjust=0, size=2)
save_plts(ancestry_PCA, "ancestry_ped_5liverandC61_PCA", w=7,h=5)



###########
## 5 liver map and C61 only
###########
samples=c(
  "P1TLH",
  "P2TLH",
  "P3TLH",
  "P5TLH",
  "C61"
) #"C102", "C86","C97", "P4TLH


PC_samples<-lapply(samples, function(smple){
  print(smple)
  df <- read.table(file=paste("../ancestry_calling/Monopogen/",smple,"_trace.ProPC.coord", sep=""),head=T,sep="\t")
  df$PopID<-"scRNAseq_liver"
  df$indivID<-smple
  df
})

study_coord<-do.call(rbind, PC_samples)

study_coord<-study_coord[,c("PopID","indivID","PC1" ,"PC2")]

dt_plt <- rbind(ref,study_coord)

dt_plt$label<-NA
dt_plt$label[which(dt_plt$indivID=="P1TLH")]<-"P1TLH\n(Not Reported)"
dt_plt$label[which(dt_plt$indivID=="P2TLH")]<-"P2TLH\n(White)"
dt_plt$label[which(dt_plt$indivID=="P3TLH")]<-"P3TLH\n(White)"
dt_plt$label[which(dt_plt$indivID=="P5TLH")]<-"P5TLH\n(Black)"
dt_plt$label[which(dt_plt$indivID=="C61")]<-"C61\n(South East Asian)"

ancestry_PCA<-ggplot() + geom_point(aes(x=PC1, y=PC2, color=PopID),dt_plt, size=1.5) +  
  scale_color_manual(values=c("#b2df8a","#fb9a99","#fdb462","#80b1d3","#8dd3c7","#cc689b","#ac8ed1","grey30")) + theme_bw() + 
  guides(colour=guide_legend(override.aes = list(alpha=0.9)))+xlim(-250, 650)+ylim(-300, 550)+
  geom_text(aes(x=PC1, y=PC2, label=label),dt_plt[which(dt_plt$indivID%in%c(  "P1TLH",
                                                                              "P2TLH",
                                                                              "P3TLH",
                                                                              "P5TLH",
                                                                              "C61")),], hjust=0, vjust=-0.25, size=2)
ancestry_PCA
save_plts(ancestry_PCA, "ancestry_5liverandC61_PCA", w=7,h=5)




###############
## scadmix output
###############
P1TLH<-read.table(here("../ancestry_calling/scadmix/liver_samples/P1TLH_ancestry.qopt"), header=T)
P1TLH$indivID<-"P1TLH"
P2TLH<-read.table(here("../ancestry_calling/scadmix/liver_samples/P2TLH_ancestry.qopt"), header=T)
P2TLH$indivID<-"P2TLH"
P3TLH<-read.table(here("../ancestry_calling/scadmix/liver_samples/P3TLH_ancestry.qopt"), header=T)
P3TLH$indivID<-"P3TLH"
P5TLH<-read.table(here("../ancestry_calling/scadmix/liver_samples/P5TLH_ancestry.qopt"), header=T)
P5TLH$indivID<-"P5TLH"
C61<-read.table(here("../ancestry_calling/scadmix/liver_samples/C61_ancestry.qopt"), header=T)
C61$indivID<-"C61"

scadmix_proportions<-rbind(P1TLH,P2TLH, P3TLH, P5TLH, C61)
scadmix_proportions<-melt(scadmix_proportions)

scadmix_proportions$label<-NA
scadmix_proportions$label[which(scadmix_proportions$indivID=="P1TLH")]<-"P1TLH\n(Not Reported)"
scadmix_proportions$label[which(scadmix_proportions$indivID=="P2TLH")]<-"P2TLH\n(White)"
scadmix_proportions$label[which(scadmix_proportions$indivID=="P3TLH")]<-"P3TLH\n(White)"
scadmix_proportions$label[which(scadmix_proportions$indivID=="P5TLH")]<-"P5TLH\n(Black)"
scadmix_proportions$label[which(scadmix_proportions$indivID=="C61")]<-"C61\n(South East Asian)"

ancestry_proportions<-ggplot(scadmix_proportions, aes(fill=variable, y=value, x=label)) + 
  geom_bar(position="fill", stat="identity", color="black")+
  scale_fill_manual(values=c("#fdb462","#b2df8a","#ac8ed1","#8dd3c7","#cc689b","#fb9a99","#80b1d3")) +
  theme_bw() +xlab("Individual ID\n(Self Reported Ethnicity)")+ylab("Proportion of genome attributed to each ancestral region")
ancestry_proportions
save_plts(ancestry_proportions, "ancestry_scadmix_5liverandC61_PCA", w=9,h=6)


