


#### emprical p value
## only tested genes expressed in more that 15% of cell express


load(here("data/pval_zonation_allsamples_combined.RData"))


xenium<-read.csv(file=here("/home/redgar/Documents/xenium_liver/data/Xenium_CombinedPanel.csv"))


sig_de_age_RR<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_RR.csv"))
sig_de_age_KC<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_KC.csv"))
sig_de_age_MHCII<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_MHCII.csv"))



###### random lists of gene should be those tested for zonation in a cell type (ie expressed) then are there more significant than the equivalent number of genes from pediatric sinificant list

pval_distance_allsamples_kc<-pval_distance_allsamples[which(pval_distance_allsamples$celltype%in%c("KC Like")),]
# 134 genes tested for zonation in KC cells

ped_associated_tested_zonation<-pval_distance_allsamples_kc$gene[which(pval_distance_allsamples_kc$gene%in%sig_de_age_KC$X)]

pval_distance_allsamples_kc_sig<-pval_distance_allsamples_kc[which(pval_distance_allsamples_kc$p_adjusted_periportal<0.001),]


peds_zonated<-length(which(ped_associated_tested_zonation%in%pval_distance_allsamples_kc_sig$gene))
peds_zonated

rnd_overlaps<-sapply(1:10000, function(x) {
  set.seed(x)
  rnd_tested_zonation<-pval_distance_allsamples_kc$gene[sample(1:nrow(pval_distance_allsamples_kc), length(ped_associated_tested_zonation))]
  length(which(rnd_tested_zonation%in%pval_distance_allsamples_kc_sig$gene))
})

hist(rnd_overlaps)

(length(which(rnd_overlaps>=peds_zonated))+1)/(10000+1)
(length(which(rnd_overlaps<=peds_zonated))+1)/(10000+1)
