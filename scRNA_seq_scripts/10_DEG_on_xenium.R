xenium<-read.csv(file=here("/home/redgar/Documents/xenium_liver/data/Xenium_CombinedPanel.csv"))


sig_de_age_RR<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_RR.csv"))
sig_de_age_KC<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_KC.csv"))
sig_de_age_MHCII<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_age_MHCII.csv"))

sig_de_age_KC[which(sig_de_age_KC$X%in%xenium$Gene),]
sig_de_age_RR[which(sig_de_age_RR$X%in%xenium$Gene),]
sig_de_age_MHCII[which(sig_de_age_MHCII$X%in%xenium$Gene),]

intersect(xenium$Gene, sig_de_age_KC$X)
intersect(xenium$Gene, sig_de_age_RR$X)
intersect(xenium$Gene, sig_de_age_MHCII$X)



sig_de_IFALD_RR<-read.csv( file=here("/home/redgar/Documents/liver_ped_map/data","differential_IFALD_RR.csv"))
sig_de_IFALD_KC<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_IFALD_KC.csv"))
sig_de_IFALD_MHCII<-read.csv(file=here("/home/redgar/Documents/liver_ped_map/data","differential_IFLAD_MHCII.csv"))


sig_de_IFALD_KC[which(sig_de_IFALD_KC$X%in%xenium$Gene),]
sig_de_IFALD_RR[which(sig_de_IFALD_RR$X%in%xenium$Gene),]
sig_de_IFALD_MHCII[which(sig_de_IFALD_MHCII$X%in%xenium$Gene),]
