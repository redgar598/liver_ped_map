BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("Orthology.eg.db")

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

mapIt <- function(mouseids, horg, morg, orth){
  mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
  mapped <- select(orth, mouseg, "Homo_sapiens","Mus_musculus")
  names(mapped) <- c("Mus_egid","Homo_egid")
  husymb <- select(horg, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  return(data.frame(Mus_symbol = mouseids,
                    mapped,
                    Homo_symbol = husymb[,2]))
}



mapIt(musGenes, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
