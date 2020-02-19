#!/Local/md_shenorr/R-3.2.0/bin/Rscript
library(readr)
library(dplyr)
panther_orthologs <- select(read_delim("/storage/md_shenorr/inbaltz/panther_orthologs.csv", delim=','),c(2,3,4))

disease = Sys.getenv("disease")
RFST = Sys.getenv("RFST")

setwd(paste0("/storage/md_shenorr/inbaltz/microarrays/",RFST,"/",disease,"/"))
AllCompPairs_GS <- read_delim("AllCompPairs_GS", delim = '\t',col_names = c('index', 'human_file','mouse_file'), skip=1)
for (comp_file in 1:nrow(AllCompPairs_GS)){
  mousedf = readRDS(paste0(AllCompPairs_GS[comp_file,]$mouse_file,"_all_orthologs.rds"))
  humandf = readRDS(paste0(AllCompPairs_GS[comp_file,]$human_file,"_all_orthologs.rds"))
  colnames(mousedf) <- paste0("MM.",colnames(mousedf))
  colnames(humandf) <- paste0("HS.",colnames(humandf))
  merged <- merge(merge(panther_orthologs,mousedf, by.x = "mouse_entrez_ID", by.y = 0),humandf,by.x = "human_entrez_ID",by.y = 0)
  saveRDS(merged,paste0("comp_",AllCompPairs_GS[comp_file,]$mouse_file,"_",AllCompPairs_GS[comp_file,]$human_file,".rds"))
}
