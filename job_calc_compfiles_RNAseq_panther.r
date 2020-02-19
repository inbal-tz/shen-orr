#!/Local/md_shenorr/R-3.2.0/bin/Rscript
library(readr)
library(dplyr)
panther_orthologs <- select(read_delim("/storage/md_shenorr/inbaltz/panther_orthologs.csv", delim=','),c(2,3,4))

disease = Sys.getenv("disease")
RFST = Sys.getenv("RFST")

setwd(paste0("/storage/md_shenorr/rachelly/Shared_Dir/GS_analysis/RNAseq/",RFST,"/",disease,"/V1.9"))
AllCompPairs_GS <- read_delim("AllCompPairs_GS", delim = '\t',col_names = c('index', 'human_file','mouse_file'), skip=1)

#these are rachelly's files for converting ensemble to entrezid
HS_ensemble_entrez <- read_delim("/storage/md_shenorr/rachelly/Gene_PWs_Analysis/AuxFiles/Orthologs/HS_Ensembl_Entrez.txt",delim='\t')
MM_ensemble_entrez <- read_delim("/storage/md_shenorr/rachelly/Gene_PWs_Analysis/AuxFiles/Orthologs/MM_Ensembl_Entrez.txt",delim='\t')


for (comp_file in 1:nrow(AllCompPairs_GS)){
  setwd(paste0("/storage/md_shenorr/rachelly/Shared_Dir/GS_analysis/RNAseq/",RFST,"/",disease,"/V1.9"))
  mousedf = select(merge(readRDS(paste0("Step_2_DEGs/sleuth_",AllCompPairs_GS[comp_file,]$mouse_file,".rds")),MM_ensemble_entrez, by.x = "target_id",by.y = "Gene ID"),c(3,4,12))
  humandf = select(merge(readRDS(paste0("Step_2_DEGs/sleuth_",AllCompPairs_GS[comp_file,]$human_file,".rds")),HS_ensemble_entrez, by.x = "target_id",by.y = "Gene ID"),c(3,4,12))
  #b column (beta) is used as Ztest!
  names(mousedf) <- c("qval","Ztest","EntrezGene ID")
  names(humandf) <- c("qval","Ztest","EntrezGene ID")
  colnames(mousedf) <- paste0("MM.",colnames(mousedf))
  colnames(humandf) <- paste0("HS.",colnames(humandf))
  merged <- merge(merge(panther_orthologs,mousedf, by.x = "mouse_entrez_ID", by.y = "MM.EntrezGene ID"),humandf,by.x = "human_entrez_ID",by.y = "HS.EntrezGene ID")
  #save in inbal's directory
  setwd(paste0("/storage/md_shenorr/inbaltz/RNAseq/",RFST,"/",disease))
  saveRDS(merged,paste0("comp_",AllCompPairs_GS[comp_file,]$mouse_file,"_",AllCompPairs_GS[comp_file,]$human_file,".rds"))
}
