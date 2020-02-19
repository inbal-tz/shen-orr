#!/Local/md_shenorr/R-3.2.0/bin/Rscript
#connect all comp files
library(readr)
library(dplyr)

###############################################  microarrays  ###############################################

supp_table = read.table("/storage/md_shenorr/inbaltz/Supplementary Table 1.csv",sep=",", header = TRUE)
supp_table = supp_table %>% filter(Technology == 'Microarrays')
RFDiseaseList = as.vector(unique(as.vector((supp_table %>% filter(RF.ST == 'RF'))['Disease'])))[[1]]
STDiseaseList = as.vector(unique(as.vector((supp_table %>% filter(RF.ST == 'ST'))['Disease'])))[[1]]

#set first 2 manually.
setwd("/storage/md_shenorr/inbaltz/microarrays/RF/Adenocarcinoma/")
tempdf <- readRDS("comp_MM_Adenocarcinoma_2_HS_Adenocarcinoma_1.rds")
merged <- readRDS("comp_MM_Adenocarcinoma_3_HS_Adenocarcinoma_1.rds")
tempdf['disease']="Adenocarcinoma"
tempdf['comp_file']= "comp_MM_Adenocarcinoma_2_HS_Adenocarcinoma_1.rds"
merged['disease']="Adenocarcinoma"
merged['comp_file']= "comp_MM_Adenocarcinoma_3_HS_Adenocarcinoma_1.rds"
merged <- rbind(merged,tempdf)
RFDiseaseList=RFDiseaseList[RFDiseaseList != "Adenocarcinoma"];

for (disease in RFDiseaseList){
  data_path = paste0("/storage/md_shenorr/inbaltz/microarrays/RF/",disease,"/")
  setwd(data_path)
  comp_files=system2(command = "ls", args = ".", stdout = TRUE)
  comp_files = comp_files[substring(comp_files, 1, 5) =="comp_"]
  for (file in comp_files){
    tempdf = readRDS(file)
    tempdf['disease']=disease
    tempdf['comp_file']= file
    merged <- rbind(merged,tempdf)
  }
}

for (disease in STDiseaseList){
  data_path = paste0("/storage/md_shenorr/inbaltz/microarrays/ST/",disease,"/")
  setwd(data_path)
  comp_files=system2(command = "ls", args = ".", stdout = TRUE)
  comp_files = comp_files[substring(comp_files, 1, 5) =="comp_"]
  for (file in comp_files){
    tempdf = readRDS(file)
    tempdf['disease']=disease
    tempdf['comp_file']= file
    merged <- rbind(merged,tempdf)
  }
}
merged["MM.qval"]=-1
merged["HS.qval"]=-1

###############################################  RNAseq  ###############################################

RFDiseaseList = c("ASPS","DownSyndrome","MyotonicDystrophy","ProstateCancer")
STDiseaseList = c("AtopicDermatitis","Autism","Huntington","MyotonicDystrophy","PulmonaryFibrosis")

for (disease in RFDiseaseList){
  data_path = paste0("/storage/md_shenorr/inbaltz/RNAseq/RF/",disease,"/")
  setwd(data_path)
  comp_files=system2(command = "ls", args = ".", stdout = TRUE)
  comp_files = comp_files[substring(comp_files, 1, 5) =="comp_"]
  for (file in comp_files){
    tempdf = readRDS(file)
    tempdf['disease']=disease
    tempdf['comp_file']= file
    tempdf["MM.FC"]=-1
    tempdf["HS.FC"]=-1
    merged <- rbind(merged,tempdf)
  }
}

for (disease in STDiseaseList){
  data_path = paste0("/storage/md_shenorr/inbaltz/RNAseq/ST/",disease,"/")
  setwd(data_path)
  comp_files=system2(command = "ls", args = ".", stdout = TRUE)
  comp_files = comp_files[substring(comp_files, 1, 5) =="comp_"]
  for (file in comp_files){
    tempdf = readRDS(file)
    tempdf['disease']=disease
    tempdf['comp_file']= file
    tempdf["MM.FC"]=-1
    tempdf["HS.FC"]=-1
    merged <- rbind(merged,tempdf)
  }
}
merged <- merged[! is.na(merged$human_entrez_ID),]  #comp_MM_ASPS_2_HS_ASPS_1.rds has many human_entrez_ID NA values, all paired w/ mouse gene 338320.
saveRDS(merged, file = "/storage/md_shenorr/inbaltz/microarrays/all_comp_files.rds")
