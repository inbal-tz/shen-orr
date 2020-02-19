#### wrapper to send jobs for ATLAS. 
#get the disease lists for the RF and ST microarray data used from the table in the paper supplamenty
#in each disease folder get the files from the subDatasets file in each folder (created by rachely)
require(dplyr)

###############################################  microarrays  ###############################################

supp_table = read.table("/storage/md_shenorr/inbaltz/Supplementary Table 1.csv",sep=",", header = TRUE)
supp_table = supp_table %>% filter(Technology == 'Microarrays')
RFDiseaseList = as.vector(unique(as.vector((supp_table %>% filter(RF.ST == 'RF'))['Disease'])))[[1]]
STDiseaseList = as.vector(unique(as.vector((supp_table %>% filter(RF.ST == 'ST'))['Disease'])))[[1]]
QUEUE = "S"
WAIT_TIME=3

for (disease in RFDiseaseList){
  data_path = paste0("/storage/md_shenorr/inbaltz/microarrays/RF/",disease,"/")
  setwd(data_path)
  job_name = paste0("RF_",disease) #job name includes the disease name and serial number
  my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/microarrays/outerr/ -e /storage/md_shenorr/inbaltz/microarrays/outerr/",' -v disease="',disease,'",RFST="RF" /storage/md_shenorr/inbaltz/scripts/job_calc_comp_files_panther.r')
  message("Executing: qsub", my_args)
  out=system2(command = "qsub", args = my_args, stdout = TRUE)
  Sys.sleep(WAIT_TIME)
}

for (disease in STDiseaseList){
  data_path = paste0("/storage/md_shenorr/inbaltz/microarrays/ST/",disease,"/")
  setwd(data_path)
  job_name = paste0("ST_",disease) #job name includes the disease name and serial number
  my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/microarrays/outerr/ -e /storage/md_shenorr/inbaltz/microarrays/outerr/",' -v disease="',disease,'",RFST="ST" /storage/md_shenorr/inbaltz/scripts/job_calc_comp_files_panther.r')
  message("Executing: qsub", my_args)
  out=system2(command = "qsub", args = my_args, stdout = TRUE)
  Sys.sleep(WAIT_TIME)
}

###############################################  RNAseq  ###############################################

#some of the slueth files in the ST category are missing. maybe ask rachelly how to fix it. 
RFDiseaseList = c("ASPS","DownSyndrome","MyotonicDystrophy","ProstateCancer")
STDiseaseList = c("AtopicDermatitis","Autism","Huntington","MyotonicDystrophy","PulmonaryFibrosis")

QUEUE = "S"
WAIT_TIME=3

for (disease in RFDiseaseList){
  data_path = paste0("/storage/md_shenorr/inbaltz/RNAseq/RF/",disease,"/")
  setwd(data_path)
  job_name = paste0("RF_",disease) #job name includes the disease name
  my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/RNAseq/outerr/ -e /storage/md_shenorr/inbaltz/RNAseq/outerr/",' -v disease="',disease,'",RFST="RF" /storage/md_shenorr/inbaltz/scripts/job_calc_compfiles_RNAseq_panther.r')
  message("Executing: qsub", my_args)
  out=system2(command = "qsub", args = my_args, stdout = TRUE)
  Sys.sleep(WAIT_TIME)
}

for (disease in STDiseaseList){
  data_path = paste0("/storage/md_shenorr/inbaltz/RNAseq/ST/",disease,"/")
  setwd(data_path)
  job_name = paste0("ST_",disease) #job name includes the disease name
  my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/RNAseq/outerr/ -e /storage/md_shenorr/inbaltz/RNAseq/outerr/",' -v disease="',disease,'",RFST="ST" /storage/md_shenorr/inbaltz/scripts/job_calc_compfiles_RNAseq_panther.r')
  message("Executing: qsub", my_args)
  out=system2(command = "qsub", args = my_args, stdout = TRUE)
  Sys.sleep(WAIT_TIME)
}
