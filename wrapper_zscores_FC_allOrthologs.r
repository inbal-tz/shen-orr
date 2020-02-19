#### wrapper to send jobs for ATLAS. 
#get the disease lists for the RF and ST microarray data used from the table in the paper supplamenty
#in each disease folder get the files from the subDatasets file in each folder (created by rachely)
require(dplyr)
supp_table = read.table("/storage/md_shenorr/inbaltz/Supplementary Table 1.csv",sep=",", header = TRUE)
supp_table = supp_table %>% filter(Technology == 'Microarrays')
RFDiseaseList = as.vector(unique(as.vector((supp_table %>% filter(RF.ST == 'RF'))['Disease'])))[[1]]
STDiseaseList = as.vector(unique(as.vector((supp_table %>% filter(RF.ST == 'ST'))['Disease'])))[[1]]

for (disease in RFDiseaseList){
  data_path = paste0("/storage/md_shenorr/rachelly/Shared_Dir/GS_analysis/Microarrays/RF/",disease,"/V1.9/")
  setwd(data_path)
  prefixes = read.table("subDatasets", nrows = 2, sep="\t", stringsAsFactors = FALSE)
  fsub_datasets = read.table(file.path("subDatasets"),skip = 2, sep="\t", stringsAsFactors = FALSE)
  for( job.id in 1: nrow(fsub_datasets)){
    species = fsub_datasets[job.id,4]
    if (species == "Homo sapiens") prefix= prefixes[1,1] else prefix = prefixes[2,1]
    RDSfile = paste0(prefix,"_", job.id,"_data.rds")
    if (!file.exists(RDSfile)) stop("Error (job_2_DE_genes_PWs): the following file doesn't exist: ", RDSfile)
    message(RDSfile)
    QUEUE = "S"
    WAIT_TIME=3
    job_name = paste0(disease,"_",job.id) #job name includes the disease name and serial number
    my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/microarrays/outerr/ -e /storage/md_shenorr/inbaltz/microarrays/outerr/",' -v RDSfile="', RDSfile,'",disease="',disease,'",RFST="RF" /storage/md_shenorr/inbaltz/scripts/job_zscores_FC_allOrthologs.r')
    message("Executing: qsub", my_args)
    out=system2(command = "qsub", args = my_args, stdout = TRUE)
    Sys.sleep(WAIT_TIME)
  }
}
 
for (disease in STDiseaseList){
  data_path = paste0("/storage/md_shenorr/rachelly/Shared_Dir/GS_analysis/Microarrays/ST/",disease,"/V1.9/")
  if (disease == "Sepsis") data_path = "/storage/md_shenorr/rachelly/Shared_Dir/GS_analysis/Microarrays/ST/Sepsis/V1.9_old/"
  setwd(data_path)
  prefixes = read.table("subDatasets", nrows = 2, sep="\t", stringsAsFactors = FALSE)
  fsub_datasets = read.table(file.path("subDatasets"),skip = 2, sep="\t", stringsAsFactors = FALSE)
  for( job.id in 1: nrow(fsub_datasets)){
    species = fsub_datasets[job.id,4]
    if (species == "Homo sapiens") prefix= prefixes[1,1] else prefix = prefixes[2,1]
    RDSfile = paste0(prefix,"_", job.id,"_data.rds")
    if (!file.exists(RDSfile)) stop("Error (job_2_DE_genes_PWs): the following file doesn't exist: ", RDSfile)
    message(RDSfile)
    #todo maybe check est here???
    QUEUE = "S"
    WAIT_TIME=3
    job_name = paste0(disease,"_",job.id) #job name includes the disease name and serial number
    my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/microarrays/outerr/ -e /storage/md_shenorr/inbaltz/microarrays/outerr/",' -v RDSfile="', RDSfile,'",disease="',disease,'",RFST="ST" /storage/md_shenorr/inbaltz/scripts/job_zscores_FC_allOrthologs.r')
    message("Executing: qsub", my_args)
    out=system2(command = "qsub", args = my_args, stdout = TRUE)
    Sys.sleep(WAIT_TIME)
  }
}

  
  