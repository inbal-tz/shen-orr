require(readr)
require(dplyr)
require(plyr)
library(gdata)

data <- readRDS("/storage/md_shenorr/inbaltz/data_no_one2one.rds")
data <- filter(data, ortholog_type != "1:1")
QUEUE = "S"
WAIT_TIME=2
HUMAN_GENES_PER_PROCESS = 40
human_genes = unique(data$human_entrez_ID)
human_genes_for_job=""
for (test_disease in unique(data$diseas)){
  for(i in 1:length(human_genes)){
    human_genes_for_job = paste0(human_genes_for_job,":",human_genes[i])
    if (i%%HUMAN_GENES_PER_PROCESS == 0){
      human_genes_for_job=substr(human_genes_for_job,2,nchar(human_genes_for_job))
      job_name = paste0(test_disease,"_",i/HUMAN_GENES_PER_PROCESS)
      my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/lasso_on_whole_groups/outerr/ -e /storage/md_shenorr/inbaltz/lasso_on_whole_groups/outerr/",' -v human_genes="',human_genes_for_job,'",test_disease="',test_disease,'",file_index="',(i/HUMAN_GENES_PER_PROCESS),'" /storage/md_shenorr/inbaltz/scripts/job_lasso_on_whole_group.R')
      message("Executing: qsub", my_args)
      out=system2(command = "qsub", args = my_args, stdout = TRUE)
      human_genes_for_job=""
      Sys.sleep(WAIT_TIME)
    }
  }
  #send last job manually:
  human_genes_for_job=substr(human_genes_for_job,2,nchar(human_genes_for_job))
  job_name = paste0(test_disease,"_",floor(i/HUMAN_GENES_PER_PROCESS)+1)
  my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/lasso_on_whole_groups/outerr/ -e /storage/md_shenorr/inbaltz/lasso_on_whole_groups/outerr/",' -v human_genes="',human_genes_for_job,'",test_disease="',test_disease,'",file_index="',(floor(i/HUMAN_GENES_PER_PROCESS)+1),'" /storage/md_shenorr/inbaltz/scripts/job_lasso_on_whole_group.R')
  message("Executing: qsub", my_args)
  out=system2(command = "qsub", args = my_args, stdout = TRUE)
  human_genes_for_job=""
  Sys.sleep(WAIT_TIME)
}


#merge all files that were created into 1 df. and save it.
job_name = "merge_all_betas_trained_as_couples"
my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/lasso_on_whole_groups/outerr/ -e /storage/md_shenorr/inbaltz/lasso_on_whole_groups/outerr/ /storage/md_shenorr/inbaltz/scripts/job_merge_betas_trained_as_couples.r")
message("Executing: qsub", my_args)
out=system2(command = "qsub", args = my_args, stdout = TRUE)

job_name = "merge_all_betas_trained_as_groups"
my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/lasso_on_whole_groups/outerr/ -e /storage/md_shenorr/inbaltz/lasso_on_whole_groups/outerr/ /storage/md_shenorr/inbaltz/scripts/job_merge_betas_trained_as_groups.r")
message("Executing: qsub", my_args)
out=system2(command = "qsub", args = my_args, stdout = TRUE)
