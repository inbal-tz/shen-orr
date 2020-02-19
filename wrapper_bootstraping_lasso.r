require(readr)
require(dplyr)
require(plyr)
library(gdata)

data <- readRDS("/storage/md_shenorr/inbaltz/all_comp_files.rds")
data <- filter(data, ortholog_type != "1:1")
saveRDS(data, "/storage/md_shenorr/inbaltz/data_no_one2one.rds")
test_comp_files = sample(unique(data$comp_file),round(length(unique(data$comp_file))*0.15))
data_train <- filter(data, !(comp_file %in% test_comp_files))
saveRDS(data_train,"/storage/md_shenorr/inbaltz/trainData.RDS")
data_test <- filter(data, comp_file %in% test_comp_files)
saveRDS(data_test,"/storage/md_shenorr/inbaltz/testData.RDS")
QUEUE = "S"
WAIT_TIME=3

human_genes = unique(data_train$human_entrez_ID)
human_genes_for_job=""
for(i in 1:length(human_genes)){
  human_genes_for_job = paste0(human_genes_for_job,":",human_genes[i])
  if (i%%20 == 0){
    human_genes_for_job=substr(human_genes_for_job,2,nchar(human_genes_for_job))
    job_name = paste0("Bootstraping_",i/20)
    my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/BootstrapingLasso/outerr/ -e /storage/md_shenorr/inbaltz/BootstrapingLasso/outerr/",' -v human_genes="',human_genes_for_job,'",index="',(i/20),'" /storage/md_shenorr/inbaltz/scripts/job_bootstraping_lasso.r')
    message("Executing: qsub", my_args)
    out=system2(command = "qsub", args = my_args, stdout = TRUE)
    human_genes_for_job=""
    Sys.sleep(WAIT_TIME)
  }
}
#send last job manually:
human_genes_for_job=substr(human_genes_for_job,2,nchar(human_genes_for_job))
job_name = paste0("Bootsraping_",floor(i/20))
my_args= paste0(" -V -q " ,QUEUE," -N ",job_name," -o /storage/md_shenorr/inbaltz/BootstrapingLasso/outerr/ -e /storage/md_shenorr/inbaltz/BootstrapingLasso/outerr/",' -v human_genes="',human_genes_for_job,'",index="',(floor(i/20)),'" /storage/md_shenorr/inbaltz/scripts/job_bootstraping_lasso.r')
message("Executing: qsub", my_args)
out=system2(command = "qsub", args = my_args, stdout = TRUE)


#merge all bootstraping files:
bootdtrapped_lasso_results = data.frame()
setwd("/storage/md_shenorr/inbaltz/BootstrapingLasso/")

files=system2(command = "ls", args = ".", stdout = TRUE)
for (file in files){
  bootdtrapped_lasso_results = rbind(bootdtrapped_lasso_results,readRDS(file))
}
saveRDS(bootdtrapped_lasso_results,"/storage/md_shenorr/inbaltz/BootstrapingLasso/bootdtrapped_lasso_all_results.RDS")


