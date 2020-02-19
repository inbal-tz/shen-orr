#!/Local/md_shenorr/R-3.2.0/bin/Rscript
require(readr)
require(dplyr)
require(plyr)

betas_trained_as_couples = data.frame()

path = "/storage/md_shenorr/inbaltz/lasso_on_whole_groups/trained_as_couples"
setwd(path)
file_list = out=system2(command = "ls", args = path, stdout = TRUE)
for (i in 1:length(file_list)){
  test_disease = strsplit(file_list[i], "_")[[1]][3]
  temp <- readRDS(file_list[i])
  if (nrow(temp)>0){
    temp$test_disease = test_disease
    betas_trained_as_couples = rbind(betas_trained_as_couples, temp)
  }
}

saveRDS(betas_trained_as_couples, "/storage/md_shenorr/inbaltz/lasso_on_whole_groups/betas_trained_as_couples.RDS")