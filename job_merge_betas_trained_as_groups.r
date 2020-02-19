#!/Local/md_shenorr/R-3.2.0/bin/Rscript
require(readr)
require(dplyr)
require(plyr)

betas_trained_as_groups = data.frame()

path = "/storage/md_shenorr/inbaltz/lasso_on_whole_groups/trained_as_groups"
setwd(path)
file_list = out=system2(command = "ls", args = path, stdout = TRUE)
for (i in 1:length(file_list)){
  temp <- readRDS(file_list[i])
  if (nrow(temp)>0){
    betas_trained_as_groups = rbind(betas_trained_as_groups, temp)
  }
}

saveRDS(betas_trained_as_groups, "/storage/md_shenorr/inbaltz/lasso_on_whole_groups/betas_trained_as_groups.RDS")