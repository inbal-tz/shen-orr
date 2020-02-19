#!/Local/md_shenorr/R-3.2.0/bin/Rscript
library(Matrix)
library(glmnet)
set.seed(1)
library(plyr)
library(methods)
library(gdata)
require(readr)
require(dplyr)
require(grid)
require(tidyr)


test_disease = Sys.getenv("test_disease")
human_genes= Sys.getenv("human_genes")
human_genes = strsplit(human_genes,":")[[1]]
file_index = Sys.getenv("file_index")
# test_disease = "IBD"
# human_genes = c(12,9)
data <- readRDS("/storage/md_shenorr/inbaltz/data_no_one2one.rds")
# data <- readRDS("/storage/md_shenorr/inbaltz/all_comp_files.rds")
# data <- readRDS("C:/Users/inbal/Google Drive/shen-orr/all_comp_files.rds")
data <- filter(data, ortholog_type != "1:1")

data_train = filter(data, disease != test_disease)
data_test = filter(data, disease == test_disease)
weighted_couple_slopes = data.frame()
betas_whole_group = data.frame()
lambda_seq <- 10^seq(2, -2, by = -.1)

for (human_gene in human_genes){
  cat(human_gene,"\n")
  for (mouse_gene in unique(filter(data_train,human_entrez_ID==human_gene)$mouse_entrez_ID)){
    temp <- filter(data_train, human_entrez_ID==human_gene & mouse_entrez_ID==mouse_gene)
    if (nrow(temp)>=30){ #only look at genes that have enough data
      for (b in 1:100){ #bootstraping 100 times
        sampled_rows <- sample(1:nrow(temp), nrow(temp),replace = TRUE) #sampling the same number of samples but with replacement, so same data sample can be sampled more than once.
        
        #train in 2 ways: A. calculate individual slopes, and the weight them by mse.
        design_mat = cbind(rep(1, nrow(temp)), temp$MM.Ztest[sampled_rows]) #the column on 1 is to regulate the beta parameter to 1.
        colnames(design_mat)=c("Intercept","MM.Ztest")
        res = temp$HS.Ztest[sampled_rows] - temp$MM.Ztest[sampled_rows] 
        glmmod <- cv.glmnet(design_mat, res, alpha = 1, intercept = FALSE, lambda = lambda_seq) #inrecept=false is to regulate the intercept to 0
        lambda <- glmmod$lambda.1se
        
        curr_coef = unmatrix(as.matrix(coef(glmmod, s = lambda)))
        mse = mean((temp$MM.Ztest*curr_coef["MM.Ztest:1"]- temp$HS.Ztest)^2)
        
        coeffs_1 = data.frame(human_entrez_ID=human_gene,mouse_entrez_ID=mouse_gene, 
                              Slope= curr_coef["MM.Ztest:1"], 
                              Intercept=curr_coef["(Intercept):1"],
                              bootstrap=b, Lambda = lambda, mse = mse, row.names = NULL)
        #save all output from the training.
        weighted_couple_slopes = rbind(weighted_couple_slopes, coeffs_1)
      }
      

    }
  }
  
  ##train in 2 ways: B. use lasso regeression on whole group.
  mouse_genes = unique(filter(data_train,human_entrez_ID==human_gene)$mouse_entrez_ID)
  if(length(mouse_genes)<=1) next
  comp_files = unique(filter(data_train,human_entrez_ID==human_gene)$comp_file)
  data_for_laso <- data.frame(matrix(ncol = length(mouse_genes)+1, nrow = length(comp_files)), row.names = comp_files)
  colnames(data_for_laso) <- c(mouse_genes,"HS.Ztest")
  
  human_gene_data=filter(data_train, human_entrez_ID==human_gene) #all of the data relevant to this current human gene
  for (row in 1:nrow(human_gene_data)){
    comp = human_gene_data[row,"comp_file"]
    mouse_gene = human_gene_data[row,"mouse_entrez_ID"]
    data_for_laso[comp,toString(mouse_gene)] = human_gene_data[row,"MM.Ztest"]
    if(is.na(data_for_laso[comp,"HS.Ztest"])) data_for_laso[comp,"HS.Ztest"]=human_gene_data[row,"HS.Ztest"]
  }
  
  data_for_laso <- data_for_laso[,which(colMeans(!is.na(data_for_laso))*nrow(data_for_laso)>=100), drop=FALSE]
  if(ncol(data_for_laso)<=2) next
  data_for_laso <- drop_na(data_for_laso) # only use samples that have all the chosen orthologs! is it right?
  
  # lasso regression.
  # the data: x is all of the orthologs that had >=100 datapoints! y is the HS.Ztest.
  
  normalized_mse = filter(weighted_couple_slopes,human_entrez_ID==human_gene & mouse_entrez_ID %in% names(data_for_laso))
  normalized_mse = normalized_mse[order(normalized_mse$mouse_entrez_ID),]
  normalized_mse = aggregate(normalized_mse, by = list(normalized_mse$mouse_entrez_ID), mean)$mse
  normalized_mse = 1/normalized_mse
  normalized_mse = normalized_mse/sum(normalized_mse)
  x_var <- data_for_laso[,-ncol(data_for_laso)]
  x_var <- x_var[,order(names(x_var))] #sort columns by gene#, so same order will be used at prediction time
  for (b in 1:100){
    # betas_whole_group_b = data.frame()
    sampled_rows <- sample(1:nrow(x_var), nrow(x_var),replace = TRUE)
    x_var_b <- x_var[sampled_rows,]
    y_var_b <- (data_for_laso$HS.Ztest)[sampled_rows]
    # y_var_b <- (data_for_laso$HS.Ztest)[sampled_rows] + as.matrix(x_var_b) %*% normalized_mse #matrix multaplication to implement the equation i want..
    
    glmmod <- cv.glmnet(as.matrix(x_var_b), y_var_b, alpha = 1, intercept = FALSE) #this is shrinked to 0 not 1
    lambda <- glmmod$lambda.min
    
    curr_coef = unmatrix(as.matrix(coef(glmmod, s = lambda)))
    curr_coef = curr_coef[2:length(curr_coef)]
    
    
    temp = data.frame()
    for (i in 1:length(curr_coef)){
      temp = rbind(temp, data.frame(human_entrez_ID = human_gene,
                                    mouse_entrez_ID = strsplit(names(curr_coef)[i], ":")[[1]][1],
                                    beta = as.vector(curr_coef)[i]
      ))
    }
    temp$lambda = lambda
    temp$b = b
    betas_whole_group = rbind(betas_whole_group, temp)
  }
}
if (nrow(betas_whole_group)>0){
  betas_whole_group$test_disease = test_disease
}
saveRDS(weighted_couple_slopes, paste0("/storage/md_shenorr/inbaltz/lasso_on_whole_groups/trained_as_couples/betasCouples_testDisease_", test_disease,"_", file_index, ".RDS"))
saveRDS(betas_whole_group, paste0("/storage/md_shenorr/inbaltz/lasso_on_whole_groups/trained_as_groups/betasGroups_testDisease_", test_disease,"_", file_index, ".RDS"))
