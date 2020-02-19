#!/Local/md_shenorr/R-3.2.0/bin/Rscript
library(Matrix)
library(glmnet)
library(plyr)
library(methods)
library(gdata)
require(readr)
require(dplyr)
require(grid)

human_genes = Sys.getenv("human_genes")
i = Sys.getenv("index")
human_genes = strsplit(human_genes,":")[[1]]

data = readRDS("/storage/md_shenorr/inbaltz/data_no_one2one.rds")
coeffs = data.frame()
for (human_gene in human_genes){ #each job is given 20 human genes.
  # cat("gene number", human_gene,'\n')
  for (mouse_gene in unique(filter(data,human_entrez_ID==human_gene)$mouse_entrez_ID)){
    temp <- filter(data, human_entrez_ID==human_gene & mouse_entrez_ID==mouse_gene)
    # cat(nrow(temp),'\n')
    if (nrow(temp)>=30){ #only look at genes that have enough data
      for (b in 1:100){ #bootstraping 100 times
        sampled_rows <- sample(1:nrow(temp), nrow(temp),replace = TRUE) #sampling the same number of samples but with replacement, so same data sample can be sampled more than once.
        design_mat = cbind(rep(1, nrow(temp)), temp$MM.Ztest[sampled_rows]) #the column on 1 is to regulate the beta parameter to 1.
        colnames(design_mat)=c("Intercept","MM.Ztest")
        res = temp$HS.Ztest[sampled_rows] - temp$MM.Ztest[sampled_rows] #why did rachelly do here HSZtest-MMZtest???
        glmmod <- cv.glmnet(design_mat, res, alpha = 1, intercept = FALSE) #inrecept=false is to regulate the intercept to 0
        lambda <- glmmod$lambda.1se
        # predglm <- predict(glmmod, cbind(rep(1, nrow(filter(data_test, human_entrez_ID==human_gene & mouse_entrez_ID==mouse_gene))), filter(data_test, human_entrez_ID==human_gene & mouse_entrez_ID==mouse_gene)$MM.Ztest), lambda)
        curr_coef = unmatrix(as.matrix(coef(glmmod, s = lambda)))
        mse = mean((temp$MM.Ztest*curr_coef["MM.Ztest:1"]- temp$HS.Ztest)^2)
        coeffs_1 = data.frame(human_entrez_ID=human_gene,mouse_entrez_ID=mouse_gene, 
                              Slope= curr_coef["MM.Ztest:1"], 
                              Intercept=curr_coef["(Intercept):1"],
                              bootstrap=b, Lambda = lambda, mse = mse, row.names = NULL)
        coeffs = rbind(coeffs, coeffs_1)
      }
    } 
    
  }
}
saveRDS(coeffs,paste0("/storage/md_shenorr/inbaltz/BootstrapingLasso/bootstrapingResults_",i,".RDS"))
