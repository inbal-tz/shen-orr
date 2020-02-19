#trained as couples
library(plyr)
require(readr)
require(dplyr)
require(ggplot2)
library(ggpubr)
data <- readRDS("C:/Users/inbal/Google Drive/shen-orr/all_comp_files_no_one2one.RDS")
betas_trained_as_couples <- readRDS("C:/Users/inbal/Google Drive/shen-orr/lasso_regressions/betas_trained_as_couples.RDS")
betas_trained_as_couples <- aggregate(list(slope=betas_trained_as_couples$Slope, lambda=betas_trained_as_couples$Lambda, mse=betas_trained_as_couples$mse), 
                                      by = list(human_entrez_ID=betas_trained_as_couples$human_entrez_ID, mouse_entrez_ID=betas_trained_as_couples$mouse_entrez_ID, test_disease=betas_trained_as_couples$test_disease), 
                                      mean)
betas_trained_as_groups <- readRDS("C:/Users/inbal/Google Drive/shen-orr/lasso_regressions/betas_trained_as_groups.RDS")
betas_trained_as_groups <- aggregate(list(beta=betas_trained_as_groups$beta, lambda=betas_trained_as_groups$lambda), 
                                     by = list(human_entrez_ID=betas_trained_as_groups$human_entrez_ID, mouse_entrez_ID=betas_trained_as_groups$mouse_entrez_ID, test_disease=betas_trained_as_groups$test_disease), 
                                     mean)

predictions_trained_as_couples = data.frame()
predictions_trained_as_groups = data.frame()
# testD="IBD"
for(testD in unique(data$disease)){
betas_test_couples <- filter(betas_trained_as_couples, test_disease==testD)
betas_test_groups <- filter(betas_trained_as_groups, test_disease==testD)
data_test <- filter(data, disease==testD)

#we have the training information 2 ways: A. training as couples: weighted_couple_slopes.
#                                         B. training as a whole group: betas_whole_group.
#test it!
for (human_gene in unique(data_test$human_entrez_ID)){
  for (comp in unique(data_test$comp_file)){
    temp <- filter(data_test, comp_file==comp & human_entrez_ID==human_gene)
    temp <- merge(temp, betas_test_couples, by = c("human_entrez_ID","mouse_entrez_ID"))
    # temp$normalized_mse = temp$mse/sum(temp$mse)
    temp = temp[order(temp$mouse_entrez_ID),]
    if( nrow(temp) > 0){
      #calculate in 2 ways:
      
      #A: using slopes from individualized couples, weighed inversly to their mse.
      pred = sum((temp$MM.Ztest * temp$slope) * (1/temp$mse))/sum(1/temp$mse)
      predictions_trained_as_couples = rbind(predictions_trained_as_couples, 
                                             data.frame(human_entrez_ID = human_gene,
                                                        mouse_entrez_ID = I(list(temp$mouse_entrez_ID)),
                                                        prediction = pred,
                                                        HS.Ztest = temp$HS.Ztest[1],
                                                        test_disease = testD,
                                                        comp_file = comp))
    }
    #B: using lasso regsression to calculate
    if(length(temp$mouse_entrez_ID) > 0 && length(temp$mouse_entrez_ID) == length(as.vector(filter(betas_test_groups, human_entrez_ID==human_gene)$mouse_entrez_ID))){
      if( all(sort(temp$mouse_entrez_ID) == as.vector(filter(betas_test_groups, human_entrez_ID==human_gene)$mouse_entrez_ID))){
        pred = sum(filter(betas_test_groups, human_entrez_ID==human_gene)$beta * temp$MM.Ztest) #- sum((temp$MM.Ztest)*((1/(temp$mse))/sum(1/(temp$mse))))
        predictions_trained_as_groups = rbind(predictions_trained_as_groups, 
                                              data.frame(human_entrez_ID = human_gene,
                                                         prediction = pred,
                                                         HS.Ztest = temp$HS.Ztest[1],
                                                         test_disease = testD,
                                                         comp_file = comp))
      }
    }
    
  }
}
prediction_couples_graph <- filter(predictions_trained_as_couples, test_disease == testD)
prediction_groups_graph <- filter(predictions_trained_as_groups, test_disease == testD)

cor = c(cor(prediction_couples_graph$HS.Ztest, prediction_couples_graph$prediction, method="spearman"),
                cor(prediction_groups_graph$HS.Ztest, prediction_groups_graph$prediction, method="spearman"),
                cor(data_test$HS.Ztest, data_test$MM.Ztest, method="spearman"))
mse = c(mean((prediction_couples_graph$prediction-prediction_couples_graph$HS.Ztest)^2),
        mean((prediction_groups_graph$prediction-prediction_groups_graph$HS.Ztest)^2),
        mean((data_test$HS.Ztest-data_test$MM.Ztest)^2))
g1 <- ggplot(data.frame(mse, category=c("trained as\ncouples", "trained as\ngroups", "original\ndata")), aes(x=category, y=mse))+geom_bar(stat="identity")+xlab("")
g2 <- ggplot(data.frame(cor, category=c("trained as\ncouples", "trained as\ngroups", "original\ndata")),aes(x=category, y=cor))+geom_bar(stat="identity")+xlab("")

cat("testD: ", testD, "\nmse: ", mse, "\ncor: ", cor, "\n")

couplePlot <- ggplot(prediction_couples_graph, aes(x=HS.Ztest, y=prediction)) + 
  geom_point()+
  ggtitle("\nTrained As Couples" ,subtitle=paste0("#human genes predicted: ", length(unique(prediction_couples_graph$human_entrez_ID))))

groupPlot <- ggplot(prediction_groups_graph, aes(x=HS.Ztest, y=prediction)) + 
  geom_point()+
  ggtitle("\nTrained As Groups", subtitle=paste0("#human genes predicted: ", length(unique(prediction_groups_graph$human_entrez_ID))))

groundTruth <- ggplot(data_test, aes(x=HS.Ztest, y=MM.Ztest))+geom_point()+
  ggtitle("original data")

correlationsplot <- ggarrange(g1, g2, nrow=2)
multiplot <- ggarrange(ggarrange(couplePlot, groupPlot, ncol=2), ggarrange(groundTruth, correlationsplot, ncol=2), nrow=2, labels = paste0("Test Disease: ", testD) )
ggsave(paste0("C:/Users/inbal/Google Drive/shen-orr/lasso_regressions/plots/",testD,".jpg"))
}




#calc error per test disease
predictions_trained_as_couples$squared_error = (predictions_trained_as_couples$prediction-predictions_trained_as_couples$HS.Ztest)^2
aggregate(list(mse=predictions_trained_as_couples$squared_error), by = list(test_disease = predictions_trained_as_couples$test_disease), mean)

predictions_trained_as_groups$squared_error = (predictions_trained_as_groups$prediction-predictions_trained_as_groups$HS.Ztest)^2
aggregate(list(mse=predictions_trained_as_groups$squared_error), by = list(test_disease = predictions_trained_as_groups$test_disease), mean)

