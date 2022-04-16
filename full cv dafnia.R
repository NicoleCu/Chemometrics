library(readxl)
library(mdatools)
library(Metrics)
library(e1071)
library(KRLS)
library(randomForest)


DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",
                                     sheet = "analysis_dafnia", range ='A1:R51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

# remove outliers 
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% 
                                c('4s','4b','8b','8s','28s','28b',"29s","29b",'30b','30s','24b')),]
prop <- round(nrow(DafniaMagna)*0.70) 
vec_PLS_full <- vector()
vec_RF_full <- vector()
vec_KRLS_full <- vector()
vec_SVM_full <- vector()


for (i in 1:nrow(DafniaMagna)) {
  m <- pls(DafniaMagna[-i,1:16], DafniaMagna[-i,17], 8, cv = 1)
  r <- predict(m, data.frame(DafniaMagna[i,1:16]))
  RMSEP_PLS_full <- rmse(r$y.pred[, r[["ncomp.selected"]],1], DafniaMagna[i,17])
  vec_PLS_full[i] <- RMSEP_PLS_full
  
  OptModelsvm <- tune(svm, Daphnia ~ ., data = DafniaMagna[-i,], kernel = 'radial',
                      ranges=list(cost = c(1,10,100),  
                                  epsilon = seq(0.01, 0.15, 0.01),
                                  gamma = c(0.0001, 0.0005, 0.0001)), 
                      tunecontrol = tune.control(cross = prop), scale = F)
  BstModel <- OptModelsvm$best.model
  pred <- predict(BstModel, DafniaMagna[i,])
  RMSEP_SVM_full <- rmse(DafniaMagna[i,17], pred)
  vec_SVM_full[i] <- RMSEP_SVM_full
  
  RF <- randomForest(Daphnia ~., data = DafniaMagna[-i,], xtest = DafniaMagna[i,1:16], ytest = DafniaMagna[i,17],  
                     mtry = 5, ntree = 500, nodesize = 5, importance = T)
  RMSEP_RF_full <- rmse(DafniaMagna[i,17], RF[["test"]][["predicted"]])
  vec_RF_full[i] <-  RMSEP_RF_full
  
  KRLS <- krls(DafniaMagna[-i,1:16], as.vector(DafniaMagna[-i,17]), vcov = T, derivative = T)
  KRLS_pred <- predict(KRLS,newdata = DafniaMagna[i,1:16], se.fit=TRUE)
  RMSEP_KRLS_full <- rmse(DafniaMagna[i,17], KRLS_pred$fit)
  vec_KRLS_full[i] <- RMSEP_KRLS_full
}


av_PLS_error_fullcv <- sum(vec_PLS_full)/length(vec_PLS_full)
av_SVM_error_fullcv <- sum(vec_SVM_full)/length(vec_SVM_full)
av_RF_error_fullcv <- sum(vec_RF_full)/length(vec_RF_full)
av_KRLS_error_fullcv <- sum(vec_KRLS_full)/length(vec_KRLS_full)
names(vec_KRLS_full) <- rownames(DafniaMagna)
names(vec_PLS_full) <- rownames(DafniaMagna)
names(vec_SVM_full) <- rownames(DafniaMagna)
names(vec_RF_full) <- rownames(DafniaMagna)

