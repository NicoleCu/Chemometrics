library(readxl)
library(mdatools)
library(Metrics)
library(e1071)
library(KRLS)
library(randomForest)

VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx", 
                                       sheet = "analysis_fisheri", range ='A1:Y27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]

# remove outliers 
VibrioFisheri <-  VibrioFisheri[!(row.names(VibrioFisheri) %in% c('w24', 'w6')),]

prop <- round(nrow(VibrioFisheri)*0.70) 
vec_PLS_full <- vector()
vec_RF_full <- vector()
vec_KRLS_full <- vector()
vec_SVM_full <- vector()



for (i in 1:nrow(VibrioFisheri)) {
  m <- pls(VibrioFisheri[-i,1:23], VibrioFisheri[-i,24], 8, cv = 1)
  r <- predict(m, data.frame(VibrioFisheri[i,1:23]))
  RMSEP_PLS_full <- rmse(r$y.pred[, r[["ncomp.selected"]],1], VibrioFisheri[i,24])
  vec_PLS_full[i] <- RMSEP_PLS_full
  
  OptModelsvm <- tune(svm, EC50 ~ ., data = VibrioFisheri[-i,], kernel = 'radial',
                      ranges=list(cost = c(1,10,100),  
                                  epsilon = seq(0.01, 0.15, 0.01),
                                  gamma = c(0.0001, 0.0005, 0.0001)), 
                      tunecontrol = tune.control(cross = prop), scale = F)
  BstModel <- OptModelsvm$best.model
  pred <- predict(BstModel, VibrioFisheri[i,])
  RMSEP_SVM_full <- rmse(VibrioFisheri[i,24], pred)
  vec_SVM_full[i] <- RMSEP_SVM_full
  
  RF <- randomForest(EC50 ~., data = VibrioFisheri[-i,], xtest = VibrioFisheri[i,1:23], ytest = VibrioFisheri[i,24],  
                     mtry = 5, ntree = 500, nodesize = 5, importance = T)
  RMSEP_RF_full <- rmse(VibrioFisheri[i,24], RF[["test"]][["predicted"]])
  vec_RF_full[i] <-  RMSEP_RF_full
  
  KRLS <- krls(VibrioFisheri[-i,1:23], as.vector(VibrioFisheri[-i,24]), vcov = T, derivative = T)
  KRLS_pred <- predict(KRLS,newdata = VibrioFisheri[i,1:23], se.fit=TRUE)
  RMSEP_KRLS_full <- rmse(VibrioFisheri[i,24], KRLS_pred$fit)
  vec_KRLS_full[i] <- RMSEP_KRLS_full
}


av_PLS_error_fullcv <- sum(vec_PLS_full)/length(vec_PLS_full)
av_SVM_error_fullcv <- sum(vec_SVM_full)/length(vec_SVM_full)
av_RF_error_fullcv <- sum(vec_RF_full)/length(vec_RF_full)
av_KRLS_error_fullcv <- sum(vec_KRLS_full)/length(vec_KRLS_full)
