library(readxl)
library(mdatools)
library(Metrics)
library(e1071)
library(KRLS)
library(randomForest)


DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                                     range ='B2:N23')) 
rownames(DafniaMagna) <- seq(0,20, 1)
DafniaMagna <- DafniaMagna[,-12]

# remove outliers 
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% c(11, 19, 18, 9, 0)),]
prop <- round(nrow(DafniaMagna)*0.70) 
vec_PLS_full <- vector()
vec_RF_full <- vector()
vec_KRLS_full <- vector()
vec_SVM_full <- vector()


for (i in 1:nrow(DafniaMagna)) {
  m <- pls(DafniaMagna[-i,1:11], DafniaMagna[-i,12], 8, cv = 1)
  r <- predict(m, data.frame(DafniaMagna[i,1:11]))
  RMSEP_PLS_full <- rmse(r$y.pred[, r[["ncomp.selected"]],1], DafniaMagna[i,12])
  vec_PLS_full[i] <- RMSEP_PLS_full
  
  OptModelsvm <- tune(svm, Daphnia..survival.rate ~ ., data = DafniaMagna[-i,], kernel = 'radial',
                      ranges=list(cost = c(1,10,100),  
                                  epsilon = seq(0.01, 0.15, 0.01),
                                  gamma = c(0.0001, 0.0005, 0.0001)), 
                      tunecontrol = tune.control(cross = prop), scale = F)
  BstModel <- OptModelsvm$best.model
  pred <- predict(BstModel, DafniaMagna[i,])
  RMSEP_SVM_full <- rmse(DafniaMagna[i,12], pred)
  vec_SVM_full[i] <- RMSEP_SVM_full
  
  RF <- randomForest(Daphnia..survival.rate ~., data = DafniaMagna[-i,], xtest = DafniaMagna[i,1:11], ytest = DafniaMagna[i,12],  
                     mtry = 5, ntree = 500, nodesize = 5, importance = T)
  RMSEP_RF_full <- rmse(DafniaMagna[i,12], RF[["test"]][["predicted"]])
  vec_RF_full[i] <-  RMSEP_RF_full
  
  KRLS <- krls(DafniaMagna[-i,1:11], as.vector(DafniaMagna[-i,12]), vcov = T, derivative = T)
  KRLS_pred <- predict(KRLS,newdata = DafniaMagna[i,1:11], se.fit=TRUE)
  RMSEP_KRLS_full <- rmse(DafniaMagna[i,12], KRLS_pred$fit)
  vec_KRLS_full[i] <- RMSEP_KRLS_full
}


av_PLS_error_fullcv <- sum(vec_PLS_full)/length(vec_PLS_full)
av_SVM_error_fullcv <- sum(vec_SVM_full)/length(vec_SVM_full)
av_RF_error_fullcv <- sum(vec_RF_full)/length(vec_RF_full)
av_KRLS_error_fullcv <- sum(vec_KRLS_full)/length(vec_KRLS_full)
