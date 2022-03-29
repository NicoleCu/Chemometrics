library(readxl)
library(mdatools)
library(Metrics)
library(e1071)
library(KRLS)


DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                                     range ='B2:N23')) 
rownames(DafniaMagna) <- seq(0,20, 1)
DafniaMagna <- DafniaMagna[,-12]

# remove outliers 
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% c(11,19, 9, 18)),]
prop <- round(nrow(DafniaMagna)*0.70) 
vec_PLS <- vector()
vec_RF <- vector()
vec_KRLS <- vector()
vec_SVM <- vector()


for (i in 1:100) {
  s <- sample(1:nrow(DafniaMagna), prop)
  m <- pls(DafniaMagna[s,1:11], DafniaMagna[s,12], 8, cv = 1)
  r <- predict(m, DafniaMagna[-s,1:11], DafniaMagna[-s,12])
  RMSEP_PLS <- rmse(r$y.pred[, r[["ncomp.selected"]],1], r$y.ref[,1])
  vec_PLS[i] <- RMSEP_PLS
  
  OptModelsvm <- tune(svm,  Daphnia..survival.rate ~ ., data = DafniaMagna[s,], kernel = 'radial',
                      ranges=list(cost = c(1,10,100),  
                                  epsilon = seq(0.01, 0.15, 0.01),
                                  gamma = c(0.0001, 0.0005, 0.0001)), 
                      tunecontrol = tune.control(cross = prop), scale = F)
  BstModel <- OptModelsvm$best.model
  pred <- predict(BstModel, DafniaMagna[-s,])
  RMSEP_SVM <- rmse(DafniaMagna[-s,12], pred)
  vec_SVM[i] <- RMSEP_SVM
  
  RF <- randomForest(Daphnia..survival.rate ~., data = DafniaMagna[s,], xtest = DafniaMagna[-s,1:11], ytest = DafniaMagna[-s,12],  
                     mtry = 5, ntree = 500, nodesize = 5, importance = T)
  RMSEP_RF <- rmse(DafniaMagna[-s,12], RF[["test"]][["predicted"]])
  vec_RF[i] <-  RMSEP_RF
  
  KRLS <- krls(DafniaMagna[s,1:11], as.vector(DafniaMagna[s,12]), vcov = T, derivative = T)
  KRLS_pred <- predict(KRLS,newdata = DafniaMagna[-s,1:11], se.fit=TRUE)
  RMSEP_KRLS <- rmse(DafniaMagna[-s,12], KRLS_pred$fit)
  vec_KRLS[i] <- RMSEP_KRLS
}


av_PLS_error <- sum(vec_PLS)/length(vec_PLS)
av_SVM_error <- sum(vec_SVM)/length(vec_SVM)
av_RF_error <- sum(vec_RF)/length(vec_RF)
av_KRLS_error <- sum(vec_KRLS)/length(vec_KRLS)
