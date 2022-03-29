library(readxl)
library(mdatools)
library(Metrics)
library(e1071)
library(KRLS)


DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",
                                     sheet = "analysis_dafnia", range ='A1:R51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

# remove outliers 
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% 
          c('4s','4b','8b','8s','28s','28b',"29s","29b",'30b','30s','24b')),]
prop <- round(nrow(DafniaMagna)*0.70) 
vec_PLS <- vector()
vec_RF <- vector()
vec_KRLS <- vector()
vec_SVM <- vector()


for (i in 1:100) {
  s <- sample(1:nrow(DafniaMagna), prop)
  m <- pls(DafniaMagna[s,1:16], DafniaMagna[s,17], 8, cv = 1)
  r <- predict(m, DafniaMagna[-s,1:16], DafniaMagna[-s,17])
  RMSEP_PLS <- rmse(r$y.pred[, r[["ncomp.selected"]],1], r$y.ref[,1])
  vec_PLS[i] <- RMSEP_PLS
  
  OptModelsvm <- tune(svm, Daphnia ~ ., data = DafniaMagna[s,], kernel = 'radial',
                      ranges=list(cost = c(1,10,100),  
                      epsilon = seq(0.01, 0.15, 0.01),
                      gamma = c(0.0001, 0.0005, 0.0001)), 
                      tunecontrol = tune.control(cross = prop), scale = F)
  BstModel <- OptModelsvm$best.model
  pred <- predict(BstModel, DafniaMagna[-s,])
  RMSEP_SVM <- rmse(DafniaMagna[-s,17], pred)
  vec_SVM[i] <- RMSEP_SVM
  
  RF <- randomForest(Daphnia ~., data = DafniaMagna[s,], xtest = DafniaMagna[-s,1:16], ytest = DafniaMagna[-s,17],  
                     mtry = 5, ntree = 500, nodesize = 5, importance = T)
  RMSEP_RF <- rmse(DafniaMagna[-s,17], RF[["test"]][["predicted"]])
  vec_RF[i] <-  RMSEP_RF
  
  KRLS <- krls(DafniaMagna[s,1:16], as.vector(DafniaMagna[s,17]), vcov = T, derivative = T)
  KRLS_pred <- predict(KRLS,newdata = DafniaMagna[-s,1:16], se.fit=TRUE)
  RMSEP_KRLS <- rmse(DafniaMagna[-s,17], KRLS_pred$fit)
  vec_KRLS[i] <- RMSEP_KRLS
}


av_PLS_error <- sum(vec_PLS)/length(vec_PLS)
av_SVM_error <- sum(vec_SVM)/length(vec_SVM)
av_RF_error <- sum(vec_RF)/length(vec_RF)
av_KRLS_error <- sum(vec_RMSEP_KRLS)/length(RMSEP_KRLS)
