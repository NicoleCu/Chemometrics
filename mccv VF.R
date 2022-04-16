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
vec_PLS <- vector()
vec_RF <- vector()
vec_KRLS <- vector()
vec_SVM <- vector()


for (i in 1:100) {
  s <- sample(1:nrow(VibrioFisheri), prop)
  m <- pls(VibrioFisheri[s,1:23], VibrioFisheri[s,'EC50'], 8, cv = 1)
  r <- predict(m, VibrioFisheri[-s,1:23], VibrioFisheri[-s,'EC50'])
  RMSEP_PLS <- rmse(r$y.pred[, r[["ncomp.selected"]],1], r$y.ref[,1])
  vec_PLS[i] <- RMSEP_PLS
  
  OptModelsvm <- tune(svm,  EC50 ~ ., data = VibrioFisheri[s,], kernel = 'radial',
                      ranges=list(cost = c(1,10,100),  
                                  epsilon = seq(0.01, 0.15, 0.01),
                                  gamma = c(0.0001, 0.0005, 0.0001)), 
                      tunecontrol = tune.control(cross = prop), scale = F)
  BstModel <- OptModelsvm$best.model
  pred <- predict(BstModel, VibrioFisheri[-s,])
  RMSEP_SVM <- rmse(VibrioFisheri[-s,'EC50'], pred)
  vec_SVM[i] <- RMSEP_SVM
  
  RF <- randomForest(EC50 ~., data = VibrioFisheri[s,], xtest = VibrioFisheri[-s,1:23], ytest = VibrioFisheri[-s,'EC50'],  
                     mtry = 5, ntree = 500, nodesize = 5, importance = T)
  RMSEP_RF <- rmse(VibrioFisheri[-s,'EC50'], RF[["test"]][["predicted"]])
  vec_RF[i] <-  RMSEP_RF
  
  KRLS <- krls(VibrioFisheri[s,1:23], as.vector(VibrioFisheri[s,'EC50']), vcov = T, derivative = T)
  KRLS_pred <- predict(KRLS,newdata = VibrioFisheri[-s,1:23], se.fit=TRUE)
  RMSEP_KRLS <- rmse(VibrioFisheri[-s,'EC50'], KRLS_pred$fit)
  vec_KRLS[i] <- RMSEP_KRLS
}


av_PLS_error <- sum(vec_PLS)/length(vec_PLS)
av_SVM_error <- sum(vec_SVM)/length(vec_SVM)
av_RF_error <- sum(vec_RF)/length(vec_RF)
av_KRLS_error <- sum(vec_KRLS)/length(vec_KRLS)
