## PLS на изначальных данных

library(readxl)
library(mdatools)
library(Metrics)

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='A1:R51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

# remove outliers 
row.names.remove <- c('4s','4b','8b','8s','28s','28b',"29s","29b",'30b','30s','24b')
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% row.names.remove),]


test <- c('7s','13s','11s','6b','15b','23s','19b','18b','20s')
col<-colnames(DafniaMagna)
DM_test <- DafniaMagna[test, col]
DM_train <- DafniaMagna[!(row.names(DafniaMagna) %in% test),]

#shapiro.test(as.integer(DM_test[,17]))


#new model
m <- pls(DM_train[,1:16], DM_train[,17], 10, cv = 1)
plot(m)
show(summary(m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)

r <- predict(m, DM_test[,1:16], DM_test[,17])
plot(r)
plotPredictions(r, show.label =T)
abline(0,1)
#summary(r)
RMSE <- rmse(r$y.pred['Comp 6'],r$y.ref)
show(rmse(r$y.pred[, 3,1], r$y.ref[,1] ))


##############################################################################
## PLS на данных с преобразованием


output <- as.data.frame(DafniaMagna[,17])
rownames(output) <- rownames(DafniaMagna)
predictors <- DafniaMagna[,1:16]
squared <- scale((predictors**2))
names <- paste(colnames(predictors), '**2')
colnames(squared) <- names

# поделим последний столбец на первый

multiplied <- matrix(1, nrow = nrow(predictors), ncol = 1)

for (i in 1:(ncol(predictors)-1)) {
    for (j in (i+1):ncol(predictors))  {
      n <- as.data.frame(predictors[,i]*predictors[,j])
      colnames(n) <- paste(as.character(i), as.character(j), sep = "*")
      multiplied <- cbind.data.frame(multiplied, n)
    }
}

rm(i, j)
multiplied <- multiplied[,-1]

divided <- matrix(1, nrow = nrow(predictors), ncol = 1)

for (i in 1:(ncol(predictors)-1)) {
  for (j in (i+1):ncol(predictors)) {
    n <- as.data.frame(predictors[,i]/predictors[,j])
    colnames(n) = paste(as.character(i), as.character(j), sep = "/")
    divided <- cbind.data.frame(divided, n)
  }
}

rm(i, j)
divided <- divided[,-1]

predictors <- cbind(scale(predictors), squared, scale(multiplied), scale(divided))

x_test <- predictors[test,]
y_test <- as.data.frame(output[test,])
rownames(y_test) <- test
colnames(y_test) <-  'Dafnia'
xy_test <- cbind(x_test, y_test)
x_train <- predictors[!(row.names(DafniaMagna) %in% test),]
y_train <- as.data.frame(output[!(row.names(DafniaMagna) %in% test),])
rownames(y_train) <- rownames(x_train)
colnames(y_train) <-  'Dafnia'
xy_train <- cbind(x_train, y_train)

#new model
new_m <- pls(x_train, y_train, 10, cv = 1, scale = F)
plot(new_m)
show(summary(new_m))
#look for the outliers
plotPredictions(m, show.labels = T)
abline(a = 0, b = 1)

r_new <- predict(new_m, x_test, y_test)
plot(r_new)
plotPredictions(r, show.label =T)
abline(0,1)
summary(r_new)
RMSE <- rmse(r_new$y.pred['Comp 3'],r_new$y.ref)
show(rmse(r_new$y.pred[, 3,1], r_new$y.ref[,1] ))



################################################################################
## RF для преобразованных данных

library(e1071)


OptModelsvm <- tune(svm, Dafnia ~., data=xy_train, kernel = 'radial',
                    ranges=list(cost = c(1:10,100),  epsilon = seq(0.01, 0.15, 0.01),gamma = c(0.0001, 0.0005, 0.0001)), 
                    tunecontrol = tune.control(cross = 30), scale = F)  
print(OptModelsvm)
#Plot the performance of SVM Regression model
#plot(OptModelsvm)
#The best model is the one with lowest MSE. 

## Select the best model out of trained models and compute RMSE
#Find out the best model
BstModel <- OptModelsvm$best.model
show(BstModel)

#Predict Y using best model
pred <- predict(BstModel, xy_test)
st <- lm(pred ~ xy_test[,17])
plot(xy_test[,17], pred, 
     main = 'SVM',
     ylim = c(0, 100), xlim = c(0, 100), 
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(xy_test[,17], pred,
     labels = (row.names(xy_test)),
     cex = 0.6, pos = 4, col = "gray")
grid()

rmsep <- rmse(xy_test[,17], pred)
bias <- bias(xy_test[,17], pred)
slope <- st$coefficients[2]

#R2 <- 1 - (sum((DM_test[,17] - pred)^2))/(sum((DM_test[,17] - mean(DM_test[,17]))^2))
show(rbind(rmsep, bias, slope))





