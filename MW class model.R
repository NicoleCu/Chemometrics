# model data, building a classification model.

library(readxl)
library(mdatools)
library(Metrics)
library(e1071)
library(randomForest)

# importing and preparing data

MW <- data.frame(read_excel("~/R/Chemometrics/ModelWaters 2.xlsx",
                           range ='B2:M23'))
MW[,12] <- factor(MW[,12])
rownames(MW) <- seq(0,20, 1)
colnames(MW) <- c(seq(1,11, 1), 'Class')

################################################################################
# PLS-DA with outliers. 

m.all <-  plsda(MW[,1:11], MW[,12], 3, cv = 1, scale = F)
summary(m.all)
getConfusionMatrix(m.all$cvres)
plotRegcoeffs(m.all, show.labels = T)
plotPredictions(m.all, show.labels = T)

m.all_n <-  plsda(MW[,c(2,5,7,9,10,11)], MW[,12], 3, cv = 1, scale = F)
summary(m.all_n)
getConfusionMatrix(m.all_n$cvres)
plotRegcoeffs(m.all_n, show.labels= T)
plotPredictions(m.all_n, show.labels = T)


################################################################################
# PLS-DA without oitliers

MW <-  MW[!(row.names(MW) %in% c('19','18', '9','11')),]


m <-  plsda(MW[,1:11], MW[,12], 3, cv = 1)
summary(m)
getConfusionMatrix(m$cvres)
plotRegcoeffs(m, show.labels = T)
plotPredictions(m, show.labels = T)



