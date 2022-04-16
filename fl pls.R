# pls с чистыми репликами

library(readxl)
library(mdatools)
library(Metrics)

fl <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "FL-Aver",
                            range ='A109:Z130')) 

fl_aver <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "FL-Aver",
                            range ='A135:Z156')) 


row.names(fl) <- c(fl[,1])
fl <-  fl[,-c(1,3)] #remove names and Pb
row.names(fl_aver) <- c(fl_aver[,1])
fl_aver <-  fl_aver[,-c(1,3)] #remove names and Pb
fl <-  fl[!(row.names(fl) %in% c('11', '5', '15', '12', '0')),]
fl_aver <- fl_aver[!(row.names(fl_aver) %in% c('11', '5', '15', '12', '0')),]
# , '5','0', '15', '12'

#new model
m <- pls(fl[,2:24], fl[,1], 10, cv = 1, scale = F)
#plot(m, ncomp = 5)
plotRMSE(m)
plotYCumVariance(m)
show(summary(m, ncomp = 7))
#, ncomp = 7
#look for the outliers
plotPredictions(m, show.labels = T, cex.lab = 3, ncomp = 7)
abline(a = 0, b = 1)
plotRegcoeffs(m, show.labels = T, ncomp = 7)


m_aver <- pls(fl_aver[,2:24],fl_aver[,1], 10, cv = 1, scale = F)
#plot(m_aver)
show(summary(m_aver))
plotRMSE(m_aver)
plotYCumVariance(m_aver)
#look for the outliers
plotPredictions(m_aver, show.labels = T, cex.lab = 3)
abline(a = 0, b = 1)
plotRegcoeffs(m_aver, show.labels = T)
