# pls с чистыми репликами

library(readxl)
library(mdatools)
library(Metrics)

ec <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "EC-Aver",
                            range ='A84:CY104')) 

ec_aver <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "EC-Aver",
                                 range ='A110:CY130')) 


row.names(ec) <- c(ec[,1])
ec <-  ec[,-c(1,2)] #remove names and Cb
row.names(ec_aver) <- c(ec_aver[,1])
ec_aver <-  ec_aver[,-c(1,2)] #remove names and Cb
ec <-  ec[!(row.names(ec) %in% c('18', '15')),]
ec_aver <- ec_aver[!(row.names(ec_aver) %in% c('18', '15')),]

#new model
m <- pls(ec[,2:101], ec[,1], 10, cv = 1, scale = F)
#plot(m, ncomp = 5)
plotRMSE(m)
plotYCumVariance(m)
show(summary(m))
#, ncomp = 7
#look for the outliers
plotPredictions(m, show.labels = T, cex.lab = 3)
abline(a = 0, b = 1)
plotRegcoeffs(m, show.labels = T, type = 'h')


m_aver <- pls(ec_aver[,2:101],ec_aver[,1], 10, cv = 1, scale = F)
#plot(m_aver)
show(summary(m_aver))
plotRMSE(m_aver)
plotYCumVariance(m_aver)
#look for the outliers
plotPredictions(m_aver, show.labels = T, cex.lab = 3)
abline(a = 0, b = 1)
#plotRegcoeffs(m_aver, show.labels = T)