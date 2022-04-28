library(readxl)
library(mdatools)
library(Metrics)



###############################################################################
# multisensory system
x <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                           range ='B2:L23'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                           range ='N2:Q23'))
y <- y[,-2]
rownames(y) <- seq(0,20, 1)
x <-  x[!(row.names(x) %in% c('0', '2', '16', '1')),]
y <-  y[!(row.names(y) %in% c('0', '2', '16', '1')),]

m_pb <- pls(x[,c(2,7,10,11)], y[,3], 10, cv = 1, scale = F)
plot(m_pb)
show(summary(m_pb, ncomp = 3))
plotYCumVariance(m_pb, type  = 'h', show.labels = TRUE)
plotXCumVariance(m_pb, type  = 'h', show.labels = TRUE)

# look for the outliers
plotRegcoeffs(m_pb, show.labels = T, ncomp = 3)
plotPredictions(m_pb, show.labels = T)
abline(a = 0, b = 1)


###############################################################################
# electrochemical sensor
# pls с чистыми репликами

ec <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "EC-Aver",
                            range ='A84:CY104')) 

row.names(ec) <- c(ec[,1])
ec <-  ec[,-c(1,2)] #remove names and Cb
ec <-  ec[!(row.names(ec) %in% c('18', '15')),]

#new model
m <- pls(ec[,2:101], ec[,1], 10, cv = 1, scale = F)
#plot(m, ncomp = 5)
plotRMSE(m)
plotYCumVariance(m)
show(summary(m))
#, ncomp = 7
#look for the outliers
plotPredictions(m, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m, show.labels = T, type = 'h')

################################################################################
###### pb data fusion

fused <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "pb data fusion",
                               range ='A3:DI23'))
row.names(fused) <- c(fused[,1])

fused <-  fused[,-c(1)] #remove names
fused <-  fused[!(row.names(fused) %in% c('2', '18')),]

#new model
m <- pls(fused[,2:112], fused[,1], 5, cv = 1, scale = T, center = T)
plotRMSE(m)
plotYCumVariance(m)
show(summary(m, ncomp =3))
#look for the outliers
plotPredictions(m, show.labels = T, lab.cex = 0.8, ncomp =3)
abline(a = 0, b = 1)
plotRegcoeffs(m, show.labels = T, type = 'h', lab.cex = 0.8, ncomp =3)


#vip = vipscores(m, ncomp = 3)
#plotVIPScores(m, type = "h", show.labels = TRUE)
#plotSelectivityRatio(m, ncomp = 3, type = "h", show.labels = TRUE)
#sr <- selratio(m, ncomp =3)

#m <- pls(fused[,2:112], fused[,1], 5, cv = 1, scale = T, center = T, exclcols = (sr < 0.02))



coef <- m[["coeffs"]][["values"]]
coef <- coef[,3,1]
#coef <- coef[,m[["ncomp.selected"]],1]

vec <- vector()

for (i in 1:length(coef)){
  if (abs(coef[i]) > 0.01) {
    vec <-  append(vec, i)
  }
}
rm(i)

#new model 
# vec +1 потому что в двтасете первым идет pb
new_m_n <- pls(fused[,vec+1], fused[,1], 10, cv = 1, scale = T)
plot(new_m_n)
show(summary(new_m_n))
#look for the outliers
plotPredictions(new_m_n, show.labels = T)
abline(a = 0, b = 1)

show(colnames(fused[,vec+1]))
plotRegcoeffs(new_m_n, type = 'h', show.labels = T)

rm(vec)

