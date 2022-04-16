library(readxl)
library(mdatools)
library(Metrics)

# data import 

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/processed data DafniaMagna.xlsx",
                                     sheet = "Sheet1", range ='A1:JN51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]
colnames(DafniaMagna) <- read_excel("~/R/Chemometrics/processed data DafniaMagna.xlsx",
                                    sheet = "Sheet1", range ='B1:JN1', col_names = F)

# remove outliers 
DafniaMagna <-  DafniaMagna[!(row.names(DafniaMagna) %in% 
                                c('4s','4b','8b','8s','28s','28b',"29s","29b",'30b','30s','24b')),]


# new dataset of predictors
#predictors <- scale(DafniaMagna[,1:272])
#colnames(DM) <- colnames(DafniaMagna)

DM_x <- DafniaMagna[,1:272]
DM_y <- DafniaMagna[,'Dafnia']



##############################################################################
## PLS на данных с преобразованием

#new model
new_m <- pls(DM_x, DM_y, 10, cv = 1, scale = T, center =T)
plot(new_m)
show(summary(new_m))    
#look for the outliers
plotPredictions(new_m, show.labels = T)
abline(a = 0, b = 1)

#plotYVariance(new_m)	

plotRegcoeffs(new_m, type = 'h', labels = T)

coef <- new_m[["coeffs"]][["values"]]
coef <- coef[,new_m[["ncomp.selected"]],1]
vec <- vector()

for (i in 1:length(coef)){
  if (abs(coef[i]) > 0.015) {
    vec <-  append(vec, i)
  }
}

#new model
new_m_n <- pls(DM_x[,vec], DM_y, 10, cv = 1, scale = T)
plot(new_m_n)
show(summary(new_m_n))
#look for the outliers
plotPredictions(new_m_n, show.labels = T)
abline(a = 0, b = 1)

show(colnames(DafniaMagna[,vec]))
plotRegcoeffs(new_m_n, type = 'h')
