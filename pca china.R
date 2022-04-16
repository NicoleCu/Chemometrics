# PCA for China data

library(readxl)
library(mdatools)

############################ FL sensor

# data import
fl <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "data of FL",
                                     range ='A1:Z85')) 
row.names(fl) <- c(fl[,1])
fl <-  fl[,-1]

m <- pca(fl[,3:25], ncomp = 3, scale = F, center = T)

#смотрю, какую часть дисперсии описывает каждая компонента
#выбор числа компонент
#plotVariance(m, type ='h', show.labels = TRUE) #expl.var
#plotCumVariance(m, type ='b', show.labels = TRUE)  # acum.var.

# строю грфик счетов и нагрузок
plotLoadings(m, c(1,2), show.labels = T) # в скобках на графике процент описываемой вариации
plotScores(m, c(1,2), show.labels = T)

# change the number of component, без нормализации
# добавлю категории точек
plotResiduals(m, show.labels = T, log = F, norm = T, cgroup = "categories")



##################### EC sensor


# data import
ec <- data.frame(read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "data of EC",
                            range ='A1:CY61')) 
row.names(ec) <- c(ec[,1])
ec <-  ec[,-1]

m_ec <- pca(ec[,3:102], ncomp = 3, scale = F, center = T)
#выбор числа компонент
#plotVariance(m, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m, type ='b', show.labels = TRUE)  # acum.var.

# строю грфик счетов и нагрузок
plotLoadings(m_ec, c(1,2), show.labels = T) # в скобках на графике процент описываемой вариации
plotScores(m_ec, c(1,2), show.labels = T)

# change the number of component, без нормализации
# добавлю категории точек
plotResiduals(m_ec, show.labels = T, log = T, norm = F, cgroup = "categories")
