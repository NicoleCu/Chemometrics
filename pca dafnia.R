#№ PCA анализ для данных, полученных с помощью сенсоров и биоиндикации
library(readxl)
library(mdatools)

DafniaMagna <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='A1:R51')) 
row.names(DafniaMagna) <- c(DafniaMagna[,1])
DafniaMagna <-  DafniaMagna[,-1]

m <- pca(DafniaMagna, ncomp =6, scale = T)

par(mfrow = c(1,2)) 

#смотрю, какую часть дисперсии описывает каждая компонента
#выбор числа компонент
plotVariance(m, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m, type ='b', show.labels = TRUE)  # acum.var.

# строю грфик счетов и нагрузок
plotLoadings(m, c(1,2), show.labels = T) # в скобках на графике процент описываемой вариации
plotScores(m, c(1,2), show.labels = T) # раскрашены по цвету

# смотрю, какие переменные в большей степени повлияли на компоненты
m$loadings

# change the number of component, без нормализации
# добавлю категории точек
plotResiduals(m, show.labels = T,  norm = F, log = T, cgroup = "categories")

# !!почему-то не все категории указаны верно
c <- categorize(m, m$res$cal)
names(c) <- rownames(DafniaMagna)




## PCA анализ для данных, полученных с помощью биоиндикации отдельно
DafniaMagna_bioind <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_dafnia",
                                     range ='R1:R51')) 
row.names(DafniaMagna_bioind) <- row.names(DafniaMagna)

m2 <- pca(DafniaMagna, ncomp = 6, center =T, scale = T)
par(mfrow = c(1,2)) 

#смотрю, какую часть дисперсии описывает каждая компонента
#выбор числа компонент
plotVariance(m2, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m2, type ='b', show.labels = TRUE)  # acum.var.

# строю грфик счетов и нагрузок
plotLoadings(m2, c(1,6), show.labels = T) # в скобках на графике процент описываемой вариации
plotScores(m2, c(1,6), show.labels = T) # раскрашены по цвету

# смотрю, какие переменные в большей степени повлияли на компоненты
m$loadings

# change the number of component, без нормализации
# добавлю категории точек
plotResiduals(m2, show.labels = T, norm = F, log = T, cgroup = "categories")
