#№ PCA анализ для данных, полученных с помощью сенсоров и биоиндикации
library(readxl)
library(mdatools)

VibrioFisheri <- data.frame(read_excel("~/R/Chemometrics/Toxicity-ET.xlsx",sheet = "analysis_fisheri",
                                     range ='A1:Y27')) 
row.names(VibrioFisheri) <- c(VibrioFisheri[,1])
VibrioFisheri <-  VibrioFisheri[,-1]


m3 <- pca(VibrioFisheri, ncomp =8, scale = T)

par(mfrow = c(1,2)) 

#смотрю, какую часть дисперсии описывает каждая компонента
#выбор числа компонент
plotVariance(m3, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m3, type ='b', show.labels = TRUE)  # acum.var.

# строю грфик счетов и нагрузок
plotLoadings(m3, c(1,5), show.labels = T) # в скобках на графике процент описываемой вариации
plotScores(m3, c(1,5), show.labels = T) # раскрашены по цвету

# смотрю, какие переменные в большей степени повлияли на компоненты
m3$loadings

# change the number of component, без нормализации
# добавлю категории точек
plotResiduals(m3, show.labels = T,  norm = F, log = T, cgroup = "categories")

# !!почему-то не все категории указаны верно
c <- categorize(m3, m3$res$cal)
names(c) <- rownames(VibrioFisheri)




