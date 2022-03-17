library(mdatools)
data(people)

X <- people

m <- pca(X, ncomp = 4, scale = TRUE)

m$loadings

par(mfrow = c(1,2)) 
plotLoadings(m, c(1,4), show.labels = T) # в скобках на графике процент описываемой вариации
# 
plotScores(m, c(1,4), show.labels = T) # раскрашены по цвету

# barploе, можно нарисовать еще много вариаций
plotLoadings(m, show.labels = T, type = 'h')

summary(m)
plotVariance(m, type ='h', show.labels = TRUE) #expl.var
plotCumVariance(m, type ='b', show.labels = TRUE)  # acum.var.

m$res
cumresidual <- 100 - m$res$cal$cumexpvar
cumresidual
barplot(as.numeric(cumresidual))

## рассто€ни€
# distance plot, автоматически использует стоько компонент, сколько задано изначально
plotResiduals(m)
# change the number of component, без нормализации
plotResiduals(m, c(1,2), show.labels =T, norm = F, log = T)
# добавлю категории точек
plotResiduals(m, c(1,2), show.labels = T, cgroup = "categories")

# compute the categories
c <- categorize(m, m$res$cal)
names(c) <- rownames(people)

# добавила 2 новые точки
P1 <- c(159, 50, +1, 36, 20, 15000, 10, 10, +1, 95, +1, 110)
P2 <- c(159, 50, +1, 36, 20, 30000, 10, 10, +1, 95, +1, 120 )
new <- rbind(P1, P2)

r <- predict(m, new) # use ready model for a new data
plotScores(m, res=list("cal" = m$res$cal, "new" = r), show.labels = TRUE)

plotResiduals(m, res=list("cal" = m$res$cal, "new" = r), show.labels = TRUE)

# смотрим, что там с residuals, потому новые данные попали в выбросы

E <- r$residuals
colnames(E) <- colnames(people)
mdaplotg(E, type ='h', show.labels = T)

# смотрим на классы 2ух новых объектов
cnew <- categorize(m, r)
cnew

# выставл€ю лимиты дл€ модели
m <- pca(people, ncomp = 6, scale = TRUE, lim.type = 'ddmoments', alpha = 0.05, gamma = 0.01)
plotResiduals(m,4)
