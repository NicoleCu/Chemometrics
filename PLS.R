# PLS

library(mdatools)
data('simdata')

xc <- simdata$spectra.c
yc <- simdata$conc.c[,2, drop = FALSE]
xt <- simdata$spectra.t
yt <- simdata$conc.t[,2, drop = FALSE]

attr(xc, "xaxis.values") <- simdata$wavelength
attr(xc, "xaxis.name") <- "Wavelength,nm"
boxplot(as.numeric(yc), as.numeric(yt))
mdaplot(xc, type='l', cgroup = as.vector(yc))
mdaplot(xt, type='l', cgroup = as.vector(yt))

#do a regression model, cv = 1 full cross-validation
# на самом деле, алгоритм сам находит оптимальное число компонент
# cv = 10 means  10 segments
# cv = list('rand', 5, 10) means random cross val, 5 segments, 10 repetitions
# cv = list("ven", 5)   venetian blinds (systematic split)

m <- pls(xc, yc, 5, cv = 1)
plot(m)
summary(m)

plotRMSE(m)
plotYVariance(m, type  = 'h', show.labels = TRUE)
plotPredictions(m)
plotRegcoeffs(m, 3, show.ci = T)   # blue lines -- доверительные интервалы
summary(m$coeffs)

#look for the outliers

plotXResiduals(m, 2)
plotXYResiduals(m,3)
plotXYScores(m, 2)


### prediction for the new dataset

# validation
r <- predict(m, xt, yt)
summary(m)
plot(r)
r$y.pred[, 3,1]     #predicted comcentrations, 3 component and 1 variable 

par(mfrow = c(2,1))
plotWeights(m, 2)
plotRegcoeffs(m,2)

