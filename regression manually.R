data("swiss")

swiss <-  as.data.frame(swiss)
cor(swiss)

### function for computing statistics
getStat <- function(y, yp) {
  e <- y - yp
  bias <- mean(e)     #literally 0    
  slope <- cov(y, yp)/var(yp)
  r2 <- 1 - sum(e^2)/sum((yp-mean(yp))^2)
  rmse <- sqrt(sum(e^2))/length(yp)
  stat <- c(bias = bias, slope = slope, r2 = r2, rmse = rmse)
  
  return(round(stat, 4))
}


###calibration and test sets
xc <- swiss[1:37, 4]
yc <-  swiss[1:37, 1]
xt <- swiss[38:47, 4]
yt <-  swiss[38:47, 1]

###find complexity
rmsec <- rep(0,8)
rmsep <- rep(0,8)
for (n in 1:8) {
  m <- lm(yc ~ poly(xc, n))
  ycp <- predict(m)
  ytp <- predict(m, data.frame(xc = xt))
  rmsec[n] <- getStat(yc, ycp)[4] 
  rmsep[n] <- getStat(yt, ytp)[4] 
  
}

plot(rmsec, type ='b', col ='blue', ylim =c(2,100))
lines(rmsep, type ='b', col ='red')

### модель полиномиальной регресии, n степень, если n = 1, то это линейна€
n <- 1
m <- lm(yc ~ poly(xc, n))
# get predictions for calibration set
ycp <- predict(m)
#predictions for the test set
ytp <- predict(m, data.frame(xc = xt))

# calculate statistics
# for cal. set and for the test set
stat <- rbind(getStat(yc, ycp), getStat(yt, ytp))
rownames (stat) <- c('Calibration', "Test")
show(stat)


### Plot part
par(mfrow = c(1,3))
# как выгл€дит распределение точек
plot(xc, yc, col = 'blue')
points(xt, yt, col ='red')

#график модели полиномиальной регресии, n степень, если n = 1, то это линейна€
plot(xc, yc, col = 'blue')
points(xt, yt, col ='red')
#чертим линию 
xm <- seq(min(xc), max(xc), length.out = 100)
ym <- predict(m, data.frame(xc = xm))
lines(xm, ym, col ='orange')

#predicted vs measured plot for calibration set 
plot(yc, ycp, main = "Predictions")
abline(a = 0, b = 1, col ='black', lty = 2) #идеальна€ лини€
abline(lm(ycp ~ yc), col = 'blue') # реальна€ лини€

#add predicted vs measured plot for the test set
points(yt, ytp, col = 'red')
abline(lm(ytp ~ yt), col = 'red') # реальна€ лини€

# cross-validation
n <- 1
ycv <- rep(0, length(yc))
for (i in 1:length(xc)) {
  xlc <- xc[-i]
  ylc <- yc[-i]
  xlt <- xc[i]

  mcv <- lm(ylc ~ poly(xlc, n))
  ycv[i] <- predict(mcv, data.frame(xlc = xlt))
}
plot(yc, ycv, col = "green")
abline(lm(ycv ~ yc), col = 'green')
abline(a=0, b=1, lty =2)

stat <- getStat(yc, ycv)
show(round(stat,4))
