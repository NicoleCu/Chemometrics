library(mdatools)
data(people)

Height <- people[,'Height']
Weight <-  people[,'Weight']

#ищу ковариацию. По свойствам оба выражения равнозначны
cov(Height, Weight)
cov(Weight, Height)

#дисперсию можно найти таким путем
cov(Height, Height)

#covariance matrix можно найти таким путем
C <- cov(people[,1:3])
C

#или таким: crossprod ~ X'X (первая транспонирована)
#tcrossprod ~ XX' (вторая транспонирована)
C2 <- crossprod(Xm)/(nrow(Xm)-1) 
C2

#mean centering
Xm <- scale(people[,1:3], center = TRUE, scale = FALSE)
head(Xm)

# проверка средних, 1 - по столбцам, 2 - по строкам
apply(people[,1:3], 2, mean)

# строю график центрированных данных
plot(Xm[,1], Xm[,2],
     main = "Center data",
     xlab = 'Height',
     ylab = 'Weight')
#провдить линию на графике, первое и второе значение intercept and slope соотвественно
# h =... линия параллельно оси абсцисс, v =... параллельно ординат
abline(h=0, col = 'Gray', lty = 2)
abline(v=0, col = 'Gray', lty = 2)
#рисую линию тренда (линейной регрессии), в lm содержится slope and intercept
abline(lm(Xm[,2]~Xm[,1]))


