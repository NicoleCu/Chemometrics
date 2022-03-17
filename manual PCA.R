## ѕостроение PC1 дл€ двумерных данных, функци€ написана вручную

library(mdatools)
data(people)

#создаю функцию, первое данные, второе кол-во компонент дл€ вычислени€, 
#третье нужно ли стандартизировать данные
nipals <- function(X, A=1, scale = FALSE){
  
  objnames <- rownames(X)  #сохран€ю название строк(объектов)
  varnames <- colnames(X)  # сохран€ю название столбцов 
  
  #mean centering, на графике сместила 0 в центр
  X <- scale(X, center = TRUE, scale = scale)
 
  #total sum of squared distances (данные стандартизованы) 
  SStot <- sum(X^2)
  #подготавливаю матрицы дл€ scores and loadings
  scores <- matrix(0, nrow(X), A)
  loadings <- matrix(0, ncol(X), A)
  
  #compute individual explained variance
  expvar <- rep(0, A)
  
  for (a in 1:A) {
    # initialisation of scores дл€ каждой из компонент
    # нахожу, в каком столбце максимальное стандратное отклонение и создаю вектор t
    col.index <- which.max(apply(X, 2, sd))
    t <- X[,col.index, drop = FALSE]
  
    for (i in 1:10) {
    
    #compute loadings
    p <- cov(X,t)
    p <- p/sqrt(sum(p^2))    # normalize
    #compute scores через матричное умножение
    t <- X%*%p
    }
    
    scores[,a] <- t
    loadings[,a] <- p
    
    #TP'
    tpt <- tcrossprod(t,p)
    X <- X - tpt  # remove the variation that are already explained
    # это действие оставл€ет остаточную вариацию, которую еще надо исследовать
    
    #explained var. %
    expvar[a] <- sum(tpt^2)/SStot
  }
  rownames(scores) <- objnames
  rownames(loadings) <- varnames
  names(expvar) <- colnames(scores) <- colnames(loadings) <- paste0("PC", 1:A)
  return(list(scores = scores, loadings = loadings, expvar = expvar))
}

X <- people[, c('Height', 'Weight', 'Region')]
m <- nipals(X, A=2, scale =TRUE)

residual <- c(1-m$expvar)
#accumulated variance
cumsum(m$expvar)

m$loadings   # можно посмотреть, кака€ компонента описывает какую переменную
par(mfrow = c(1,3)) # соединю 3 нижних графика
barplot(m$expvar)
plot(m$loadings, col = 'Blue', xlim = c(-1,1), ylim = c(-1,1))
abline(h=0, col = 'Gray', lty = 1)
abline(v=0, col = 'Gray', lty = 1)
grid()
text(m$loadings[,1], m$loadings[,2], rownames(m$loadings), cex = 0.6, pos = 2)

plot(m$scores, col = 'Blue')
abline(h=0, col = 'Gray', lty = 1)
abline(v=0, col = 'Gray', lty = 1)
grid()
text(m$scores[,1], m$scores[,2], rownames(m$scores), cex = 0.55, pos = 1)


#как нарисовать единичный вектор и PC1 в реальных координатах
#arrows(0,0, p[1], p[2]) #можно умножить на 100 просто, чтобы его было видно
#abline(a=0, b = p[2]/p[1], col ='red')

#compute coordinates of projection on PC1
#tpt <- tcrossprod(t,p)
#points(tpt, pch = 3)