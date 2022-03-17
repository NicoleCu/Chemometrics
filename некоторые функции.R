d <- pnorm(q = 70, mean = 100, sd = 15, lower.tail = FALSE ) - pnorm(q = 112, mean = 100, sd = 15, lower.tail = FALSE) 
d

# стандартная ошибка
sd <- 5
n <- 100
SE <- sd/sqrt(n)
mean <- 10
#вероятность 99%
q <- qnorm(0.99508)
left <- mean - q*SE
right <- mean + q*SE
left
right
