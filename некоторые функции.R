d <- pnorm(q = 70, mean = 100, sd = 15, lower.tail = FALSE ) - pnorm(q = 112, mean = 100, sd = 15, lower.tail = FALSE) 
d

# ����������� ������
sd <- 5
n <- 100
SE <- sd/sqrt(n)
mean <- 10
#����������� 99%
q <- qnorm(0.99508)
left <- mean - q*SE
right <- mean + q*SE
left
right
