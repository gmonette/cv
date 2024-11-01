library(cv)

phi.x <- 0.5
sigma.eps.x <- 1
phi.1 <- 0.5
phi.2 <- 0.3
sigma <- 1
sigma.eps <- 1
b1 <- 1
b2 <- 1
a <- 1
b.x <- 0.01
n <- 1e3

set.seed(94962)
x <- numeric(n)
x[1] <- rnorm(1, 0, sigma.eps.x)
for (i in 2:n){
  x[i] <- phi.x*x[i - 1] + rnorm(1, 0, sigma.eps.x)
}
acf(x)
pacf(x)
x <- x + b.x*1:n
plot(x, type="l")

eps <- numeric(n)
eps[1] <- rnorm(1, 0, sigma)
eps[2] <- phi.1*eps[1] + rnorm(1, 0, sigma)
for (i in 3:n){
  eps[i] <- phi.1*eps[i - 1] + phi.2*eps[i - 2] + rnorm(1, 0, sigma)
}
acf(eps)
pacf(eps)

y <- a + b1*x + b2*x^2 + eps
plot(y, type="l")

D <- data.frame(x, y)

m.ls <- lm(y ~ poly(x, 2, raw=TRUE), data=D)
summary(m.ls)
acf(residuals(m.ls))
pacf(residuals(m.ls))

summary(Arima(~ x, order=c(1, 0, 0), data=D))

m.1 <- Arima(y ~ x, order=c(2, 0, 0), data=D)
summary(m.1)
plot(m.1)
testArima(m.1)
plot(effects::Effect("x", m.1, residuals=TRUE),
     partial.residuals=list(span=0.125))

m.2 <- Arima(y ~ poly(x, 2, raw=TRUE),
             order=c(2, 0, 0), data=D)
summary(m.2)
plot(m.2)
testArima(m.2)
plot(effects::Effect("x", m.2, residuals=TRUE),
     partial.residuals=list(span=0.125))

m.3 <- Arima(y ~ poly(x, 3, raw=TRUE),
             order=c(2, 0, 0), data=D)
summary(m.3)
plot(m.3)
testArima(m.3)
plot(effects::Effect("x", m.3, residuals=TRUE),
     partial.residuals=list(span=0.125))


system.time(cv.m <- cv(models(linear=m.1, quadratic=m.2, cubic=m.3),
           lead=1:5, data=D))
summary(cv.m)
plot(cv.m, legend=list(x=3, y=3.5))

system.time(cv.m.p <- cv(models(linear=m.1, quadratic=m.2, cubic=m.3),
                       lead=1:5, data=D, ncores=2))
all.equal(cv.m, cv.m.p)

system.time(cv.m.c <- cv(models(linear=m.1, quadratic=m.2, cubic=m.3),
           lead=1:5, data=D, fold.type="cumulative"))
summary(cv.m.c)
plot(cv.m.c, legend=list(x=3, y=40))

system.time(cv.m.c.100 <- cv(models(linear=m.1, quadratic=m.2, cubic=m.3),
                             lead=1:5, data=D, fold.type="cumulative", k=100))
summary(cv.m.c.100)
plot(cv.m.c.100, legend=list(x=3, y=40))

system.time(cv.m.pr <- cv(models(linear=m.1, quadratic=m.2, cubic=m.3),
                         lead=1:5, data=D, fold.type="preceding", k=10))
summary(cv.m.pr)
plot(cv.m.pr, legend=list(x="topright"))
