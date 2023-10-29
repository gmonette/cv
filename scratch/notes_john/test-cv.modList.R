library(car)
library(cv)

m1 <- lm(prestige ~ income + education, data=Duncan)
m2 <- lm(prestige ~ income + education + type, data=Duncan)
m3 <- glm(prestige ~ income, data=Duncan)
m4 <- lm(prestige ~ income + education, data=Prestige)
m5 <- lm(prestige ~ income + education + type, data=Prestige)
m6 <- lm(women ~ income + education, data=Prestige)


lst1 <- models(m1=m1, m2=m2)
lst2 <- models(m1, m2)
lst3 <- models(m1, m3)
lst4 <- models(m1, m4)
lst5 <- models(m4, m5)
lst6 <- models(m4, m6)

cv(models(m1=m1, m2=m2), data=Duncan)

cv(models(m1=m1, m2=m2), data=Duncan, k="loo")


# ---- application to Auto data

data("Auto", package="ISLR2")
for (p in 1:10){
  assign(paste0("m.", p),
         lm(mpg ~ poly(horsepower, p), data=Auto))
}

cv.auto.10 <- cv(models(m.1, m.2, m.3, m.4, m.5,
                     m.6, m.7, m.8, m.9, m.10),
              data=Auto, seed=2120)
cv.auto.10[[2]]
cv.auto.loo <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 data=Auto, k="loo")
cv.auto.loo[[2]]
cv.mse.10 <- sapply(cv.auto.10, function(x) x[["adj CV crit"]])
cv.mse.loo <- sapply(cv.auto.loo, function(x) x[["CV crit"]])
plot(c(1, 10), range(cv.mse.10, cv.mse.loo), type="n",
     xlab="Degree of polynomial, p",
     ylab="Cross-Validated MSE")
lines(1:10, cv.mse.10, lwd=2, lty=1, col=2, pch=16, type="b")
lines(1:10, cv.mse.loo, lwd=2, lty=2, col=3, pch=17, type="b")
legend("topright", inset=0.02,
       legend=c("10-Fold CV", "LOO CV"),
       lwd=2, lty=2:1, col=3:2, pch=17:16)
