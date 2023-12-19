data(Auto, package="ISLR2")

for (p in 1:10){
  assign(paste0("m.", p),
         lm(mpg ~ poly(horsepower, p), data=Auto))
}

cv.auto.10.noci <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 data=Auto, seed=2120)

plot(cv.auto.10.noci)

cv.auto.10 <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 data=Auto, seed=2120, confint=TRUE,
                 level=.50)

plot(cv.auto.10)

cv.auto.10.someci <- cv.auto.10
cv.auto.10.someci[[2]] <- cv.auto.10.noci[[2]]
cv.auto.10.someci[[7]] <- cv.auto.10.noci[[7]]
plot(cv.auto.10.someci)

cv.auto.10.reps <- cv(models(m.1, m.2, m.3, m.4, m.5,
                             m.6, m.7, m.8, m.9, m.10),
                      data=Auto, seed=2120, reps=10,
                      confint=TRUE, level=0.50)

plot(cv.auto.10.reps)

data(Duncan, package="carData")
m1 <- lm(prestige ~ income + education, data=Duncan)
m2 <- lm(prestige ~ income + education + type, data=Duncan)
m3 <- lm(prestige ~ (income + education)*type, data=Duncan)
(cv.models <- cv(models(m1=m1, m2=m2, m3=m3),
                  data=Duncan, seed=7949, reps=5))
plot(cv.models)

(cv.models.ci <- cv(models(m1=m1, m2=m2, m3=m3),
                 data=Duncan, seed=5962, confint=TRUE,
                 level=0.50))
plot(cv.models.ci) # nb: n too small for accurate CIs


cv.auto.loo <- cv(models(m.1, m.2, m.3, m.4, m.5,
                         m.6, m.7, m.8, m.9, m.10),
                  data=Auto, k="loo")
cv.auto.loo
plot(cv.auto.loo)
