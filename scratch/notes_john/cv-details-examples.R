data("Auto", package="ISLR2")
m.auto <- lm(mpg ~ poly(horsepower, 2), data=Auto)
cv.lm.loo <- cv(m.auto,  k="loo")
cv.lm.loo$details
cv.lm.wood <- cv(m.auto, seed=1234)
cv.lm.wood$details
cv.lm.naive <- cv(m.auto, seed=1234, method="naive")
cv.lm.naive$details
all.equal(cv.lm.wood$details, cv.lm.naive$details)
max(mapply(cv.lm.wood$details$coefficients,
       cv.lm.naive$details$coefficients,
       FUN=function(x, y) max(abs(x - y)/x)))
cbind(unlist(cv.lm.wood$details$coefficients),
      unlist(cv.lm.naive$details$coefficients))

data("Mroz", package="carData")
m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
cv.glm.naive <- cv(m.mroz, criterion=BayesRule, seed=123)
cv.glm.naive$details
cv.glm.wood <- cv(m.mroz, criterion=BayesRule, seed=123,
                  method="Woodbury")
cv.glm.wood$details
all.equal(cv.glm.wood$details, cv.glm.naive$details)
all.equal(lapply(cv.glm.wood$details$coefficients, round),
          lapply(cv.glm.naive$details$coefficients, round))
max(abs(mapply("-", cv.glm.wood$details$coefficients,
       cv.glm.naive$details$coefficients)))

library("lme4")
# from ?lmer:
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
cvloo <- cv(fm1, clusterVariables="Subject") # LOO CV of clusters
cvcase <- cv(fm1, seed=447) # 10-fold CV of cases
cvclus <- cv(fm1, clusterVariables="Subject", k=5,
   seed=834) # 5-fold CV of clusters
cvcase$details
cvclus$details
lapply(cvcase$details$coefficients,
       function(x) rownames(x$Subject))


