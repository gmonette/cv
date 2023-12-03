data("Auto", package="ISLR2")
lm.fit <- lm(mpg ~ horsepower, data=Auto)
glm.fit <- glm(mpg ~ horsepower,data=Auto)

cv(lm.fit, k="loo", method="naive")
cv(lm.fit, k="loo", method="naive")

cv(lm.fit, k=10, seed=123)
cv(lm.fit, k=10, seed=123, ncores=2)

cv(lm.fit, k=10, seed=123, criterion=rmse)
sapply(cv(lm.fit, k=10, seed=123, criterion=rmse)[1:3],
       function(x) x^2)

cv(lm.fit, k=10, seed=123, method="naive")
cv(lm.fit, k=10, seed=123, ncores=2, method="naive")

cv(glm.fit, k=10, seed=123)
cv(glm.fit, k=10, seed=123, ncores=2)
cv(glm.fit, k=10, seed=123, method="Woodbury")
cv(glm.fit, k=10, seed=123, method="Woodbury", ncores=2)

library(lme4)
data("sleepstudy", package="lme4")
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

cv(fm1, clusterVariables="Subject")
cv(fm1, clusterVariables="Subject", ncores=2)

cv(fm1, clusterVariables="Subject", k=5, seed=123)
cv(fm1, clusterVariables="Subject", k=5, seed=123, ncores=2)


library(MASS)
data(Auto, package="ISLR2")
m <- lm(mpg ~ . - name - origin, data=Auto)

cvSelect(selectStepAIC, Auto, k=5, seed=123,
         model=m, method="naive")
cvSelect(selectStepAIC, Auto, k=5, seed=123,
         model=m, method="Woodbury")
res <- cvSelect(selectStepAIC, Auto, k=5, seed=123,
         model=m, method="Woodbury", ncores=2)
res
coef(res)

cvSelect(selectStepAIC, Auto, k=5, seed=123,
         model=m, method="naive",
         response.expr=expression(log(mpg)))
