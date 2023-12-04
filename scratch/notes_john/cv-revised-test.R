# commented output is from previous version of cv
# should match for criterion=mse

data("Auto", package="ISLR2")
lm.fit <- lm(mpg ~ horsepower, data=Auto)
glm.fit <- glm(mpg ~ horsepower,data=Auto)

cv(lm.fit, k="loo", method="naive")
# n-Fold Cross Validation
# method: naive
# criterion: mse
# cross-validation criterion = 24.23151
# bias-adjusted cross-validation criterion = 24.23114
# full-sample criterion = 23.94366
cv(lm.fit, k="loo", method="naive")

cv(lm.fit, k=10, seed=123)
# R RNG seed set to 123
# 10-Fold Cross Validation
# method: Woodbury
# criterion: mse
# cross-validation criterion = 24.44621
# bias-adjusted cross-validation criterion = 24.41958
# full-sample criterion = 23.94366
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
# n-Fold Cross Validation based on 18 {Subject} clusters
# cross-validation criterion = 2460.604
# bias-adjusted cross-validation criterion = 2454.627
# full-sample criterion = 2251.398
cv(fm1, clusterVariables="Subject", ncores=2)

cv(fm1, clusterVariables="Subject", k=5, seed=123)
# R RNG seed set to 123
# 5-Fold Cross Validation based on 18 {Subject} clusters
# cross-validation criterion = 2397.44
# bias-adjusted cross-validation criterion = 2379.979
# full-sample criterion = 2251.398
cv(fm1, clusterVariables="Subject", k=5, seed=123, ncores=2)


library(MASS)
data(Auto, package="ISLR2")
m <- lm(mpg ~ . - name - origin, data=Auto)

cvSelect(selectStepAIC, Auto, k=5, seed=123,
         model=m, method="naive")
# R RNG seed set to 123
# 5-Fold Cross Validation
# cross-validation criterion = 12.06025
# bias-adjusted cross-validation criterion = 12.0452
# full-sample criterion = 11.65549
cvSelect(selectStepAIC, Auto, k=5, seed=123,
         model=m, method="Woodbury")
res <- cvSelect(selectStepAIC, Auto, k=5, seed=123,
         model=m, method="Woodbury", ncores=2)
res
compareFolds(res)

# need to debug() to observe:
cvSelect(selectStepAIC, Auto, k=5, seed=123,
         model=m, method="naive",
         response.expr=expression(log(mpg)))

data(Prestige, package="carData")
m <- lm(prestige ~ income + education, data=Prestige)
(res <- cvSelect(selectTrans, data=Prestige, model=m, seed=123,
                predictors=c("income", "education"),
                response="prestige", family="bcPower"))
# R RNG seed set to 123
# 10-Fold Cross Validation
# cross-validation criterion = 55.02292
# bias-adjusted cross-validation criterion = 54.88635
# full-sample criterion = 50.75109
compareFolds(res)
