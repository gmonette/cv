
data("Auto", package="ISLR2")
m.auto <- lm(mpg ~ . - name - origin, data=Auto)
cv(selectStepAIC, Auto, seed=123, working.model=m.auto)

data("Prestige", package="carData")
m.pres <- lm(prestige ~ income + education + women,
             data=Prestige)
cvt <- cv(selectTrans, data=Prestige, working.model=m.pres, seed=123,
                predictors=c("income", "education", "women"),
                response="prestige", family="yjPower")
cvt
compareFolds(cvt)
coef(cvt, average=median, NAs=1) # NAs not really needed here
