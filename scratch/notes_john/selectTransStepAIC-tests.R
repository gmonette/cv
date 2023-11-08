library(cv)
library(car)
library(effects)
source("./scratch/notes_john/selectTransStepAIC.R")

data("Auto", package="ISLR2")

Auto$Year <- as.factor(Auto$year)
Auto$origin <- factor(Auto$origin,
                      labels=c("America", "Europe", "Japan"))
rownames(Auto) <- make.names(Auto$name, unique=TRUE)
Auto$name <- NULL

m.auto <- lm(mpg ~ . - year, data=Auto)
cvs <- cvSelect(selectTransStepAIC, data=Auto, seed=76692, model=m.auto,
                predictors=c("cylinders", "displacement", "horsepower",
                             "weight", "acceleration"),
                response="mpg", AIC=FALSE)
cvs
compareFolds(cvs)

predictors <- c("cylinders", "displacement", "horsepower",
                "weight", "acceleration")
tr.x <- powerTransform(Auto[, predictors])
summary(tr.x)
A <- Auto
powers <- tr.x$roundlam
for (pred in predictors){
  A[, pred] <- bcPower(A[, pred], lambda=powers[pred])
}
head(A)

m <- update(m.auto, data=A)
summary(powerTransform(m))
m <- update(m, log(mpg) ~ .)

m.step <- MASS::stepAIC(m, k=log(nrow(A)), trace=FALSE)
summary(m.step)
Anova(m.step)
plot(predictorEffects(m.step, residuals=TRUE))

m.cv <- lm(log(mpg) ~ log(horsepower) + log(weight) + Year,
           data=Auto)
summary(m.cv)
Anova(m.cv)
plot(predictorEffects(m.cv, residuals=TRUE))

anova(m.cv, m.step)

mse(getResponse(m.cv), fitted(m.cv))

