library(cv)
library(car)
library(effects)
source("./scratch/notes_john/selectTransStepAIC.R")

data("Auto", package="ISLR2")
xtabs(~ year, data=Auto)
xtabs(~ cylinders, data="Auto")

Auto$year <- as.factor(Auto$year)
Auto$origin <- factor(Auto$origin,
                      labels=c("America", "Europe", "Japan"))
Auto$cylinders <- factor(Auto$cylinders,
                         labels=c("3.4", "3.4", "5.6", "5.6", "8"))
rownames(Auto) <- make.names(Auto$name, unique=TRUE)
Auto$name <- NULL

scatterplotMatrix(~ mpg  + displacement + horsepower
                  + weight + acceleration, data=Auto)

m.auto <- lm(mpg ~ ., data = Auto)
summary(m.auto)
Anova(m.auto)
crPlots(m.auto)

predictors <- c("displacement", "horsepower", "weight", "acceleration")
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

scatterplotMatrix(~ log(mpg) + displacement + horsepower
                  + weight + acceleration, data=A)

m <- update(m, log(mpg) ~ .)
summary(m)
Anova(m)

m.step <- MASS::stepAIC(m, k=log(nrow(A)), trace=FALSE)
summary(m.step)
Anova(m.step)
plot(predictorEffects(m.step, residuals=TRUE))

cvs <- cvSelect(selectTransStepAIC, data=Auto, seed=76692, model=m.auto,
                predictors=c("displacement", "horsepower", "weight",
                             "acceleration"),
                response="mpg", AIC=FALSE)
cvs
mse(Auto$mpg, exp(fitted(m.step))) # check
cv(m.auto) # pre-transformation & selection
compareFolds(cvs)
