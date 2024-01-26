## -------------------------------------------------------
## Introduction to the cv (cross-validation) package for R
## John Fox and Georges Monette
## York SCS
## February 2023
## -------------------------------------------------------

options(digits=5, show.signif.stars=FALSE)
palette(car::carPalette())

library("cv")

# the cv() function
methods("cv")

# Example: Polynomial regression for the Auto data
#  (from James, Witten, Hastie, and Tibshirani,
#   An Introduction to Statistical Learning with Applications in R)
data("Auto", package="ISLR2")
summary(Auto)

plot(mpg ~ horsepower, data=Auto)
horsepower <- with(Auto,
                   seq(min(horsepower), max(horsepower),
                       length=1000))
for (p in 1:5){
  m <- lm(mpg ~ poly(horsepower,p), data=Auto)
  mpg <- predict(m, newdata=data.frame(horsepower=horsepower))
  lines(horsepower, mpg, col=p + 1, lty=p, lwd=4)
}
legend("topright", legend=1:5, col=2:6, lty=1:5, lwd=4,
       title="Degree", inset=0.02)

var <- mse <- numeric(10)
for (p in 1:10){
  m <- lm(mpg ~ poly(horsepower, p), data=Auto)
  mse[p] <- mse(Auto$mpg, fitted(m))
  var[p] <- summary(m)$sigma^2
}
plot(c(1, 10), range(mse, var), type="n",
     xlab="Degree of polynomial, p",
     ylab="Estimated Squared Error")
lines(1:10, mse, lwd=2, lty=1, col=2, pch=16, type="b")
lines(1:10, var, lwd=2, lty=2, col=3, pch=17, type="b")
legend("topright", inset=0.02,
       legend=c(expression(hat(sigma)^2), "MSE"),
       lwd=2, lty=2:1, col=3:2, pch=17:16)


args(cv:::cv.lm) # the "lm" method for cv()

  # e.g., quadratic model:
m.auto <- lm(mpg ~ poly(horsepower, 2), data=Auto)
summary(m.auto)

cv(m.auto, confint=TRUE) # default 10-fold CV
cv(m.auto, k="loo") # leave-one-out CV
cv(m.auto, k="loo", method="naive", confint=TRUE)
cv(m.auto, k="loo", method="Woodbury", confint=TRUE)

  # using cv.glm() from the boot package:
m.auto.glm <- glm(mpg ~ poly(horsepower, 2), data=Auto)
boot::cv.glm(Auto, m.auto.glm)$delta

  # comparative timings
microbenchmark::microbenchmark(
  hatvalues = cv(m.auto, k="loo"),
  Woodbury = cv(m.auto, k="loo", method="Woodbury"),
  naive = cv(m.auto, k="loo", method="naive"),
  cv.glm = boot::cv.glm(Auto, m.auto.glm),
  times=10
)

# Unit: milliseconds
# expr           min       lq     mean   median       uq
# hatvalues   1.1096   1.2872   1.6809   1.3638   1.4972
# Woodbury   10.2632  10.3265  10.5335  10.5025  10.7472
# naive     218.7086 219.2291 250.1460 222.5466 292.8995
# cv.glm    387.9930 391.6168 420.6690 395.0409 458.9040
#      max neval cld
#   4.7041    10 a
#  10.8839    10 a
# 301.3010    10  b
# 525.2368    10   c

  # polynomial models of degree 1 to 10:
mlist <- vector(10, mode="list")
for (p in 1:10) mlist[[p]] <- lm(mpg ~ poly(horsepower, p), data = Auto)
names(mlist) <- paste0("m.", 1:10)
mlist[2] # e.g., the quadratic fit

  # 10-fold CV
mlist <- do.call(models, mlist) # create "modList" object
cv.auto.10 <- cv(mlist, data=Auto, seed=2120)
cv.auto.10[2] # e.g., for quadratic model

  # LOO CV
cv.auto.loo <- cv(mlist, data=Auto, k="loo")
cv.auto.loo[2] # e.g., for quadratic model

par(mfrow=c(1, 2))
plot(cv.auto.10, main="Polynomial Regressions, 10-Fold CV",
     axis.args=list(labels=1:10), xlab="Degree of Polynomial, p")
plot(cv.auto.loo, main="Polynomial Regressions, LOO CV",
     axis.args=list(labels=1:10), xlab="Degree of Polynomial, p")

names(Auto) # recall
xtabs(~ year, data=Auto)
xtabs(~ origin, data=Auto)
xtabs(~ cylinders, data=Auto)

Auto$cylinders <- factor(Auto$cylinders,
                         labels=c("3.4", "3.4", "5.6", "5.6", "8"))
Auto$year <- as.factor(Auto$year)
Auto$origin <- factor(Auto$origin,
                      labels=c("America", "Europe", "Japan"))
rownames(Auto) <- make.names(Auto$name, unique=TRUE)
Auto$name <- NULL
library("car")
some(Auto)

scatterplotMatrix(~ mpg + displacement + horsepower + weight + acceleration,
                  smooth=list(spread=FALSE), data=Auto, pch=".")

m.auto <- lm(mpg ~ ., data = Auto)
crPlots(m.auto, pch=".")

  # transform predictors towards multivariate normality:
num.predictors <- c("displacement", "horsepower", "weight",
                    "acceleration")
tr.x <- powerTransform(Auto[, num.predictors])
summary(tr.x)

A <- Auto
powers <- tr.x$roundlam
for (pred in num.predictors){
  A[, pred] <- bcPower(A[, pred], lambda=powers[pred])
}
m <- update(m.auto, data=A)

  # Box-Cox transformation of response:
summary(powerTransform(m))
m <- update(m, log(mpg) ~ .)

scatterplotMatrix(~ log(mpg) + displacement + horsepower + weight
                  + acceleration,
                  smooth=list(spread=FALSE), data=A, pch=".")

crPlots(m, pch=".")

  # stepwise predictor selection using BIC:
library("MASS")
m.step <- stepAIC(m, k=log(nrow(A)), trace=FALSE)
summary(m.step)

mse(Auto$mpg, exp(fitted(m.step)))
mse(Auto$mpg, fitted(m.auto))

  # cross-validation of transformations and predictor selection:
args(cvSelect)
args(selectTransStepAIC)

num.predictors
cvs <- cvSelect(selectTransStepAIC, data=Auto, seed=76692, model=m.auto,
                predictors=num.predictors,
                response="mpg", AIC=FALSE)
cvs

compareFolds(cvs)

# Extending the cv package
#   Example: to multinomial logistic regression

data("BEPS", package="carData")
summary(BEPS)

library("nnet")
m.beps <- multinom(vote ~ age + gender + economic.cond.national +
                       economic.cond.household + Blair + Hague + Kennedy +
                       Europe*political.knowledge, data=BEPS, trace=FALSE)

plot(effects::Effect(c("Europe", "political.knowledge"), m.beps,
            xlevels=list(Europe=1:11, political.knowledge=0:3),
            fixed.predictors=list(given.values=c(gendermale=0.5))),
     lines=list(col=c("blue", "red", "orange")),
     axes=list(x=list(rug=FALSE), y=list(style="stacked")))

BayesRule # for binary (e.g., logistic) regression

head(BEPS$vote)
yhat <- predict(m.beps, type="class")
head(yhat)

  # cost criterion for polytomous response:
BayesRuleMulti <- function(y, yhat){
  result <- mean(y != yhat)
  attr(result, "casewise loss") <- "y != yhat"
  result
}

BayesRuleMulti(BEPS$vote, yhat)
xtabs(~ vote, data=BEPS)/nrow(BEPS)

cv(m.beps, seed=3465, criterion=BayesRuleMulti) # fails!

GetResponse.multinom <- function(model, ...) {
  insight::get_response(model)
}

head(GetResponse(m.beps))

cv(m.beps, seed=3465, criterion=BayesRuleMulti) # still fails!

cv.multinom <- function (model, data, criterion=BayesRuleMulti, k, reps,
                         seed, ...){
  NextMethod(type="class", criterion=criterion)
}

cv(m.beps, seed=3465)
