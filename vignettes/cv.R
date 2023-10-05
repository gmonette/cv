## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = TRUE,
  warning = TRUE,
  fig.align = "center",
  fig.height = 6,
  fig.width = 7,
  fig.path = "fig/",
  dev = "png",
  comment = "#>" #,
  # eval = nzchar(Sys.getenv("REBUILD_VIGNETTES"))
)

# save some typing
knitr::set_alias(w = "fig.width",
                 h = "fig.height",
                 cap = "fig.cap")

# colorize text: use inline as `r colorize(text, color)`
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
      x)
  } else x
}


.opts <- options(digits = 5)

## ----Auto---------------------------------------------------------------------
data("Auto", package="ISLR2")
head(Auto)
dim(Auto)

## ----mpg-horsepower-scatterplot-----------------------------------------------
plot(mpg ~ horsepower, data=Auto)

## ----mpg-horsepower-scatterplot-polynomials-----------------------------------
plot(mpg ~ horsepower, data=Auto)
horsepower <- with(Auto, 
                   seq(min(horsepower), max(horsepower), 
                       length=1000))
for (p in 1:5){
  m <- lm(mpg ~ poly(horsepower,p), data=Auto)
  mpg <- predict(m, newdata=data.frame(horsepower=horsepower))
  lines(horsepower, mpg, col=p + 1, lty=p, lwd=2)
}
legend("topright", legend=1:5, col=2:6, lty=1:5, lwd=2,
       title="Degree", inset=0.02)

## ----mpg-horsepower-MSE-se2---------------------------------------------------
library("cv") # for mse() and other functions

se <- mse <- numeric(10)
for (p in 1:10){
  m <- lm(mpg ~ poly(horsepower, p), data=Auto)
  mse[p] <- mse(Auto$mpg, fitted(m))
  se[p] <- summary(m)$sigma
}

plot(c(1, 10), range(mse, se^2), type="n",
     xlab="Degree of polynomial, p",
     ylab="Estimated Squared Error")
lines(1:10, mse, lwd=2, lty=1, col=2, pch=16, type="b")
lines(1:10, se^2, lwd=2, lty=2, col=3, pch=17, type="b")
legend("topright", inset=0.02,
       legend=c(expression(hat(sigma)^2), "MSE"),
       lwd=2, lty=2:1, col=3:2, pch=17:16)

## ----cv-lm-1------------------------------------------------------------------
m.auto <- lm(mpg ~ poly(horsepower, 2), data=Auto)
summary(m.auto)

cv(m.auto)

## ----cv.lm-2`-----------------------------------------------------------------
cv(m.auto, k="loo")

## ----cv.lm-3------------------------------------------------------------------
cv(m.auto, k="loo", method="naive")

cv(m.auto, k="loo", method="Woodbury")

## ----polyomial-models---------------------------------------------------------
for (p in 1:10){
  assign(paste0("m.", p),
         lm(mpg ~ poly(horsepower, p), data=Auto))
}
objects(pattern="m\\.[0-9]")
summary(m.2) # for example, the quadratic fit

## ----polynomial-regression-CV-------------------------------------------------
# 10-fold CV
cv.auto.10 <- cv(models(m.1, m.2, m.3, m.4, m.5,
                     m.6, m.7, m.8, m.9, m.10),
              data=Auto, seed=2120)
cv.auto.10[1:2] # for the linear and quadratic models

# LOO CV
cv.auto.loo <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 data=Auto, k="loo")
cv.auto.loo[1:2] # linear and quadratic models

## ----polynomial-regression-CV-graph-------------------------------------------
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

## ----polynomial-regression-CV-graph-2, fig.show="hold"------------------------
plot(cv.auto.10, main="Polynomial Regressions, 10-Fold CV",
     axis.args=list(labels=1:10), xlab="Degree of Polynomial, p")
plot(cv.auto.loo, main="Polynomial Regressions, LOO CV",
     axis.args=list(labels=1:10), xlab="Degree of Polynomial, p")

## ----Mroz-data----------------------------------------------------------------
data("Mroz", package="carData")
head(Mroz, 3)
tail(Mroz, 3)

## ----Mroz-logistic-regresion--------------------------------------------------
m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
summary(m.mroz)

BayesRule(ifelse(Mroz$lfp == "yes", 1, 0), 
          fitted(m.mroz, type="response"))

## ----cv-Mroz-10-fold----------------------------------------------------------
cv(m.mroz, criterion=BayesRule, seed=248)

cv(m.mroz, criterion=BayesRule, seed=248, method="Woodbury")

## ----cv-Mroz-LOO--------------------------------------------------------------
cv(m.mroz, k="loo", criterion=BayesRule)

cv(m.mroz, k="loo", criterion=BayesRule, method="Woodbury")

cv(m.mroz, k="loo", criterion=BayesRule, method="hatvalues")

## ----HSB-data-----------------------------------------------------------------
data("MathAchieve", package="nlme")
dim(MathAchieve)
head(MathAchieve, 3)
tail(MathAchieve, 3)

data("MathAchSchool", package="nlme")
dim(MathAchSchool)
head(MathAchSchool, 2)
tail(MathAchSchool, 2)

## ----mroz-reps----------------------------------------------------------------
cv(m.mroz, criterion=BayesRule, seed=248, reps=5, 
   method="Woodbury")

## ----model-comparison-with-reps-----------------------------------------------
cv.auto.reps <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 data=Auto, seed=8004, reps=5)
plot(cv.auto.reps)

## ----generate-selection-data--------------------------------------------------
set.seed(24361) # for reproducibility
D <- data.frame(
  y = rnorm(1000, mean=10),
  X = matrix(rnorm(1000*100), 1000, 100)
)
head(D[, 1:6])

## ----omnibus-F----------------------------------------------------------------
m.full <- lm(y ~ ., data=D)
m.null <- lm(y ~ 1, data=D)
anova(m.null, m.full)

summary(m.null)

## ----forward-selection--------------------------------------------------------
m.select <- MASS::stepAIC(m.null,
                          direction="forward", trace=FALSE,
                          scope=list(lower=~1, upper=formula(m.full)))
summary(m.select)
mse(D$y, fitted(m.select))

## ----cv-selectedModel---------------------------------------------------------
cv(m.select, seed=2529)

## ----compare-selected-models--------------------------------------------------
compareFolds(cv.select)

## ----recall-Mroz-regression---------------------------------------------------
summary(m.mroz)

## ----mroz-selection-----------------------------------------------------------
m.mroz.sel <- MASS::stepAIC(m.mroz, k=log(nrow(Mroz)),
                            trace=FALSE)
summary(m.mroz.sel)
BayesRule(Mroz$lfp == "yes",
          predict(m.mroz.sel, type="response"))

## ----cv-mroz-regression-------------------------------------------------------
cv(m.mroz.sel, criterion=BayesRule, seed=345266)

## ----cv-mroz-selection--------------------------------------------------------
m.mroz.sel.cv <- cvSelect(selectStepAIC, Mroz, 
                          seed=6681,
                          criterion=BayesRule,
                          model=m.mroz,
                          AIC=FALSE)
m.mroz.sel.cv

## ----compare-selected-models-mroz---------------------------------------------
compareFolds(m.mroz.sel.cv)

## ----Prestige-data------------------------------------------------------------
data("Prestige", package="carData")
head(Prestige)
summary(Prestige)

## ----scatterplot-matrix-------------------------------------------------------
library("car")
scatterplotMatrix(~ prestige + income + education + women,
                  data=Prestige, smooth=list(spread=FALSE))

## ----power-transform-Prestige-------------------------------------------------
trans <- powerTransform( cbind(income, education, women) ~ 1,
                         data=Prestige, family="yjPower")
summary(trans)

## ----transformed-predictors---------------------------------------------------
P <- Prestige[, c("prestige", "income", "education", "women")]
(lambdas <- trans$roundlam)
names(lambdas) <- c("income", "education", "women")
for (var in c("income", "education", "women")){
  P[, var] <- yjPower(P[, var], lambda=lambdas[var])
}
summary(P)

scatterplotMatrix(~ prestige + income + education + women,
                  data=P, smooth=list(spread=FALSE))

## ----prestige-regressions-----------------------------------------------------
m.pres <- lm(prestige ~ income + education + women, data=Prestige)
m.pres.trans <- lm(prestige ~ income + education + women, data=P)
mse(Prestige$prestige, fitted(m.pres))
mse(P$prestige, fitted(m.pres.trans))

## ----CR-plots-untransformed---------------------------------------------------
crPlots(m.pres)

## ----CR-plots-transformed-----------------------------------------------------
crPlots(m.pres.trans)

## ----transform-response-------------------------------------------------------
summary(powerTransform(m.pres.trans))

## ----selectTrans--------------------------------------------------------------
selectTrans(data=Prestige, model=m.pres,
            predictors=c("income", "education", "women"),
            response="prestige", family="yjPower")

## ----cv-select-transformations------------------------------------------------
cvs <- cvSelect(selectTrans, data=Prestige, model=m.pres, seed=1463,
                predictors=c("income", "education", "women"),
                response="prestige",
                family="yjPower")
cvs

cv(m.pres, seed=1463) # untransformed model with same folds

compareFolds(cvs)

## ----parallel-computation-----------------------------------------------------
system.time(m.mroz.sel.cv <- cvSelect(selectStepAIC, Mroz,
                          seed=6681,
                          criterion=BayesRule,
                          model=m.mroz,
                          AIC=FALSE))

system.time(m.mroz.sel.cv.p <- cvSelect(selectStepAIC, Mroz,
                          seed=6681,
                          criterion=BayesRule,
                          model=m.mroz,
                          AIC=FALSE,
                          ncores=2))
all.equal(m.mroz.sel.cv, m.mroz.sel.cv.p)

