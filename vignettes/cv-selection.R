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

## ----generate-selection-data--------------------------------------------------
set.seed(24361) # for reproducibility
D <- data.frame(y = rnorm(1000, mean = 10),
                X = matrix(rnorm(1000 * 100), 1000, 100))
head(D[, 1:6])

## ----omnibus-F----------------------------------------------------------------
m.full <- lm(y ~ ., data = D)
m.null <- lm(y ~ 1, data = D)
anova(m.null, m.full)

summary(m.null)

## ----forward-selection--------------------------------------------------------
library("MASS")  # for stepAIC()
m.select <- stepAIC(
  m.null,
  direction = "forward",
  trace = FALSE,
  scope = list(lower =  ~ 1, upper = formula(m.full))
)
summary(m.select)

library("cv")
mse(D$y, fitted(m.select))

## ----cv-selectedModel---------------------------------------------------------
library("cv")

summary(cv(m.select, seed = 2529))

## ----compare-selected-models--------------------------------------------------
compareFolds(cv.select)

## ----polynomial-regression-CV-graph-duplicated, echo=FALSE--------------------
data("Auto", package="ISLR2")
for (p in 1:10) {
  command <- paste0("m.", p, "<- lm(mpg ~ poly(horsepower, ", p,
                    "), data=Auto)")
  eval(parse(text = command))
}
# 10-fold CV
cv.auto.10 <- cv(
  models(m.1, m.2, m.3, m.4, m.5,
         m.6, m.7, m.8, m.9, m.10),
  data = Auto,
  seed = 2120
)

# LOO CV
cv.auto.loo <- cv(models(m.1, m.2, m.3, m.4, m.5,
                         m.6, m.7, m.8, m.9, m.10),
                  data = Auto,
                  k = "loo")

cv.mse.10 <- as.data.frame(cv.auto.10, 
                           rows="cv",             
                           columns="criteria"
                           )$adjusted.criterion
cv.mse.loo <- as.data.frame(cv.auto.loo, 
                           rows="cv",             
                           columns="criteria"
                           )$criterion
plot(
  c(1, 10),
  range(cv.mse.10, cv.mse.loo),
  type = "n",
  xlab = "Degree of polynomial, p",
  ylab = "Cross-Validated MSE"
)
lines(
  1:10,
  cv.mse.10,
  lwd = 2,
  lty = 1,
  col = 2,
  pch = 16,
  type = "b"
)
lines(
  1:10,
  cv.mse.loo,
  lwd = 2,
  lty = 2,
  col = 3,
  pch = 17,
  type = "b"
)
legend(
  "topright",
  inset = 0.02,
  legend = c("10-Fold CV", "LOO CV"),
  lwd = 2,
  lty = 2:1,
  col = 3:2,
  pch = 17:16
)

## ----recursive-CV-polynomials-------------------------------------------------
recursiveCV.auto <- cv(
  selectModelList,
  Auto,
  working.model = models(m.1, m.2, m.3, m.4, m.5,
                         m.6, m.7, m.8, m.9, m.10),
  save.model = TRUE,
  seed = 2120
)
summary(recursiveCV.auto)
(m.sel <- cvInfo(recursiveCV.auto, "selected model"))
cv(m.sel, seed = 2120) # same seed for same folds

## ----recursive-cv-alt---------------------------------------------------------
recursiveCV.auto.alt <- cv(
  models(m.1, m.2, m.3, m.4, m.5,
         m.6, m.7, m.8, m.9, m.10),
  data = Auto,
  seed = 2120,
  recursive = TRUE,
  save.model = TRUE
)
all.equal(recursiveCV.auto, recursiveCV.auto.alt)

## ----recall-Mroz-regression---------------------------------------------------
data("Mroz", package = "carData")
m.mroz <- glm(lfp ~ ., data = Mroz, family = binomial)
summary(m.mroz)

## ----mroz-selection-----------------------------------------------------------
m.mroz.sel <- stepAIC(m.mroz, k = log(nrow(Mroz)),
                      trace = FALSE)
summary(m.mroz.sel)
BayesRule(Mroz$lfp == "yes",
          predict(m.mroz.sel, type = "response"))

## ----cv-mroz-regression-------------------------------------------------------
summary(cv(m.mroz.sel, criterion = BayesRule, seed = 345266))

## ----cv-mroz-selection--------------------------------------------------------
m.mroz.sel.cv <- cv(
  selectStepAIC,
  Mroz,
  seed = 6681,
  criterion = BayesRule,
  working.model = m.mroz,
  AIC = FALSE
)
summary(m.mroz.sel.cv)

## ----compare-selected-models-mroz---------------------------------------------
compareFolds(m.mroz.sel.cv)

## ----Prestige-data------------------------------------------------------------
data("Prestige", package = "carData")
head(Prestige)
summary(Prestige)

## ----scatterplot-matrix-------------------------------------------------------
library("car")
scatterplotMatrix(
  ~ prestige + income + education + women,
  data = Prestige,
  smooth = list(spread = FALSE)
)

## ----power-transform-Prestige-------------------------------------------------
trans <- powerTransform(cbind(income, education, women) ~ 1,
                        data = Prestige,
                        family = "yjPower")
summary(trans)

## ----transformed-predictors---------------------------------------------------
P <- Prestige[, c("prestige", "income", "education", "women")]
(lambdas <- trans$roundlam)
names(lambdas) <- c("income", "education", "women")
for (var in c("income", "education", "women")) {
  P[, var] <- yjPower(P[, var], lambda = lambdas[var])
}
summary(P)

scatterplotMatrix(
  ~ prestige + income + education + women,
  data = P,
  smooth = list(spread = FALSE)
)

## ----prestige-regressions-----------------------------------------------------
m.pres <- lm(prestige ~ income + education + women, data = Prestige)
m.pres.trans <- lm(prestige ~ income + education + women, data = P)
mse(Prestige$prestige, fitted(m.pres))
mse(P$prestige, fitted(m.pres.trans))

## ----CR-plots-untransformed---------------------------------------------------
crPlots(m.pres)

## ----CR-plots-transformed-----------------------------------------------------
crPlots(m.pres.trans)

## ----transform-response-------------------------------------------------------
summary(powerTransform(m.pres.trans))

## ----selectTrans--------------------------------------------------------------
selectTrans(
  data = Prestige,
  model = m.pres,
  predictors = c("income", "education", "women"),
  response = "prestige",
  family = "yjPower"
)

## ----cv-select-transformations------------------------------------------------
cvs <- cv(
  selectTrans,
  data = Prestige,
  working.model = m.pres,
  seed = 1463,
  predictors = c("income", "education", "women"),
  response = "prestige",
  family = "yjPower"
)
summary(cvs)

cv(m.pres, seed = 1463) # untransformed model with same folds

compareFolds(cvs)

## ----Auto-redux---------------------------------------------------------------
summary(Auto)
xtabs( ~ year, data = Auto)
xtabs( ~ origin, data = Auto)
xtabs( ~ cylinders, data = Auto)

## ----Auto-explore-------------------------------------------------------------
Auto$cylinders <- factor(Auto$cylinders,
                         labels = c("3.4", "3.4", "5.6", "5.6", "8"))
Auto$year <- as.factor(Auto$year)
Auto$origin <- factor(Auto$origin,
                      labels = c("America", "Europe", "Japan"))
rownames(Auto) <- make.names(Auto$name, unique = TRUE)
Auto$name <- NULL

scatterplotMatrix(
  ~ mpg + displacement + horsepower + weight + acceleration,
  smooth = list(spread = FALSE),
  data = Auto
)

## ----Auto-working-model-------------------------------------------------------
m.auto <- lm(mpg ~ ., data = Auto)
summary(m.auto)

Anova(m.auto)

crPlots(m.auto)

## ----Auto-transform-----------------------------------------------------------
num.predictors <-
  c("displacement", "horsepower", "weight", "acceleration")
tr.x <- powerTransform(Auto[, num.predictors])
summary(tr.x)

## ----Auto-with-transformed-predictors-----------------------------------------
A <- Auto
powers <- tr.x$roundlam
for (pred in num.predictors) {
  A[, pred] <- bcPower(A[, pred], lambda = powers[pred])
}
head(A)

m <- update(m.auto, data = A)

## ----Auto-Box-Cox-------------------------------------------------------------
summary(powerTransform(m))

m <- update(m, log(mpg) ~ .)
summary(m)

Anova(m)

## ----Auto-transformed-scatterplot-matrix--------------------------------------
scatterplotMatrix(
  ~ log(mpg) + displacement + horsepower + weight
  + acceleration,
  smooth = list(spread = FALSE),
  data = A
)

## ----Auto-CR-plots-transformed------------------------------------------------
crPlots(m)

## -----------------------------------------------------------------------------
m.step <- stepAIC(m, k=log(nrow(A)), trace=FALSE)
summary(m.step)

Anova(m.step)

## ----MSE-whole-selected-model-------------------------------------------------
mse(Auto$mpg, exp(fitted(m.step)))

## ----MSE-working-model--------------------------------------------------------
mse(Auto$mpg, fitted(m.auto))

## ----Auto-median-absolute-error-----------------------------------------------
medAbsErr(Auto$mpg, exp(fitted(m.step)))
medAbsErr(Auto$mpg, fitted(m.auto))

## ----Auto-transform-and-select------------------------------------------------
num.predictors
cvs <- cv(
  selectTransStepAIC,
  data = Auto,
  seed = 76692,
  working.model = m.auto,
  predictors = num.predictors,
  response = "mpg",
  AIC = FALSE,
  criterion = medAbsErr
)
summary(cvs)

compareFolds(cvs)

## ----coda, include = FALSE----------------------------------------------------
options(.opts)

