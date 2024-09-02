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

## ----Auto---------------------------------------------------------------------
data("Auto", package="ISLR2")
head(Auto)
dim(Auto)

## ----mpg-horsepower-scatterplot-----------------------------------------------
plot(mpg ~ horsepower, data=Auto)

## ----mpg-horsepower-scatterplot-polynomials-----------------------------------
plot(mpg ~ horsepower, data = Auto)
horsepower <- with(Auto,
                   seq(min(horsepower), max(horsepower),
                       length = 1000))
for (p in 1:5) {
  m <- lm(mpg ~ poly(horsepower, p), data = Auto)
  mpg <- predict(m, newdata = data.frame(horsepower = horsepower))
  lines(horsepower,
        mpg,
        col = p + 1,
        lty = p,
        lwd = 2)
}
legend(
  "topright",
  legend = 1:5,
  col = 2:6,
  lty = 1:5,
  lwd = 2,
  title = "Degree",
  inset = 0.02
)

## ----mpg-horsepower-MSE-se----------------------------------------------------
library("cv") # for mse() and other functions

var <- mse <- numeric(10)
for (p in 1:10) {
  m <- lm(mpg ~ poly(horsepower, p), data = Auto)
  mse[p] <- mse(Auto$mpg, fitted(m))
  var[p] <- summary(m)$sigma ^ 2
}

plot(
  c(1, 10),
  range(mse, var),
  type = "n",
  xlab = "Degree of polynomial, p",
  ylab = "Estimated Squared Error"
)
lines(
  1:10,
  mse,
  lwd = 2,
  lty = 1,
  col = 2,
  pch = 16,
  type = "b"
)
lines(
  1:10,
  var,
  lwd = 2,
  lty = 2,
  col = 3,
  pch = 17,
  type = "b"
)
legend(
  "topright",
  inset = 0.02,
  legend = c(expression(hat(sigma) ^ 2), "MSE"),
  lwd = 2,
  lty = 2:1,
  col = 3:2,
  pch = 17:16
)

## ----cv-lm-1------------------------------------------------------------------
m.auto <- lm(mpg ~ poly(horsepower, 2), data = Auto)
summary(m.auto)

cv(m.auto)

## ----cv-summary---------------------------------------------------------------
summary(cv(m.auto))

## ----cv.lm-2`-----------------------------------------------------------------
summary(cv(m.auto, k = "loo"))

## ----cv.lm-3------------------------------------------------------------------
summary(cv(m.auto, k = "loo", method = "naive"))

summary(cv(m.auto, k = "loo", method = "Woodbury"))

## ----polyomial-models---------------------------------------------------------
for (p in 1:10) {
  command <- paste0("m.", p, "<- lm(mpg ~ poly(horsepower, ", p,
                    "), data=Auto)")
  eval(parse(text = command))
}
objects(pattern = "m\\.[0-9]")
summary(m.2) # for example, the quadratic fit

## ----polynomial-regression-CV-------------------------------------------------
# 10-fold CV
cv.auto.10 <- cv(
  models(m.1, m.2, m.3, m.4, m.5,
         m.6, m.7, m.8, m.9, m.10),
  data = Auto,
  seed = 2120
)
cv.auto.10[1:2] # for the linear and quadratic models
summary(cv.auto.10)

# LOO CV
cv.auto.loo <- cv(models(m.1, m.2, m.3, m.4, m.5,
                         m.6, m.7, m.8, m.9, m.10),
                  data = Auto,
                  k = "loo")
cv.auto.loo[1:2] # linear and quadratic models
summary(cv.auto.loo)

## ----polynomial-regression-CV-graph-------------------------------------------
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

## ----polynomial-regression-CV-graph-2, fig.show="hold"------------------------
plot(cv.auto.10, main="Polynomial Regressions, 10-Fold CV",
     axis.args=list(labels=1:10), xlab="Degree of Polynomial, p")
plot(cv.auto.loo, main="Polynomial Regressions, LOO CV",
     axis.args=list(labels=1:10), xlab="Degree of Polynomial, p")

## ----Mroz-data----------------------------------------------------------------
data("Mroz", package = "carData")
head(Mroz, 3)
tail(Mroz, 3)

## ----Mroz-logistic-regresion--------------------------------------------------
m.mroz <- glm(lfp ~ ., data = Mroz, family = binomial)
summary(m.mroz)

BayesRule(ifelse(Mroz$lfp == "yes", 1, 0),
          fitted(m.mroz, type = "response"))

## ----cv-Mroz-10-fold----------------------------------------------------------
summary(cv(m.mroz, criterion = BayesRule, seed = 248))

summary(cv(m.mroz,
   criterion = BayesRule,
   seed = 248,
   method = "Woodbury"))

## ----cv-Mroz-LOO--------------------------------------------------------------
summary(cv(m.mroz, k = "loo", criterion = BayesRule))

summary(cv(m.mroz,
   k = "loo",
   criterion = BayesRule,
   method = "Woodbury"))

summary(cv(m.mroz,
   k = "loo",
   criterion = BayesRule,
   method = "hatvalues"))

## ----mroz-reps----------------------------------------------------------------
summary(cv(
  m.mroz,
  criterion = BayesRule,
  seed = 248,
  reps = 5,
  method = "Woodbury"
))

## ----model-comparison-with-reps-----------------------------------------------
cv.auto.reps <- cv(
  models(m.1, m.2, m.3, m.4, m.5,
         m.6, m.7, m.8, m.9, m.10),
  data = Auto,
  seed = 8004,
  reps = 5
)
plot(cv.auto.reps)

## ----recall-cv.auto.reps------------------------------------------------------
cv.auto.reps
class(cv.auto.reps)

## ----as.data.frame------------------------------------------------------------
D <- as.data.frame(cv.auto.reps)
dim(D)
class(D)

## -----------------------------------------------------------------------------
head(D)

## -----------------------------------------------------------------------------
D <- as.data.frame(cv.auto.reps, columns="criteria")
head(D)
head(subset(D, fold == 0))

## ----summary.cvDataFrame------------------------------------------------------
summary(D, adjusted.criterion ~ model + rep) # fold "0" only
summary(D, criterion ~ model + rep, 
        include="folds") # mean over folds
summary(D, criterion ~ model + rep, fun=sd, 
        include="folds")

## ----coda, include = FALSE----------------------------------------------------
options(.opts)

