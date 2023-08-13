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
  comment = "#>"
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

## ----loadpackages-------------------------------------------------------------
library(cv)    # 

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
m.select <- MASS::stepAIC(m.null, direction="forward", trace=FALSE,
                     scope=list(lower=~1, upper=formula(m.full)))
summary(m.select)
mse(D$y, fitted(m.select))

## ----cv-selectedModel---------------------------------------------------------
cv(m.select, seed=2529)

