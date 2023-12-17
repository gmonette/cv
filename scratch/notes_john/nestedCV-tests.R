library(cv)
library(nestedcv)
source("~/Documents/R-package-sources/cv/scratch/notes_john/nestedCV.R")

data("Duncan", package="carData")
m <- lm(prestige ~ income + education, data=Duncan)
cv(m, k=5, criterion=mse, seed=123)

summary(nestedCV(m, seed=123, k=5, reps=10))

all.equal(
  summary(nestedCV.default(m, seed=123, k=5, reps=10)),
  summary(nestedCV(m, seed=123, k=5, reps=10))
)

summary(nestedCV(m, seed=123, k=10, reps=200))

microbenchmark::microbenchmark(
  update=nestedCV.default(m, seed=123, k=10, reps=200),
  Woodbury=nestedCV(m, seed=123, k=10, reps=200),
  times=5
)

D <- Duncan[-2, ]
rownames(D) <- NULL
mm <- lm(prestige ~ income + education, data=D)
summary(nestedCV(mm, seed=123, k=5, reps=10))

# time-consuming:
set.seed(37427)
n <- 100
p <- 20
reps <- 100
result <- vector(reps, mode="list")
for (i in 1:reps){
  X <- matrix(rnorm(n*p), n, p)
  y <- rnorm(n)
  D <- data.frame(X, y)
  ms <- lm(y ~ . - 1, data=D)
  result[[i]] <- nestedCV(ms)
}

lo <- sapply(result, function(x) x["ci.lower.ncv"])
hi <- sapply(result, function(x) x["ci.upper.ncv"])
wid.ncv <- hi - lo
lo <- sapply(result, function(x) x["ci.lower.cv"])
hi <- sapply(result, function(x) x["ci.upper.cv"])
wid.cv <- hi - lo
hist(wid.ncv/wid.cv) # cf. supplement fig. F3, ~ 1 - 1.75
summary(wid.ncv/wid.cv)
plot(wid.cv, wid.ncv)
abline(ls.line <- lm(wid.ncv ~ wid.cv))
summary(ls.line)

## -------

source("~/temp/nested_cv.R")
source("~/temp/nested_cv_helper.R")
library(cv)

se_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
}

fitter_ols <- function(X, Y, idx = NA, funcs_params = NA) {
  if(sum(is.na(idx)) > 0) {idx <- 1:nrow(X)}
  fit <- lm(Y[idx] ~ X[idx, ])

  fit
}

predictor_ols <- function(fit, X_new, funcs_params = NA) {
  X_new %*% fit$coefficients[-1] + fit$coefficients[1]
}

ols_funs <- list(fitter = fitter_ols,
                 predictor = predictor_ols,
                 loss = se_loss,
                 name = "ols")

set.seed(123)
n <- 100
p <- 20
X <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)
D <- data.frame(X, y)
ms <- lm(y ~ . - 1, data=D)

system.time(summary(nestedCV(ms, seed=432)))
nestedcv::naive_cv(X, y, ols_funs)[1:3]

system.time(print(nested_cv(X, y, ols_funs, reps=200)[1:7]))
