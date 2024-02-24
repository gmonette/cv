
folds(100, 100)

folds(102, 5)
folds(102, 20)

library(cv)
data("Auto", package="ISLR2")
Auto$origin <- factor(Auto$origin)
Auto$name <- NULL
m.auto <- lm(mpg ~ ., data=Auto)
cv <- cv(m.auto, k="loo")
cv1 <- cv.lm1(m.auto, k="loo")
cv2 <- cv.lm2(m.auto, k="loo")
all.equal(cv, cv1)
all.equal(cv, cv2)

microbenchmark::microbenchmark(
  cv1 = cv.lm1(m.auto, k="loo"),
  cv2 = cv.lm2(m.auto, k="loo")
)

library(cv)
source("~/temp/folds.R")
set.seed(123)
n <- 1e6
p <- 100
X <- matrix(rnorm(n*p), n, p)
y <- X %*% rep(1, p) + rnorm(n)
D <- data.frame(X, y)
m <- lm(y ~ X, data=D)

system.time(cv <- cv(m, k="loo")) # cv package
system.time(cv1 <- cv.lm1(m, k="loo")) # current code
system.time(cv2 <- cv.lm2(m, k="loo")) # uses folds()
all.equal(cv, cv1)
all.equal(cv, cv2)

microbenchmark::microbenchmark(
  cv = cv(m, k="loo"),
  cv1 = cv.lm1(m, k="loo"),
  cv2 = cv.lm2(m, k="loo"),
  times=10
)

set.seed(123)
n <- 1e4
p <- 100
X <- matrix(rnorm(n*p), n, p)
y <- X %*% rep(1, p) + rnorm(n)
D <- data.frame(X, y)
m <- lm(y ~ X, data=D)
getLossFn <- cv:::getLossFn
system.time(cvk1 <- cv.lm1(m, k=10, seed=1234)) # current code
system.time(cvk2 <- cv.lm2(m, k=10, seed=1234)) # uses folds()
all.equal(cvk1, cvk2)


ffs <- folds(102, 5)
ffs
fold(ffs, 1)
fold(ffs, 2)
