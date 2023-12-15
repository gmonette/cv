data("Duncan", package="carData")
m <- lm(prestige ~ income + education, data=Duncan)
cv(m, k=5, criterion=mse, rep=20)

nestedCV(m, seed=123, k=5, reps=10)

nestedCV(m, seed=123, k=10, reps=200)

D <- Duncan[-2, ]
rownames(D) <- NULL
mm <- lm(prestige ~ income + education, data=D)
nestedCV(mm, seed=123, k=5, reps=10)


# (seed <- sample(1e5, 1))

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
