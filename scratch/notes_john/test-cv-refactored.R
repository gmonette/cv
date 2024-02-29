library(cv)
source("~/Documents/R-package-sources/cv/scratch/notes_john/cv-refactored.R")

data("Mroz", package="carData")
m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
old <- cv:::cv.default(m.mroz, criterion=BayesRule, seed=123)
new <- cv.default(m.mroz, criterion=BayesRule, seed=123)
old
new
all.equal(old, new)

old <- cv:::cv.default(m.mroz, criterion=BayesRule, seed=123, ncores=2)
new <- cv.default(m.mroz, criterion=BayesRule, seed=123, ncores=2)
old
new
all.equal(old, new)

old <- cv:::cv.default(m.mroz, criterion=BayesRule, seed=123, reps=2)
new <- cv.default(m.mroz, criterion=BayesRule, seed=123, reps=2)
old
new
all.equal(old, new)


data(Duncan, package="carData")
m <- lm(prestige ~ income + education, data=Duncan)

old <- cv:::cv.lm(m, seed=123)
new <- cv.lm(m, seed=123)
old
new
all.equal(old, new)

old <- cv:::cv.lm(m, seed=123, ncores=2)
new <- cv.lm(m, seed=123, ncores=2)
old
new
all.equal(old, new)

old <- cv:::cv.lm(m, seed=123, reps=2)
new <- cv.lm(m, seed=123, reps=2)
old
new
isTRUE(all.equal(old, new))

microbenchmark::microbenchmark(
  old = cv:::cv.lm(m),
  new = cv.lm(m)
)


data("Mroz", package="carData")
m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
old <- cv:::cv.glm(m.mroz, criterion=BayesRule, seed=123)
new <- cv.glm(m.mroz, criterion=BayesRule, seed=123)
old
new
all.equal(old, new)

old <- cv:::cv.glm(m.mroz, criterion=BayesRule, seed=123,
                   method="Woodbury")
new <- cv.glm(m.mroz, criterion=BayesRule, seed=123,
             method="Woodbury")
old
new
all.equal(old, new)

microbenchmark::microbenchmark(
  old = cv:::cv.glm(m.mroz, criterion=BayesRule, method="Woodbury"),
  new = cv.glm(m.mroz, criterion=BayesRule, method="Woodbury")
)

old <- cv:::cv.glm(m.mroz, criterion=BayesRule, seed=123,
                   method="Woodbury", ncores=2)
new <- cv.glm(m.mroz, criterion=BayesRule, seed=123,
             method="Woodbury", ncores=2)
old
new
all.equal(old, new)

old <- cv:::cv.glm(m.mroz, criterion=BayesRule, seed=123,
                   method="Woodbury", reps=3)
new <- cv.glm(m.mroz, criterion=BayesRule, seed=123,
             method="Woodbury", reps=3)
old
new
all.equal(old, new)

set.seed(123)
X <- matrix(rnorm(1e4*100), 1e4, 100)
y <- X %*% rep(1, 100) + rnorm(1e4)
D <- data.frame(X, y)
m <- lm(y ~ X, data=D)

library(profvis)

profvis({
  old <- cv:::cv.lm(m, method="Woodbury", seed=321)
  new <- cv.lm(m, method="Woodbury", seed=321)
})
all.equal(old, new)

system.time(cv:::cv.lm(m, method="Woodbury", seed=321))
system.time(cv.lm(m, method="Woodbury", seed=321))

profvis({
  old <- cv:::cv.lm(m, method="Woodbury", seed=321, reps=3)
  new <- cv.lm(m, method="Woodbury", seed=321, reps=3)
})
isTRUE(all.equal(old, new))

system.time(cv:::cv.lm(m, method="Woodbury", seed=321, reps=3))
system.time(cv.lm(m, method="Woodbury", seed=321, reps=3))
