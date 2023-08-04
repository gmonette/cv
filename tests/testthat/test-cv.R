# test that lm and glm methods produce the same results for a linear model

data("Auto", package="ISLR")
lm.fit <- lm(mpg ~ horsepower, data=Auto)
glm.fit <- glm(mpg ~ horsepower,data=Auto)

test_that("lm vs glm loo", {
  expect_equal(cv(glm.fit, k="loo")$"CV crit",
               cv(lm.fit, k="loo")$"CV crit")
})

test_that("lm vs glm k=10", {
  expect_equal(cv(lm.fit, k=10, seed=123)[1:3],
  cv(glm.fit, k=10, seed=123)[1:3])
})

# test different algorithms for a linear model

test_that("lm hatvalues vs naive loo", {
  expect_equal(cv(lm.fit, k="loo", method="hatvalues")$"CV crit",
               cv(lm.fit, k="loo", method="naive")$"CV crit")
})

test_that("lm Woodbury vs naive loo", {
  expect_equal(cv(lm.fit, k="loo", method="Woodbury")[1:3],
               cv(lm.fit, k="loo", method="naive")[1:3])
})

test_that("lm Woodbury vs naive k=10", {
  expect_equal(cv(lm.fit, k=10, method="Woodbury", seed=123)[1:3],
               cv(lm.fit, k=10, method="naive", seed=123)[1:3])
})

# test that parallel computations work for a linear model

test_that("parallel computations lm loo", {
  expect_equal(cv(lm.fit, k="loo"),
  cv(lm.fit, k="loo", parallel=TRUE, ncores=2))
})

test_that("parallel computations lm k=10", {
  expect_equal(cv(lm.fit, k=10, seed=123),
               cv(lm.fit, k=10, seed=123, parallel=TRUE, ncores=2))
})

# test that parallel computations work for a generalized linear model

data("Caravan", package="ISLR")
Caravan <- Caravan[1:500, c(1:10, 86)]
m.caravan <- glm(Purchase ~ ., data=Caravan, family=binomial)

test_that("parallel computations glm loo", {
  expect_equal(cv(m.caravan, k="loo", criterion=BayesRule, seed=123),
               cv(m.caravan, k="loo", criterion=BayesRule, seed=123,
                  parallel=TRUE, ncores=2))
})

test_that("parallel computations glm k=10", {
expect_equal(cv(m.caravan, k=10, criterion=BayesRule, seed=123),
             cv(m.caravan, k=10, criterion=BayesRule, seed=123,
                parallel=TRUE, ncores=2))
})
