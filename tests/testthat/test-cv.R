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
