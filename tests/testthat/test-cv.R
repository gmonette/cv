# test that lm and glm methods produce the same results for a linear model

data("Auto", package="ISLR2")
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

# test different algorithms for a Gaussian generalized linear model

test_that("glm hatvalues vs exact loo", {
  expect_equal(cv(glm.fit, k="loo", method="hatvalues")$"CV crit",
               cv(glm.fit, k="loo", method="exact")$"CV crit")
})

test_that("glm Woodbury vs exact loo", {
  expect_equal(cv(glm.fit, k="loo", method="Woodbury")[1:3],
               cv(glm.fit, k="loo", method="exact")[1:3])
})

test_that("glm Woodbury vs exact k=10", {
  expect_equal(cv(glm.fit, k=10, method="Woodbury", seed=123)[1:3],
               cv(glm.fit, k=10, method="exact", seed=123)[1:3])
})

# longer tests

# test that parallel computations work for a linear model

test_that("parallel computations lm loo", {
  skip_on_cran()
  expect_equal(cv(lm.fit, k="loo"),
               cv(lm.fit, k="loo", ncores=2))
})

test_that("parallel computations lm k=10", {
  skip_on_cran()
  expect_equal(cv(lm.fit, k=10, seed=123),
               cv(lm.fit, k=10, seed=123, ncores=2))
})

# test that parallel computations work for a generalized linear model

data("Caravan", package="ISLR2")
assign("Cara", Caravan[1:500, c(1:10, 86)], envir=.GlobalEnv)
m.caravan <- glm(Purchase ~ ., data=Cara, family=binomial)

test_that("parallel computations glm k=10", {
  skip_on_cran()
  expect_equal(cv(m.caravan, k=10, criterion=BayesRule, seed=123),
               cv(m.caravan, k=10, criterion=BayesRule, seed=123,
                  ncores=2))
})

test_that("parallel computations glm loo", {
  skip_on_cran()
  expect_equal(cv(m.caravan, k="loo", criterion=BayesRule, seed=123),
               cv(m.caravan, k="loo", criterion=BayesRule, seed=123,
                  ncores=2))
})

# test vs boot:cv.glm()

if (require(boot)){

  test_that("glm method for linear model matches boot::cv.glm()", {
    skip_on_cran()
    expect_equal(boot::cv.glm(Auto, glm.fit)$delta,
                 as.vector(unlist(cv(glm.fit, k="loo")[1:2])))
  })

  test_that("lm method matches boot::cv.glm()", {
    skip_on_cran()
    expect_equal(boot::cv.glm(Auto, glm.fit)$delta,
                 as.vector(unlist(cv(lm.fit, k="loo", method="Woodbury")[1:2])))
  })

  test_that("glm method for GLM matches boot::cv.glm()", {
    skip_on_cran()
    expect_equal(boot::cv.glm(Cara, m.caravan)$delta,
                 as.vector(unlist(cv(m.caravan, k="loo")[1:2])))
  })
}
