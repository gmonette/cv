## The tests in this file run only if the environment variable
##   RUN_ALL_TESTS is set to true, in which case the tests create
##   three data sets in the global environment, Auto, Cara, and
##   Caravan.

if (Sys.getenv("RUN_ALL_TESTS") == "true"){

# test that different algorithms produce the same results
#  using selectStepAIC() for a lm


data(Auto, package="ISLR2")
m <- lm(mpg ~ . - name - origin, data=Auto)

test_that("cvSelect naive vs Woodbury lm", {
  expect_equal(cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m, method="naive"),
               cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m, method="Woodbury")
  )
})

# test that parallel computations work correctly using selectStepAIC()
#  for a lm

test_that("cvSelect parallel lm", {
  expect_equal(cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m),
               cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m, ncores=2)
               )
})

# test that different algorithms produce the same results
#  using selectStepAIC() for a Gaussian glm

m.glm <- lm(mpg ~ . - name - origin, data=Auto)

test_that("cvSelect exact vs Woodbury glm", {
  expect_equal(cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m.glm, method="exact"),
               cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m.glm, method="Woodbury")
  )
})


# longer tests

test_that("cvSelect exact vs hatvalues loo glm", {
  expect_equal(cvSelect(selectStepAIC, Auto, k="loo",
                        model=m.glm, method="exact"),
               cvSelect(selectStepAIC, Auto, k="loo",
                        model=m.glm, method="hatvalues")
  )
})

test_that("cvSelect Woodbury vs hatvalues loo glm", {
  expect_equal(cvSelect(selectStepAIC, Auto, k="loo",
                        model=m.glm, method="Woodbury"),
               cvSelect(selectStepAIC, Auto, k="loo",
                        model=m.glm, method="hatvalues")
  )
})

# test that parallel computations work correctly using selectStepAIC()
#  for a glm

data("Caravan", package="ISLR2")
assign("Cara", Caravan[1:500, c(1:10, 86)], envir=.GlobalEnv)
m.caravan <- glm(Purchase ~ ., data=Cara, family=binomial)

test_that("cvSelect parallel glm", {
  expect_equal(cvSelect(selectStepAIC, Cara, k=5, seed=123,
                        model=m.caravan, criterion=BayesRule),
               cvSelect(selectStepAIC, Cara, k=5, seed=123,
                        model=m.caravan, ncores=2, criterion=BayesRule)
  )
})

}
