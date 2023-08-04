# test that different algorithms produce the same results
#  using selectStepAIC() for a lm

data(Auto, package="ISLR")
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
                        model=m, parallel=TRUE, ncores=2)
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

## FIXME!!! The following works when run directly but fails
## when run in test suite ???

# test that parallel computations work correctly using selectStepAIC()
#  for a glm

# data("Caravan", package="ISLR")
# Caravan <- Caravan[1:500, c(1:10, 86)]
# m.caravan <- glm(Purchase ~ ., data=Caravan, family=binomial)
#
# test_that("cvSelect parallel glm", {
#   expect_equal(cvSelect(selectStepAIC, Caravan, k=5, seed=123,
#          model=m.caravan),
#          cvSelect(selectStepAIC, Caravan, k=5, seed=123,
#          model=m.caravan, parallel=TRUE, ncores=2)
#   )
# })
