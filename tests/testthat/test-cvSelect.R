# test that different algorithms produce the same results
#  using selectStepAIC()

data(Auto, package="ISLR")
m <- lm(mpg ~ . - name - origin, data=Auto)

test_that("cvSelect naive vs Woodbury", {
  expect_equal(cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m, method="naive"),
               cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m, method="Woodbury")
  )
})

# test that parallel computations work correctly using selectStepAIC()

test_that("cvSelect naive vs Woodbury", {
  expect_equal(cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m),
               cvSelect(selectStepAIC, Auto, k=5, seed=123,
                        model=m, parallel=TRUE, ncores=2)
               )
})
