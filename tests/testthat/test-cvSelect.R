## The tests in this file run only if the environment variable
##   RUN_ALL_CV_TESTS is set to true, in which case the tests create
##   three data sets in the global environment, Auto, Cara,
##   Caravan, and Prestige.

if (Sys.getenv("RUN_ALL_CV_TESTS") == "true"){

# test that different algorithms produce the same results
#  using selectStepAIC() for a lm


data(Auto, package="ISLR2")
m <- lm(mpg ~ . - name - origin, data=Auto)

test_that("cvSelect naive vs Woodbury lm", {
  expect_equal(cv(selectStepAIC, Auto, k=5, seed=123,
                        working.model=m, method="naive"),
               cv(selectStepAIC, Auto, k=5, seed=123,
                        working.model=m, method="Woodbury")
  )
})

# test that parallel computations work correctly using selectStepAIC()
#  for a lm

test_that("cvSelect parallel lm", {
  expect_equal(cv(selectStepAIC, Auto, k=5, seed=123,
                        working.model=m),
               cv(selectStepAIC, Auto, k=5, seed=123,
                        working.model=m, ncores=2)
               )
})

# test that different algorithms produce the same results
#  using selectStepAIC() for a Gaussian glm

m.glm <- lm(mpg ~ . - name - origin, data=Auto)

test_that("cvSelect exact vs Woodbury glm", {
  expect_equal(cv(selectStepAIC, Auto, k=5, seed=123,
                        working.model=m.glm, method="exact"),
               cv(selectStepAIC, Auto, k=5, seed=123,
                        working.model=m.glm, method="Woodbury")
  )
})


# longer tests

test_that("cvSelect exact vs hatvalues loo glm", {
  expect_equal(cv(selectStepAIC, Auto, k="loo",
                        working.model=m.glm, method="exact"),
               cv(selectStepAIC, Auto, k="loo",
                        working.model=m.glm, method="hatvalues")
  )
})

test_that("cvSelect Woodbury vs hatvalues loo glm", {
  expect_equal(cv(selectStepAIC, Auto, k="loo",
                        working.model=m.glm, method="Woodbury"),
               cv(selectStepAIC, Auto, k="loo",
                        working.model=m.glm, method="hatvalues")
  )
})

# test that parallel computations work correctly using selectStepAIC()
#  for a glm

data("Caravan", package="ISLR2")
assign("Cara", Caravan[1:500, c(1:10, 86)], envir=.GlobalEnv)
m.caravan <- glm(Purchase ~ ., data=Cara, family=binomial)

test_that("cvSelect parallel glm", {
  expect_equal(cv(selectStepAIC, Cara, k=5, seed=123,
                        working.model=m.caravan, criterion=BayesRule),
               cv(selectStepAIC, Cara, k=5, seed=123,
                        working.model=m.caravan, ncores=2, criterion=BayesRule)
  )
})

data("Prestige", package="carData")
m.pres <- lm(prestige ~ income + education + women, data=Prestige)

test_that("cvSelect parallel selectTrans", {
  expect_equal(
    cv(selectTrans, data=Prestige,
       working.model=m.pres, seed=1463,
       predictors=c("income", "education", "women"),
       response="prestige", family="yjPower"),
    cv(selectTrans, data=Prestige,
       working.model=m.pres, seed=1463,
       predictors=c("income", "education", "women"),
       response="prestige", family="yjPower",
       ncores=2)
  )
})

data("Auto", package="ISLR2")
Auto$cylinders <- factor(Auto$cylinders,
                         labels=c("3.4", "3.4", "5.6", "5.6", "8"))
Auto$year <- as.factor(Auto$year)
Auto$origin <- factor(Auto$origin,
                      labels=c("America", "Europe", "Japan"))
rownames(Auto) <- make.names(Auto$name, unique=TRUE)
Auto$name <- NULL
m.auto <- lm(mpg ~ ., data = Auto)
num.predictors <- c("displacement", "horsepower", "weight", "acceleration")

test_that("cvSelect parallel selectTransStepAIC", {
  expect_equal(
    cv(selectTransStepAIC, data=Auto, seed=76692,
       working.model=m.auto, predictors=num.predictors,
       response="mpg", AIC=FALSE, criterion=medAbsErr),
    cv(selectTransStepAIC, data=Auto, seed=76692,
       working.model=m.auto, predictors=num.predictors,
       response="mpg", AIC=FALSE, criterion=medAbsErr,
       ncores=2)
  )
})
}
