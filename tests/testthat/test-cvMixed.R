## The tests in this file run only if the environment variable
##   RUN_ALL_CV_TESTS is set to true, in which case the tests create
##   a number of objects in the global environment, and the tests are slow.
##   Use, e.g. Sys.setenv(RUN_ALL_CV_TESTS = "true")
##   and Sys.unsetenv("RUN_ALL_CV_TESTS") to set the environment variable.

# test parallel computations for mixed models

if (Sys.getenv("RUN_ALL_CV_TESTS") == "true"){

library("lme4")
data("sleepstudy", package="lme4")
fm1 <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

test_that("parallel computations lmer loo clusters", {
  expect_equal(cv(fm1, clusterVariables="Subject"),
               cv(fm1, clusterVariables="Subject", ncores=2))
})

test_that("parallel computations lmer k-fold clusters", {
  expect_equal(cv(fm1, clusterVariables="Subject", k=5, seed=123),
               cv(fm1, clusterVariables="Subject", k=5, seed=123,
                  ncores=2))
})

test_that("parallel computations lmer k-fold cases", {
  expect_equal(suppressWarnings(cv(fm1, k="loo")),
               cv(fm1, k="loo", ncores=2))
})

test_that("parallel computations lmer k-fold cases", {
  expect_equal(cv(fm1, k=5, seed=321),
               cv(fm1, k=5, ncores=2, seed=321))
})

require("nlme")
data("Orthodont", package="nlme")
fm2 <- nlme::lme(distance ~ age + Sex, data = Orthodont,
           random = ~ 1 | Subject)

test_that("parallel computations lme loo clusters", {
  expect_equal(cv(fm2, clusterVariables="Subject"),
               cv(fm2, clusterVariables="Subject", ncores=2))
})

test_that("parallel computations lme k-fold clusters", {
  expect_equal(cv(fm2, clusterVariables="Subject", k=5, seed=123),
               cv(fm2, clusterVariables="Subject", k=5, seed=123,
                  ncores=2))
})

test_that("parallel computations lme k-fold cases", {
  expect_equal(cv(fm2, k=5, seed=321),
               cv(fm2, k=5, seed=321, ncores=2))
})

test_that("parallel computations lme LOO cases", {
  expect_equal(cv(fm2, k="loo"),
               cv(fm2, k="loo", ncores=2))
})

data("Salamanders", package="glmmTMB")
m1 <- glmmTMB::glmmTMB(count ~ mined + (1|site),
                       zi=~mined,
                       family=poisson, data=Salamanders)

test_that("parallel computations glmmTMB k-fold cases", {
  expect_equal(cv(m1, seed=123, k=5),
               cv(m1, seed=123, k=5, ncores=2))
})

m1p <- update(m1, data=Salamanders[1:100, ])
test_that("parallel computations glmmTMB LOO cases", {
  expect_equal(suppressWarnings(cv(m1p, k="loo")),
               cv(m1p, k="loo", ncores=2))
})

test_that("parallel computations glmmTMB k-fold clusrters", {
  expect_equal(cv(m1, clusterVariables="site", k=5, seed=123),
               cv(m1, clusterVariables="site", k=5, seed=123,
                  ncores=2))
})

test_that("parallel computations glmmTMB LOO clusters", {
  expect_equal(cv(m1, clusterVariables="site"),
               cv(m1, clusterVariables="site", ncores=2))
})

}
