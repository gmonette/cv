library(cv)
library(lme4)
library(nlme)
source("./__notes_john/cvMixed.R")

example("cv.merMod")

cv(fm1, clusterVariables="Subject")
cv:::cv.merMod(fm1, clusterVariables="Subject")

cv(fm1, clusterVariables="Subject", k=5, seed=123)
cv:::cv.merMod(fm1, clusterVariables="Subject", k=5, seed=123)

cv(fm1, k="loo")
cv:::cv.merMod(fm1, k="loo")

cv(fm1, seed=123)
cv:::cv.merMod(fm1, seed=123)


cv(fm2, clusterVariables="Subject")
cv:::cv.lme(fm2, clusterVariables="Subject")

cv(fm2, clusterVariables="Subject", k=5, seed=123)
cv:::cv.lme(fm2, clusterVariables="Subject", k=5, seed=123)

microbenchmark::microbenchmark(
  new=cv(fm1, clusterVariables="Subject"),
  old=cv:::cv.merMod(fm1, clusterVariables="Subject"),
  times=10
)


library(glmmTMB)

fm1TMB <- glmmTMB(Reaction ~ Days + (Days | Subject),
                  data=sleepstudy)

cv(fm1TMB, clusterVariables="Subject")

cv(fm1TMB, clusterVariables="Subject", k=5, seed=123)

cv(fm1TMB, k="loo")

cv(fm1TMB, seed=123)

