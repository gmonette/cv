library(cv)
#source("~/Documents/R-package-sources/cv/scratch/notes_John/as.data.frame-methods.R")
source(here::here("scratch/notes_John/as.data.frame-methods.R"))

data("Auto", package="ISLR2")
m.auto <- lm(mpg ~ horsepower, data=Auto)
(cv.auto <- cv(m.auto, seed=1234, confint=TRUE))
as.data.frame(cv.auto)

(cv.auto.loo <- cv(m.auto,  k="loo", method="Woodbury", details=TRUE))
D <- as.data.frame(cv.auto.loo)
car::some(D)
head(D)
dim(D)
1 + nrow(Auto)

(cv.auto.reps <- cv(m.auto, seed=1234, reps=3))
D <- as.data.frame(cv.auto.reps)
head(D)
car::some(D)
dim(D)
(10 + 1)*3

data("Duncan", package="carData")
m1 <- lm(prestige ~ income + education, data=Duncan)
m2 <- lm(prestige ~ income + education + type, data=Duncan)
m3 <- lm(prestige ~ (income + education)*type, data=Duncan)
(cv.models <- cv(models(m1=m1, m2=m2, m3=m3),
                 data=Duncan, seed=7949, reps=5))
D <- as.data.frame(cv.models)
dim(D)
(k=10 + 1)*3*5
car::some(D)
head(D, 20)
tail(D, 20)

(cv.models.loo <- cv(models(m1=m1, m2=m2, m3=m3),
                 data=Duncan, k="loo", details=TRUE,
                 method="Woodbury"))
D <- as.data.frame(cv.models.loo)
dim(D)
(n=45 + 1)*3
car::some(D)
head(D)
tail(D)

library("lme4")
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
(cv.clusters <- cv(fm1, clusterVariables="Subject",
                   details=TRUE)) # LOO CV of clusters
as.data.frame(cv.clusters)

(cv.cases <- cv(fm1, seed=447)) # 10-fold CV of cases
D <- as.data.frame(cv.cases) # *** FIXME


data("Auto", package="ISLR2")
m.auto <- lm(mpg ~ . - name - origin, data=Auto)
cv.sel <- cv(selectStepAIC, Auto, seed=123, working.model=m.auto,
             save.model=TRUE)
D <- as.data.frame(cv.sel)
D

cv.sel.reps <- cv(selectStepAIC, Auto, seed=123, working.model=m.auto,
   AIC=FALSE, k=5, reps=3, save.model=TRUE) # via BIC
D <- as.data.frame(cv.sel.reps)
dim(D)
3*(5 + 1)
D

data("Duncan", package="carData")
m1 <- lm(prestige ~ income + education, data=Duncan)
m2 <- lm(prestige ~ income + education + type, data=Duncan)
m3 <- lm(prestige ~ (income + education)*type, data=Duncan)
cv.rec <- cv(selectModelList, data=Duncan, seed=5962,
   working.model=models(m1, m2, m3), save.model=TRUE) # recursive CV
D <- as.data.frame(cv.rec)
D

Auto$year <- as.factor(Auto$year)
Auto$origin <- factor(Auto$origin,
                      labels=c("America", "Europe", "Japan"))
rownames(Auto) <- make.names(Auto$name, unique=TRUE)
Auto$name <- NULL
m.auto <- lm(mpg ~ . , data=Auto)
cvs <- cv(selectTransStepAIC, data=Auto, seed=76692, working.model=m.auto,
          criterion=medAbsErr,
          predictors=c("cylinders", "displacement", "horsepower",
                       "weight", "acceleration"),
          response="mpg", AIC=FALSE, save.model=TRUE)
cvs
D <- as.data.frame(cvs)
D

data("Prestige", package="carData")
m.pres <- lm(prestige ~ income + education + women,
             data=Prestige)
cvt <- cv(selectTrans, data=Prestige, working.model=m.pres, seed=123,
          predictors=c("income", "education", "women"),
          response="prestige", family="yjPower", save.model=TRUE)
D <- as.data.frame(cvt)
D
