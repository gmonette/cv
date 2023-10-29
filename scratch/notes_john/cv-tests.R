
library(ISLR2)
library(car)
library(cv)

lm.fit <- lm(mpg ~ horsepower, data=Auto)
glm.fit <- glm(mpg ~ horsepower,data=Auto)
compareCoefs(lm.fit, glm.fit)


boot::cv.glm(Auto, glm.fit)$delta
cv(glm.fit, k="loo")
cv(lm.fit, k="loo")
cv(glm.fit, k=10, seed=1234)
cv(lm.fit, k=10, seed=1234)
cv(lm.fit, k=10, seed=1234, parallel=TRUE)

library(microbenchmark)
microbenchmark(boot=boot::cv.glm(Auto, glm.fit), mine.glm=cv(glm.fit, k="loo"),
               mine.lm=cv(lm.fit, k="loo"),
               mine.glm.approx=cv(lm.fit, k="loo", method="hatvalues"),
               parallel.glm=cv(glm.fit, k="loo", parallel=TRUE),
               parallel.lm=cv(lm.fit, k="loo", parallel=TRUE),
               times=10)
               # parallel is slower when n is small

# rank-deficient model

mr <- lm(prestige ~ income + education + I(income + education), data=Duncan)
cv(mr, k="loo")
boot::cv.glm(Duncan, mr)$delta

# a larger problem:

set.seed(4321)
D <- data.frame(x=rnorm(1e4))
D$y <- 1 + 2*D$x + rnorm(1e4)
mod <- glm(y ~ x, data=D) # use glm() for boot.cv.glm()
mod.lm <- lm(y ~ x, data=D)
# n-fold CV:
system.time(print(boot::cv.glm(D, mod)$delta)) # takes awhile
system.time(print(cv(mod, k="loo")))                    # takes awhile
system.time(print(cv(mod, k="loo", parallel=TRUE)))     # possibly much faster (depending on n of cores)
system.time(print(cv(mod.lm, k="loo"))) # much faster
system.time(print(cv(mod, k="loo", method="hatvalues"))) # similar
system.time(print(cv(mod.lm, k="loo", method="hatvalues", parallel=TRUE))) # slower for this problem


# Even bigger problem

set.seed(4321)
D <- data.frame(x=rnorm(1e5))
D$y <- 1 + 2*D$x + rnorm(1e5)
mod.lm <- lm(y ~ x, data=D)
system.time(print(cv(mod.lm, k="loo", method="hatvalues")))
system.time(print(cv(mod.lm, k="loo", method="hatvalues", parallel=TRUE)))



# a logistic regression problem

m.caravan <- glm(Purchase ~ ., data=Caravan, family=binomial)
system.time(print(cv(m.caravan, k=10, criterion=BayesRule, seed=123)))
system.time(print(cv(m.caravan, k=10, criterion=BayesRule, seed=123, parallel=TRUE)))
set.seed(123)
system.time(print(boot::cv.glm(Caravan, m.caravan, K=10,
                               cost=function(y, yhat) mean(y != round(yhat)))$delta))

system.time(print(cv(m.caravan, criterion=BayesRule, k="loo", method="Woodbury")))

m.caravan.1 <- glm(Purchase ~ ., data=Caravan[1:1000, ], family=binomial)
system.time(print(cv(m.caravan.1, k="loo", criterion=BayesRule)))
system.time(print(cv(m.caravan.1, k="loo", criterion=BayesRule, parallel=TRUE)))
system.time(print(cv(m.caravan.1, k="loo", criterion=BayesRule, method="Woodbury")))
system.time(print(boot::cv.glm(Caravan[1:1000, ], m.caravan.1,
                               cost=function(y, yhat) mean(y != round(yhat)))$delta))

# test approximation for GLM

data("Mroz", package="carData")

m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
cv(m.mroz, k="loo", criterion=BayesRule)
boot::cv.glm(Mroz, m.mroz,
             cost=function(y, yhat) mean(y != round(yhat)))$delta
cv(m.mroz, k="loo", criterion=BayesRule, method="Woodbury")
cv(m.mroz, k="loo", criterion=BayesRule, method="hatvalues")

unlist(cv(m.mroz, k="loo", criterion=BayesRule)[1:3]) ==
  unlist(cv(m.mroz, k="loo", criterion=BayesRule, method="Woodbury")[1:3])

unlist(cv(m.mroz, k="loo", criterion=BayesRule)[1:3]) -
  unlist(cv(m.mroz, k="loo", criterion=BayesRule, method="Woodbury")[1:3])

microbenchmark::microbenchmark(
  exact=cv(m.mroz, k="loo", criterion=BayesRule),
  wood=cv(m.mroz, k="loo", criterion=BayesRule, method="Woodbury"),
  hat=cv(m.mroz, k="loo", criterion=BayesRule, method="hatvalues"),
   # boot=boot::cv.glm(Mroz, m.mroz,
   #             cost=function(y, yhat) mean(y != round(yhat))),
  times=10)

(cv.e <- cv(m.mroz, criterion=BayesRule, k=10, seed=248))
(cv.a <- cv(m.mroz, criterion=BayesRule, k=10, seed=248, method="Woodbury"))
unlist(cv.e[1:3]) - unlist(cv.a[1:3])
unlist(cv.e[1:3]) == unlist(cv.a[1:3])

# mixed model

library(lme4)
  # from ?lmer:
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
cv(fm1, clusterVariables="Subject")
cv(fm1, seed=447)
cv(fm1, k="loo")

cv(fm1, clusterVariables="Subject", k=3, seed=123)

cv(fm1, clusterVariables="Subject", parallel=TRUE)
cv(fm1, clusterVariables="Subject", k=3, seed=123, parallel=TRUE)

gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data =cbpp, family = binomial)
cv(gm1, clusterVariables="herd")
cv(gm1)

#: various algorithms for lm method

data(Duncan, package="carData")
m <- lm(prestige ~ income + education, data=Duncan)
m1 <- update(m, subset = -1)
(Duncan$prestige - predict(m1, newdata=Duncan))[1]
(residuals(m)/(1 - hatvalues(m)))[1]

cv(m, k="loo") # uses method="hatvalues"
cv(m, k="loo", method="hatvalues")
cv(m, k="loo", method="Woodbury")
cv(m, k="loo", method="naive") # uses default method

microbenchmark::microbenchmark(
  hat=cv(m, k="loo"),
  woodbury=cv(m, k="loo", method="Woodbury"),
  naive=cv(m, k="loo", method="naive")
)

mw <- lm(prestige ~ income + education, weights=income, data=Duncan)
mw1 <- update(mw, subset = -1)
(Duncan$prestige - predict(mw1, newdata=Duncan))[1]
(residuals(mw)/(1 - hatvalues(mw)))[1]

cv(mw, k="loo")
cv(mw, k="loo", method="hatvalues")
cv(mw, k="loo", method="Woodbury")
cv(mw, k="loo", method="naive")

# example from ISLR2 Sec. 5.1

library("cv")
data("Auto", package="ISLR2")

se <- mse <- numeric(10)
for (p in 1:10){
  m <- lm(mpg ~ poly(horsepower, p), data=Auto)
  mse[p] <- mse(Auto$mpg, fitted(m))
  se[p] <- summary(m)$sigma
}

plot(c(1, 10), range(mse, se^2), type="n",
     xlab="Degree of Polynomial, p",
     ylab="Estimated Squared Error")
lines(1:10, mse, lwd=2, lty=1, col=2, pch=16, type="b")
lines(1:10, se^2, lwd=2, lty=2, col=3, pch=17, type="b")
legend("topright", inset=0.02,
       legend=c(expression(hat(sigma)^2), "MSE"),
       lwd=2, lty=2:1, col=3:2, pch=17:16)

m.auto <- lm(mpg ~ poly(horsepower, 2), data=Auto)
summary(m.auto)
cv(m.auto)
cv(m.auto, k="loo")
cv(m.auto, k="loo", method="naive")
cv(m.auto, k="loo", method="Woodbury")

microbenchmark::microbenchmark(
  hatvalues = cv(m.auto, k="loo"),
  Woodbury = cv(m.auto, k="loo", method="Woodbury"),
  naive = cv(m.auto, k="loo", method="naive"),
  times=10
)


cv.loo.mse <- numeric(10)
for (p in 1:10){
  m <- lm(mpg ~ poly(horsepower, p), data=Auto)
  cv.loo.mse[p] <- cv(m, k="loo")$"CV crit"
}

plot(1:10, cv.loo.mse, type="b",
     ylab="Estimated MSE based on LOO CV",
     xlab="Degree of Polynomial, p",
     main="Polynomial Regression of mpg on horsepower")

cv.k.mse <- matrix(0, 10, 100)
set.seed(7951) # for replicability
for (rep in 1:100){
  for (p in 1:10){
    m <- lm(mpg ~ poly(horsepower, p), data=Auto)
    suppressMessages(cv.k.mse[p, rep] <-
                       cv(m, k=10)$"adj CV crit")
  }
}


matplot(1:10, cv.k.mse, type="b", col="gray", lty=2, pch="",
        ylab="Estimated MSE based on 10 reps of 10-fold CV",
        xlab="Degree of Polynomial, p",
        main="Polynomial Regression of mpg on horsepower")
average <- rowMeans(cv.k.mse)
lines(average, lwd=2, col="magenta")

# Caravan ISLR2 example

data("Caravan", package="ISLR2")
Caravan$Purchase <- ifelse(Caravan$Purchase == "Yes", 1, 0)
m.caravan <- glm(Purchase ~ ., data=Caravan, family=binomial)

# system.time(print(cv(m.caravan, k="loo", method="exact")))

cv(m.caravan, k="loo", criterion=BayesRule, method="hatvalues")

cv(m.caravan, k="loo", criterion=BayesRule, method="Woodbury")

cv(m.caravan, k=10, criterion=BayesRule, method="Woodbury", seed=1590)
cv(m.caravan, k=10, criterion=BayesRule, method="exact", seed=1590)

m.caravan.null <- glm(Purchase ~ 1, data=Caravan, family=binomial)
m.caravan.sel <- MASS::stepAIC(m.caravan.null,
              direction="forward",
              scope=list(upper = formula(m.caravan), lower= ~ 1),
              trace=FALSE, k=log(nrow(Caravan)))
summary(m.caravan.sel)
BayesRule(Caravan$Purchase,
          predict(m.caravan.sel, type="response"))

system.time(print(
  cvSelect(selectStepAIC, data=Caravan, seed=9896,
           parallel=TRUE, ncores=2,
           model=m.caravan.null,
           direction="forward", criterion=BayesRule,
           scope=list(upper = formula(m.caravan),
                      lower= ~ 1),
           k.=log(nrow(Caravan)))
))

cvSelect(selectStepAIC, data=Caravan, k=5, seed=9896,
           parallel=TRUE, ncores=5,
           model=m.caravan.null,
           direction="forward", criterion=BayesRule,
           scope=list(upper = formula(m.caravan),
                      lower= ~ 1))

# Auto data

Auto$origin <- factor(Auto$origin,
                      labels=c("America", "Europe", "Japan"))
m.auto <- lm(mpg ~ . - name, data=Auto)
summary(m.auto)
m.auto.sel <- MASS::stepAIC(m.auto, k=log(nrow(Auto)),
                            trace=FALSE)
summary(m.auto.sel)
cv(m.auto.sel, seed=9896)
cvSelect(selectStepAIC, Auto, seed=123,
         model=m.auto, k.=log(nrow(Auto)))

# simulated selection faux pas

set.seed(24361)
D <- data.frame(
  y = rnorm(1000, mean=10),
  X = matrix(rnorm(1000*100), 1000, 100)
)
head(D[, 1:6])

m.full <- lm(y ~ ., data=D)
m.null <- lm(y ~ 1, data=D)
anova(m.null, m.full)
summary(m.null)

m.select <- MASS::stepAIC(m.null, direction="forward", trace=FALSE,
                     scope=list(lower=~1, upper=formula(m.full)))
summary(m.select)

(upper.form <- as.formula(paste("~", paste(paste0("X.", 1:100),
                               collapse=" + "))))

summary(m)

s.m <- MASS::stepAIC(m, direction="forward", trace=FALSE,
                     scope=list(lower=~1, upper=upper.form))
s.m <- MASS::stepAIC(m, direction="forward", trace=FALSE,
               scope=list(lower=~1, upper=upper.form))
summary(s.m)
mse(D$y, fitted(s.m))

cv(s.m)

res <- cvSelect(selectStepAIC, data=D, seed=3791,
         model=m.select,
         direction="forward",
         scope=list(lower=~1, upper=upper.form))


set.seed(123)
y <- rnorm(1000, mean=10)
X <- matrix(rnorm(1000*100), 1000, 100)
colnames(X) <- 1:100
r <- as.vector(cor(X, y))
names(r) <- 1:100
(best5 <- sort(r, decreasing=TRUE)[1:5])

m <- lm(y ~ X[, names(best5)])
summary(m)

D <- data.frame(y, X)

selectBest <- function(D, n.best, y, x){
  y <- D[ , y]
  X <- as.matrix(D[ , x])
  colnames(X) <- 1:ncol(X)
  r <- as.vector(cor(X, y))
  names(r) <- 1:ncol(X)
  (best <- sort(r, decreasing=TRUE)[1:n.best])
  best <- names(best)
  X <- X[, best]
  D <- data.frame(y, X)
  lm(y ~ X, data=D)
  }

selectBest(D, 5, 1, 2:101)

cvSelectBest <- function(data, indices, n.best, y, x){
  if (missing(indices)) {
    model.i <- selectBest(data, n.best=n.best, y=y, x=x)
    fit.o.i <- fitted(model.i)
    return(mse(data[, y], fit.o.i))
  }
  data.minus.fold <- data[-indices, ]
  model.i <- selectBest(data.minus.fold, n.best, y, x)
  X <- cbind(1, as.matrix(data[, names(coef(model.i))[-1]]))
  fit.o.i <- as.vector(X %*% coef(model.i))
  fit.i <- fit.o.i[indices]
  c(mse(data[indices, y], fit.i), mse(data[, y], fit.o.i))
}

cvSelectBest(D, n.best=5, y=1, x=2:101)

cvSelect(cvSelectBest, data=D, seed=321,
         n.best=5, y=1, x=2:101)

data("Mroz", package="carData")
# Mroz$lfp <- as.numeric(Mroz$lfp == "yes")

m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
summary(m.mroz)

m.mroz.sel <- MASS::stepAIC(m.mroz, k=log(nrow(Mroz)),
                            trace=FALSE)
summary(m.mroz.sel)
BayesRule(Mroz$lfp == "yes",
          predict(m.mroz.sel, type="response"))

cv(m.mroz.sel, criterion=BayesRule, seed=345266)

sel <- cvSelect(selectStepAIC, Mroz, seed=6681,
         criterion=BayesRule,
         model=m.mroz, k.=log(nrow(Mroz)))
sel

do.call(car::compareCoefs, sel$models)


## -- mixed models



data("MathAchieve", package="nlme")
data("MathAchSchool", package="nlme")

library("dplyr")
MathAchieve %>% group_by(School) %>%
  summarize(mean.ses = mean(SES)) -> Temp
Temp <- merge(MathAchSchool, Temp, by="School")
HSB <- merge(Temp[, c("School", "Sector", "mean.ses")],
             MathAchieve[, c("School", "SES", "MathAch")], by="School")
names(HSB) <- tolower(names(HSB))

HSB$cses <- with(HSB, ses - mean.ses)

library("lme4")
hsb.lmer.1 <- lmer(mathach ~ mean.ses*cses + sector*cses
                   + (cses | school), data=HSB)
summary(hsb.lmer.1)
cv(hsb.lmer.1, k=5, clusterVariables="school", seed=123)
cv(hsb.lmer.1, k=5, clusterVariables="school", seed=123,
   includeRandom=FALSE)
cv(hsb.lmer.1, k=5, seed=123)

# cv(hsb.lmer.1, clusterVariables="school")

library(lme4)
data("MplsStops", package="carData")
data("MplsDemo", package="carData")
Mpls <- merge(MplsStops, MplsDemo, by="neighborhood")
Mpls <- subset(Mpls,
               subset = problem == "suspicious" & MDC == "MDC",
               select=c("neighborhood", "race", "gender", "black",
                        "personSearch"))
Mpls$race <- car::Recode(Mpls$race,
                    ' "White" = "White"; "Black" = "Black";
        "Native American" = "Native American"; else=NA ')
Mpls$gender[Mpls$gender == "Unknown"] <- NA
Mpls$gender <- droplevels(Mpls$gender)
Mpls <- na.omit(Mpls)
Mpls$neighborhood <- droplevels(Mpls$neighborhood)
Mpls$race <- factor(Mpls$race,
                    levels=c("White", "Black", "Native American"))
nrow(Mpls)
mod.mpls <- glmer(personSearch ~ race*gender + black
                  + (1 | neighborhood), data=Mpls, family=binomial)

# Mpls$y <- with(Mpls, as.numeric(personSearch == "YES"))
# mod.mpls <- glmer(y ~ race*gender + black
#                   + (1 | neighborhood), data=Mpls, family=binomial)

cv(mod.mpls, k=5, clusterVariables="neighborhood", criterion=BayesRule,
   seed=123)
cv(mod.mpls, k=5, clusterVariables="neighborhood", criterion=BayesRule,
   seed=123, includeRandom=FALSE)
cv(mod.mpls, k=5, criterion=BayesRule, seed=123)
