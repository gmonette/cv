## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = TRUE,
  warning = TRUE,
  fig.align = "center",
  fig.height = 6,
  fig.width = 7,
  fig.path = "fig/",
  dev = "png",
  comment = "#>" #,
)

# save some typing
knitr::set_alias(w = "fig.width",
                 h = "fig.height",
                 cap = "fig.cap")

# colorize text: use inline as `r colorize(text, color)`
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
      x)
  } else x
}


.opts <- options(digits = 5)

## ----Auto---------------------------------------------------------------------
data("Auto", package="ISLR2")
head(Auto)
dim(Auto)

## ----mpg-horsepower-scatterplot-----------------------------------------------
plot(mpg ~ horsepower, data=Auto)

## ----mpg-horsepower-scatterplot-polynomials-----------------------------------
plot(mpg ~ horsepower, data=Auto)
horsepower <- with(Auto, 
                   seq(min(horsepower), max(horsepower), 
                       length=1000))
for (p in 1:5){
  m <- lm(mpg ~ poly(horsepower,p), data=Auto)
  mpg <- predict(m, newdata=data.frame(horsepower=horsepower))
  lines(horsepower, mpg, col=p + 1, lty=p, lwd=2)
}
legend("topright", legend=1:5, col=2:6, lty=1:5, lwd=2,
       title="Degree", inset=0.02)

## ----mpg-horsepower-MSE-se----------------------------------------------------
library("cv") # for mse() and other functions

var <- mse <- numeric(10)
for (p in 1:10){
  m <- lm(mpg ~ poly(horsepower, p), data=Auto)
  mse[p] <- mse(Auto$mpg, fitted(m))
  var[p] <- summary(m)$sigma^2
}

plot(c(1, 10), range(mse, var), type="n",
     xlab="Degree of polynomial, p",
     ylab="Estimated Squared Error")
lines(1:10, mse, lwd=2, lty=1, col=2, pch=16, type="b")
lines(1:10, var, lwd=2, lty=2, col=3, pch=17, type="b")
legend("topright", inset=0.02,
       legend=c(expression(hat(sigma)^2), "MSE"),
       lwd=2, lty=2:1, col=3:2, pch=17:16)

## ----cv-lm-1------------------------------------------------------------------
m.auto <- lm(mpg ~ poly(horsepower, 2), data=Auto)
summary(m.auto)

cv(m.auto)

## ----cv.lm-2`-----------------------------------------------------------------
cv(m.auto, k="loo")

## ----cv.lm-3------------------------------------------------------------------
cv(m.auto, k="loo", method="naive")

cv(m.auto, k="loo", method="Woodbury")

## ----polyomial-models---------------------------------------------------------
for (p in 1:10){
  command <- paste0("m.", p, "<- lm(mpg ~ poly(horsepower, ", p,
                    "), data=Auto)")
  eval(parse(text=command))
}
objects(pattern="m\\.[0-9]")
summary(m.2) # for example, the quadratic fit

## ----polynomial-regression-CV-------------------------------------------------
# 10-fold CV
cv.auto.10 <- cv(models(m.1, m.2, m.3, m.4, m.5,
                     m.6, m.7, m.8, m.9, m.10),
              data=Auto, seed=2120)
cv.auto.10[1:2] # for the linear and quadratic models

# LOO CV
cv.auto.loo <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 data=Auto, k="loo")
cv.auto.loo[1:2] # linear and quadratic models

## ----polynomial-regression-CV-graph-------------------------------------------
cv.mse.10 <- sapply(cv.auto.10, function(x) x[["adj CV crit"]])
cv.mse.loo <- sapply(cv.auto.loo, function(x) x[["CV crit"]])
plot(c(1, 10), range(cv.mse.10, cv.mse.loo), type="n",
     xlab="Degree of polynomial, p",
     ylab="Cross-Validated MSE")
lines(1:10, cv.mse.10, lwd=2, lty=1, col=2, pch=16, type="b")
lines(1:10, cv.mse.loo, lwd=2, lty=2, col=3, pch=17, type="b")
legend("topright", inset=0.02,
       legend=c("10-Fold CV", "LOO CV"),
       lwd=2, lty=2:1, col=3:2, pch=17:16)

## ----polynomial-regression-CV-graph-2, fig.show="hold"------------------------
plot(cv.auto.10, main="Polynomial Regressions, 10-Fold CV",
     axis.args=list(labels=1:10), xlab="Degree of Polynomial, p")
plot(cv.auto.loo, main="Polynomial Regressions, LOO CV",
     axis.args=list(labels=1:10), xlab="Degree of Polynomial, p")

## ----Mroz-data----------------------------------------------------------------
data("Mroz", package="carData")
head(Mroz, 3)
tail(Mroz, 3)

## ----Mroz-logistic-regresion--------------------------------------------------
m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
summary(m.mroz)

BayesRule(ifelse(Mroz$lfp == "yes", 1, 0), 
          fitted(m.mroz, type="response"))

## ----cv-Mroz-10-fold----------------------------------------------------------
cv(m.mroz, criterion=BayesRule, seed=248)

cv(m.mroz, criterion=BayesRule, seed=248, method="Woodbury")

## ----cv-Mroz-LOO--------------------------------------------------------------
cv(m.mroz, k="loo", criterion=BayesRule)

cv(m.mroz, k="loo", criterion=BayesRule, method="Woodbury")

cv(m.mroz, k="loo", criterion=BayesRule, method="hatvalues")

## ----HSB-data-----------------------------------------------------------------
data("MathAchieve", package="nlme")
dim(MathAchieve)
head(MathAchieve, 3)
tail(MathAchieve, 3)

data("MathAchSchool", package="nlme")
dim(MathAchSchool)
head(MathAchSchool, 2)
tail(MathAchSchool, 2)

## ----data---------------------------------------------------------------------
# Parameters:
set.seed(9693) 
Nb <- 100     # number of groups
Nw <- 5       # number of individuals within groups
Bb <- 0       # between-group regression coefficient on group mean
SDre <- 2.0   # between-group SD of random level relative to group mean of x
SDwithin <- 0.5  # within group SD
Bw <- 1          # within group effect of x
Ay <- 10         # intercept for response
Ax <- 20         # starting level of x
Nx <- Nw*10      # number of distinct x values

Data <- data.frame(
  group = factor(rep(1:Nb, each=Nw)),
  x = Ax + rep(1:Nx, length.out = Nw*Nb)
) |>
  within(
    {
      xm  <- ave(x, group, FUN = mean) # within-group mean
      y <- Ay +
        Bb * xm +                    # contextual effect
        Bw * (x - xm) +              # within-group effect
        rnorm(Nb, sd=SDre)[group] +  # random level by group
        rnorm(Nb*Nw, sd=SDwithin)    # random error within groups
    }
  )

## ----plot1--------------------------------------------------------------------
library("lattice")
library("latticeExtra")
plot <- xyplot(y ~ x, data=Data[1:Nx, ], group=group,
               ylim=c(4, 16),
               par.settings=list(superpose.symbol=list(pch=1, cex=0.7))) +
    layer(panel.ellipse(..., center.cex=0))
plot # display graph

## -----------------------------------------------------------------------------
summary(lm(y ~ x, data=Data))

## ----include=FALSE, echo=FALSE------------------------------------------------
library(lme4) # necessary for some reason to knit vignette in RStudio, harmless otherwise

## -----------------------------------------------------------------------------
# random intercept only:
mod.0 <- lmer(y ~ 1 + (1 | group), Data)
summary(mod.0)

## -----------------------------------------------------------------------------
# effect of x and random intercept:
mod.1 <- lmer(y ~ x + (1 | group), Data)

# effect of x, contextual (student) mean of x, and random intercept:
mod.2 <- lmer(y ~ x + xm + (1 | group), Data)
        # equivalent to y ~ I(x - xm) + xm + (1 | group)

# model generating the data (where Bb = 0)
mod.3 <- lmer(y ~ I(x - xm) + (1 | group), Data)

## -----------------------------------------------------------------------------
Data <- within(Data, {
  fit_mod0.fe <- predict(mod.0, re.form = ~ 0) # fixed effects only
  fit_mod0.re <- predict(mod.0) # fixed and random effects (BLUPs)
  fit_mod1.fe <- predict(mod.1, re.form = ~ 0)
  fit_mod1.re <- predict(mod.1)
  fit_mod2.fe <- predict(mod.2, re.form = ~ 0)
  fit_mod2.re <- predict(mod.2)
  fit_mod3.fe <- predict(mod.3, re.form = ~ 0)
  fit_mod3.re <- predict(mod.3)
})

## -----------------------------------------------------------------------------
Data_long <- reshape(Data[1:Nx, ], direction = "long", sep = ".", 
              timevar = "effect", varying = grep("\\.", names(Data[1:Nx, ])))
Data_long$id <- 1:nrow(Data_long)
Data_long <- reshape(Data_long, direction = "long", sep = "_", 
              timevar = "modelcode",  varying = grep("_", names(Data_long)))
Data_long$model <- factor(
  c("~ 1", "~ 1 + x", "~ 1 + x + xm", "~ 1 + I(x - xm)")
  [match(Data_long$modelcode, c("mod0", "mod1", "mod2", "mod3"))]
)

## ----plot-fits-mod0-----------------------------------------------------------
(plot +
  xyplot(fit ~ x, subset(Data_long, modelcode == "mod0" & effect == "fe"),
         groups=group, type="l", lwd=2) +
  xyplot(fit ~ x, subset(Data_long, modelcode == "mod0" &  effect == "re"),
         groups=group, type="l", lwd=2, lty=3)
) |> update(
  main="Model: y ~ 1 + (1 | group)",
  key=list(
    corner=c(0.05, 0.05),
    text=list(c("fixed effects only","fixed and random")),
    lines=list(lty=c(1, 3))))

## -----------------------------------------------------------------------------
summary(mod.1)

## ----plot-fits-mod1-----------------------------------------------------------
(plot +
  xyplot(fit ~ x, subset(Data_long, modelcode == "mod1" & effect == "fe"),
         groups=group, type="l", lwd=2) +
  xyplot(fit ~ x, subset(Data_long, modelcode == "mod1" & effect == "re"),
         groups=group, type="l", lwd=2, lty=3)
) |> update(
  main="Model: y ~ 1 + x + (1 | group)",
  ylim=c(-15, 35),
  key=list(
    corner=c(0.95, 0.05),
    text=list(c("fixed effects only","fixed and random")),
    lines=list(lty=c(1, 3))))

## -----------------------------------------------------------------------------
summary(mod.2)

## ----plot-fits-mod2-----------------------------------------------------------
(plot +
  xyplot(fit ~ x, subset(Data_long, modelcode == "mod2" & effect == "fe"),
         groups=group, type="l", lwd=2) +
  xyplot(fit ~ x, subset(Data_long, modelcode == "mod2" & effect == "re"),
         groups=group, type="l", lwd=2, lty=3)
) |> update(
  main="Model: y ~ 1 + x + xm + (1 | group)",
  ylim=c(4, 16),
  key=list(
    corner=c(0.05, 0.05),
    text=list(c("fixed effects only","fixed and random")),
    lines=list(lty=c(1, 3))))

## -----------------------------------------------------------------------------
summary(mod.3)

## ----plot-fits-mod3-----------------------------------------------------------
(plot +
  xyplot(fit ~ x, subset(Data_long, modelcode == "mod3" & effect == "fe"),
         groups=group, type="l", lwd=2) +
  xyplot(fit ~ x, subset(Data_long, modelcode == "mod3" & effect == "re"),
         groups=group, type="l", lwd=2, lty=3)
) |> update(
  main="Model: y ~ 1 + I(x - xm) + (1 | group)",
  ylim=c(4, 16),
  key=list(
    corner=c(0.05, 0.05),
    text=list(c("fixed effects only","fixed and random")),
    lines=list(lty=c(1, 3))))

## ----cross-validation-clusters------------------------------------------------
modlist <- models("~ 1"=mod.0, "~ 1 + x"=mod.1, 
                  "~ 1 + x + xm"=mod.2, "~ 1 + I(x - xm)"=mod.3)
cvs_clusters <- cv(modlist, data=Data, cluster="group", k=10, seed=6449)
plot(cvs_clusters, main="Model Comparison, Cluster-Based CV")

## ----cross-validation-cases---------------------------------------------------
cvs_cases <- cv(modlist, data=Data, seed=9693)
plot(cvs_cases, main="Model Comparison, Case-Based CV")

## ----pigs---------------------------------------------------------------------
head(Pigs, 9)
head(xtabs(~ id + week, data=Pigs), 3)
tail(xtabs(~ id + week, data=Pigs), 3)

## ----pigs-graph---------------------------------------------------------------
plot(weight ~ week, data=Pigs, type="n")
for (i in unique(Pigs$id)){
  with(Pigs, lines(x=1:9, y=Pigs[id == i, "weight"],
                   col="gray"))
}
abline(lm(weight ~ week, data=Pigs), col="blue", lwd=2)
lines(with(Pigs, loess.smooth(week, weight, span=0.5)),
      col="magenta", lty=2, lwd=2)

## ----pigs-lmer----------------------------------------------------------------
m.p <- lmer(weight ~ week + (1 | id) + (1 | week),
            data=Pigs, REML=FALSE, # i.e., ML
            control=lmerControl(optimizer="bobyqa"))
summary(m.p)

## ----pigs-cv------------------------------------------------------------------
cv(m.p, clusterVariables="id")

cv(m.p, clusterVariables="week")

cv(m.p, clusterVariables=c("id", "week"), k=10, seed=8469)

## ----pigs-cv-cases------------------------------------------------------------
cv(m.p, k=10, seed=8469)

## ----mroz-reps----------------------------------------------------------------
cv(m.mroz, criterion=BayesRule, seed=248, reps=5, 
   method="Woodbury")

## ----model-comparison-with-reps-----------------------------------------------
cv.auto.reps <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 data=Auto, seed=8004, reps=5)
plot(cv.auto.reps)

## ----generate-selection-data--------------------------------------------------
set.seed(24361) # for reproducibility
D <- data.frame(
  y = rnorm(1000, mean=10),
  X = matrix(rnorm(1000*100), 1000, 100)
)
head(D[, 1:6])

## ----omnibus-F----------------------------------------------------------------
m.full <- lm(y ~ ., data=D)
m.null <- lm(y ~ 1, data=D)
anova(m.null, m.full)

summary(m.null)

## ----forward-selection--------------------------------------------------------
library("MASS")  # for stepAIC()
m.select <- stepAIC(m.null,
                    direction="forward", trace=FALSE,
                    scope=list(lower=~1, upper=formula(m.full)))
summary(m.select)
mse(D$y, fitted(m.select))

## ----cv-selectedModel---------------------------------------------------------
cv(m.select, seed=2529)

## ----compare-selected-models--------------------------------------------------
compareFolds(cv.select)

## ----polynomial-regression-CV-graph-duplicated, echo=FALSE--------------------
cv.mse.10 <- sapply(cv.auto.10, function(x) x[["adj CV crit"]])
cv.mse.loo <- sapply(cv.auto.loo, function(x) x[["CV crit"]])
plot(c(1, 10), range(cv.mse.10, cv.mse.loo), type="n",
     xlab="Degree of polynomial, p",
     ylab="Cross-Validated MSE")
lines(1:10, cv.mse.10, lwd=2, lty=1, col=2, pch=16, type="b")
lines(1:10, cv.mse.loo, lwd=2, lty=2, col=3, pch=17, type="b")
legend("topright", inset=0.02,
       legend=c("10-Fold CV", "LOO CV"),
       lwd=2, lty=2:1, col=3:2, pch=17:16)

## ----nested-CV-polynomials----------------------------------------------------
nestedCV.auto <- cv(selectModelList, Auto,
                    working.model=models(m.1, m.2, m.3, m.4, m.5,
                                        m.6, m.7, m.8, m.9, m.10),
                    save.model=TRUE,
                    seed=2120)
nestedCV.auto
nestedCV.auto$selected.model
cv(m.7, seed=2120) # same seed for same folds

## ----recall-Mroz-regression---------------------------------------------------
summary(m.mroz)

## ----mroz-selection-----------------------------------------------------------
m.mroz.sel <- stepAIC(m.mroz, k=log(nrow(Mroz)),
                      trace=FALSE)
summary(m.mroz.sel)
BayesRule(Mroz$lfp == "yes",
          predict(m.mroz.sel, type="response"))

## ----cv-mroz-regression-------------------------------------------------------
cv(m.mroz.sel, criterion=BayesRule, seed=345266)

## ----cv-mroz-selection--------------------------------------------------------
m.mroz.sel.cv <- cv(selectStepAIC, Mroz, 
                    seed=6681,
                    criterion=BayesRule,
                    working.model=m.mroz,
                    AIC=FALSE)
m.mroz.sel.cv

## ----compare-selected-models-mroz---------------------------------------------
compareFolds(m.mroz.sel.cv)

## ----Prestige-data------------------------------------------------------------
data("Prestige", package="carData")
head(Prestige)
summary(Prestige)

## ----scatterplot-matrix-------------------------------------------------------
library("car")
scatterplotMatrix(~ prestige + income + education + women,
                  data=Prestige, smooth=list(spread=FALSE))

## ----power-transform-Prestige-------------------------------------------------
trans <- powerTransform( cbind(income, education, women) ~ 1,
                         data=Prestige, family="yjPower")
summary(trans)

## ----transformed-predictors---------------------------------------------------
P <- Prestige[, c("prestige", "income", "education", "women")]
(lambdas <- trans$roundlam)
names(lambdas) <- c("income", "education", "women")
for (var in c("income", "education", "women")){
  P[, var] <- yjPower(P[, var], lambda=lambdas[var])
}
summary(P)

scatterplotMatrix(~ prestige + income + education + women,
                  data=P, smooth=list(spread=FALSE))

## ----prestige-regressions-----------------------------------------------------
m.pres <- lm(prestige ~ income + education + women, data=Prestige)
m.pres.trans <- lm(prestige ~ income + education + women, data=P)
mse(Prestige$prestige, fitted(m.pres))
mse(P$prestige, fitted(m.pres.trans))

## ----CR-plots-untransformed---------------------------------------------------
crPlots(m.pres)

## ----CR-plots-transformed-----------------------------------------------------
crPlots(m.pres.trans)

## ----transform-response-------------------------------------------------------
summary(powerTransform(m.pres.trans))

## ----selectTrans--------------------------------------------------------------
selectTrans(data=Prestige, model=m.pres,
            predictors=c("income", "education", "women"),
            response="prestige", family="yjPower")

## ----cv-select-transformations------------------------------------------------
cvs <- cv(selectTrans, data=Prestige, 
          working.model=m.pres, seed=1463,
          predictors=c("income", "education", "women"),
          response="prestige", family="yjPower")
cvs

cv(m.pres, seed=1463) # untransformed model with same folds

compareFolds(cvs)

## ----Auto-redux---------------------------------------------------------------
summary(Auto)
xtabs(~ year, data=Auto)
xtabs(~ origin, data=Auto)
xtabs(~ cylinders, data=Auto)

## ----Auto-explore-------------------------------------------------------------
Auto$cylinders <- factor(Auto$cylinders,
                         labels=c("3.4", "3.4", "5.6", "5.6", "8"))
Auto$year <- as.factor(Auto$year)
Auto$origin <- factor(Auto$origin,
                      labels=c("America", "Europe", "Japan"))
rownames(Auto) <- make.names(Auto$name, unique=TRUE)
Auto$name <- NULL

scatterplotMatrix(~ mpg + displacement + horsepower + weight + acceleration, 
                  smooth=list(spread=FALSE), data=Auto)

## ----Auto-working-model-------------------------------------------------------
m.auto <- lm(mpg ~ ., data = Auto)
summary(m.auto)

Anova(m.auto)

crPlots(m.auto)

## ----Auto-transform-----------------------------------------------------------
num.predictors <- c("displacement", "horsepower", "weight", "acceleration")
tr.x <- powerTransform(Auto[, num.predictors])
summary(tr.x)

## ----Auto-with-transformed-predictors-----------------------------------------
A <- Auto
powers <- tr.x$roundlam
for (pred in num.predictors){
  A[, pred] <- bcPower(A[, pred], lambda=powers[pred])
}
head(A)

m <- update(m.auto, data=A)

## ----Auto-Box-Cox-------------------------------------------------------------
summary(powerTransform(m))

m <- update(m, log(mpg) ~ .)
summary(m)

Anova(m)

## ----Auto-transformed-scatterplot-matrix--------------------------------------
scatterplotMatrix(~ log(mpg) + displacement + horsepower + weight 
                  + acceleration, 
                  smooth=list(spread=FALSE), data=A)

## ----Auto-CR-plots-transformed------------------------------------------------
crPlots(m)

## -----------------------------------------------------------------------------
m.step <- stepAIC(m, k=log(nrow(A)), trace=FALSE)
summary(m.step)

Anova(m.step)

## ----MSE-whole-selected-model-------------------------------------------------
mse(Auto$mpg, exp(fitted(m.step)))

## ----MSE-working-model--------------------------------------------------------
mse(Auto$mpg, fitted(m.auto))

## ----Auto-median-absolute-error-----------------------------------------------
medAbsErr(Auto$mpg, exp(fitted(m.step)))
medAbsErr(Auto$mpg, fitted(m.auto))

## ----Auto-transform-and-select------------------------------------------------
num.predictors
cvs <- cv(selectTransStepAIC, data=Auto, seed=76692, 
          working.model=m.auto, predictors=num.predictors,
          response="mpg", AIC=FALSE, criterion=medAbsErr)
cvs

compareFolds(cvs)

## ----parallel-computation-----------------------------------------------------
system.time(m.mroz.sel.cv <- cv(selectStepAIC, Mroz,
                          seed=6681,
                          criterion=BayesRule,
                          working.model=m.mroz,
                          AIC=FALSE))

system.time(m.mroz.sel.cv.p <- cv(selectStepAIC, Mroz,
                          seed=6681,
                          criterion=BayesRule,
                          working.model=m.mroz,
                          AIC=FALSE,
                          ncores=2))
all.equal(m.mroz.sel.cv, m.mroz.sel.cv.p)

## ----coda, include = FALSE----------------------------------------------------
options(.opts)

