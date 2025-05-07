# Replication Script for Fox and Monette,
# cv: An R Package for Cross-Validating Regression Models

# packages used: cv, car, microbenchmark, lme4, glmmTMB, lattice, latticeExtra,
#    MASS, carData, nnet, effects, rsample, modeldata, purrr, ggplot2, caret

# session info at the start of the script:

sessionInfo()

## R version 4.5.0 (2025-04-11)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.4.1
##
## Matrix products: default
## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
##
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
## time zone: America/Toronto
## tzcode source: internal
##
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods
## [7] base
##
## loaded via a namespace (and not attached):
##  [1] digest_0.6.37     R6_2.6.1          fastmap_1.2.0
##  [4] xfun_0.52         cachem_1.1.0      rticles_0.27
##  [7] knitr_1.50        htmltools_0.5.8.1 rmarkdown_2.29
## [10] lifecycle_1.0.4   cli_3.6.5         sass_0.4.10
## [13] jquerylib_0.1.4   compiler_4.5.0    rstudioapi_0.17.1
## [16] tools_4.5.0       evaluate_1.0.3    bslib_0.9.0
## [19] yaml_2.3.10       jsonlite_2.0.0    rlang_1.1.6

# set some options
options(digits=5)
palette(car::carPalette())

# code for Sec. 1: Introduction

library("cv")
methods("cv")

# code for Sec. 2: Preliminary Example: Polynomial regression

  # the Auto Data

data("Auto", package = "ISLR2")
summary(Auto)

  # polynomial regressions fit to the auto data

  # Fig. 1 (a)

plot(mpg ~ horsepower, data=Auto)

horsepower <- with(Auto,
                   seq(min(horsepower), max(horsepower),
                       length=1000))
for (p in 1:5){
  m <- lm(mpg ~ poly(horsepower,p), data=Auto)
  mpg <- predict(m, newdata=data.frame(horsepower=horsepower))
  lines(horsepower, mpg, col=p + 1, lty=p, lwd=4)
}
legend("topright", legend=1:5, col=2:6, lty=1:5, lwd=4,
       title="Degree", inset=0.02)

# (preferable to polynomial regression:)

plot(mpg ~ horsepower, data=Auto, log="xy")
abline(lm(log10(mpg) ~ log10(horsepower), data=Auto), lwd=2)

var <- mse <- numeric(10)
for (p in 1:10){
  m <- lm(mpg ~ poly(horsepower, p), data=Auto)
  mse[p] <- mse(Auto$mpg, fitted(m))
  var[p] <- summary(m)$sigma^2
}

  # Fig. 1 (b)

plot(c(1, 10), range(mse, var), type="n",
     xlab="Degree of polynomial, p",
     ylab="Estimated Squared Error")
lines(1:10, mse, lwd=2, lty=1, col=2, pch=16, type="b")
lines(1:10, var, lwd=2, lty=2, col=3, pch=17, type="b")
legend("topright", inset=0.02,
       legend=c(expression(hat(sigma)^2), "MSE"),
       lwd=2, lty=2:1, col=3:2, pch=17:16)


  # the "lm" methods for cv()
args(cv:::cv.lm)

# to make results below reproducible (hidden in the paper)
set.seed(74595)

library("car") # for brief() and other functions

  # CV for second-degree polynomial

m.auto <- lm(mpg ~ poly(horsepower, 2), data = Auto)
brief(m.auto)

(cv.auto <- cv(m.auto, confint = TRUE))
summary(cv.auto)

summary(cv(m.auto, k = "loo"))
summary(cv(m.auto, k = "loo", method = "naive", confint = TRUE))

  # some relative timings

m.auto.glm <- glm(mpg ~ poly(horsepower, 2), data = Auto)
boot::cv.glm(Auto, m.auto.glm)$delta # MSE, biased-corrected MSE

set.seed(19412)
# set times = 100 for greater precision (and patience)
print(microbenchmark::microbenchmark(
  hatvalues = cv(m.auto, k = "loo"),
  Woodbury = cv(m.auto, k = "loo", method = "Woodbury"),
  naive = cv(m.auto, k = "loo", method = "naive"),
  cv.glm = boot::cv.glm(Auto, m.auto.glm),
  times = 10, unit = "relative"), signif = 3)


  # Sec. 2.1: Comparing competing models

mlist <- vector(10, mode = "list")
for (p in 1:10) mlist[[p]] <- lm(mpg ~ poly(horsepower, p), data = Auto)
names(mlist) <- paste0("m.", 1:10)
mlist[2] # 2nd degree polyomial
mlist <- models(mlist)

      # 10-fold CV

cv.auto.10 <- cv(mlist, data = Auto, seed = 2120)
cv.auto.10[2] # 2nd degree polyomial

      # LOO CV

cv.auto.loo <- cv(mlist, data = Auto, k = "loo")
cv.auto.loo[2] # 2nd degree polyomial

      # Fig. 2 (a)

plot(cv.auto.10, main = "Polynomial Regressions, 10-Fold CV",
     axis.args = list(labels = 1:10), xlab = "Degree of Polynomial, p")

      # Fig. 2 (b)

plot(cv.auto.loo, main = "Polynomial Regressions, LOO CV",
     axis.args = list(labels = 1:10), xlab = "Degree of Polynomial, p")

     # re-run with confidence intervals:

cv.auto.10.ci <- cv(mlist, data = Auto, seed = 2120, confint = TRUE)
cv.auto.loo.ci <- cv(mlist, data = Auto, k = "loo", method="Woodbury",
                     confint = TRUE)
plot(cv.auto.10.ci, main = "Polynomial Regressions, 10-Fold CV",
     axis.args = list(labels = 1:10), xlab = "Degree of Polynomial, p")
plot(cv.auto.loo.ci, main = "Polynomial Regressions, LOO CV",
     axis.args = list(labels = 1:10), xlab = "Degree of Polynomial, p")

# Sec. 3: Cross-validating mixed-effects models

  # Sec. 3.1: High-School and Beyond data

data("MathAchieve", package = "nlme")
dim(MathAchieve)
head(MathAchieve, 3)
tail(MathAchieve, 3)

data("MathAchSchool", package = "nlme")
dim(MathAchSchool)
head(MathAchSchool, 2)
tail(MathAchSchool, 2)

    # data management

HSB <- MathAchieve
HSB <- merge(MathAchSchool[, c("School", "Sector")],
             HSB[, c("School", "SES", "MathAch")], by = "School")
names(HSB) <- tolower(names(HSB))
HSB <- within(HSB, {
  mean.ses <- ave(ses, school)
  cses <- ses - mean.ses
})

    # fit mixed model

library("lme4")
hsb.lmer <- lmer(mathach ~ mean.ses*cses + sector*cses
                   + (cses | school), data = HSB)
S(hsb.lmer, brief = TRUE)

    # cluster-based CV

summary(cv(hsb.lmer, k = 10, clusterVariables = "school", seed = 5240))

    # case-based CV

summary(cv(hsb.lmer, seed = 1575)) # note one convergence failure


  # Sec. 3.2. Example: Contrasting cluster-based and case-based CV

library("glmmTMB")

    # generate simulated data

      # Parameters:

Nb <- 20     # number of patients
Nw <- 5      # number of occasions for each patient
Bb <- 1.0    # between-patient regression coefficient on patient means
Bw <- -0.5   # within-patient effect of x

SD_between <- c(0, 5, 6, 8)               # SD between patients
SD_within <- rep(2.5, length(SD_between)) # SD within patients

Nv <- length(SD_within)       # number of variance profiles
SD_ratio <- paste0('SD ratio = ', SD_between,' / ',SD_within)
SD_ratio <- factor(SD_ratio, levels = SD_ratio)

set.seed(833885)

Data_template <- expand.grid(patient = 1:Nb, obs = 1:Nw) |>
  within({
    xw <- seq(-2, 2, length.out = Nw)[obs]
    x <- patient + xw
    xm  <- ave(x, patient)   # within-patient mean
    # Scaled random error within each SD_ratio_i group
    re_std <- scale(resid(lm(rnorm(Nb*Nw) ~ x)))
    re_between <- ave(re_std, patient)
    re_within <- re_std - re_between
    re_between <- scale(re_between)/sqrt(Nw)
    re_within <- scale(re_within)
  })
Data <- do.call(
  rbind,
  lapply(
    1:Nv,
    function(i) {
      cbind(Data_template, SD_ratio_i = i)
    }
  )
)
Data <- within(
  Data,
  {
    SD_within_ <- SD_within[SD_ratio_i]
    SD_between_ <- SD_between[SD_ratio_i]
    SD_ratio <- SD_ratio[SD_ratio_i]
    y <- 10 +
      Bb * xm +                  # contextual effect
      Bw * (x - xm) +            # within-patient effect
      SD_within_ * re_within +   # within patient random effect
      SD_between_ * re_between   # adjustment to between patient random effect
  }
)


      # plot simulated data sets

library("lattice")
library("latticeExtra")
layer <- latticeExtra::layer # to avert conflicts with ggplot2::layer

      # Fig. 3

xyplot(y ~ x | SD_ratio, data = Data, group = patient,
       ylab = list('symptoms (y)', cex = 0.7),
       xlab = list('dosage (x)', cex = 0.7),
       par.strip.text = list(cex = 0.7),
       scales = list(cex = 0.7, x = list(alternating = FALSE)),
       layout = c(Nv, 1),
       par.settings = list(superpose.symbol = list(pch = 1, cex = 0.7))) +
  layer(panel.ellipse(..., center.pch = 16, center.cex = 0.5,
                      level = 0.5),
        panel.abline(a = 10, b = 1))


      # fit mixed models to simulated data sets

model.formulas <- c(
  ' ~ 1'             =  y ~ 1 + (1 | patient),
  '~ 1 + x'          =  y ~ 1 + x + (1 | patient),
  '~ 1 + x + xm'     =  y ~ 1 + x + xm + (1 | patient),
  '~ 1 + I(x - xm)'  =  y ~ 1 + I(x - xm) + (1 | patient)
)
fits <- lapply(split(Data, ~ SD_ratio),
               function(d) {
                 lapply(model.formulas, function(form) {
                   glmmTMB(form, data = d)
                 })
               })

      # predictions based on fixed effects only and on BLUPs

pred.fixed <- lapply(fits, lapply, predict, re.form = ~0)
pred.BLUPs <- lapply(fits, lapply, predict)

      # prepare data and predictions for plotting

Dataf <- lapply(split(Data, ~ SD_ratio),
                function(d) {
                  lapply(names(model.formulas),
                         function(form) cbind(d, formula = form))
                }) |>
  lapply(function(dlist) do.call(rbind, dlist)) |>
  do.call(rbind, args = _) |>
  within(
    {
      pred.fixed <- unlist(pred.fixed)
      pred.BLUPs <- unlist(pred.BLUPs)
      panel <- factor(formula, levels = c(names(model.formulas), 'data'))
    }
  )
Datap <- reshape(Dataf, direction = "long", sep = '.',
                 varying = c('pred.fixed', 'pred.BLUPs'),
                 timevar = 'prediction_type')
Datap$panel2 <- with(Datap, paste0(formula,":",prediction_type))
flevels <- c(unique(Datap$panel2), 'data')
Datap$panel2 <- factor(Datap$panel2, levels = flevels)
Data$panel <- factor('data', levels = c(names(model.formulas), 'data'))
Data$panel2 <- factor('data', levels = flevels)

levs3 <- levels(Data$panel2)
levs3 <- sub(':fixed','', levs3)
levs3 <- sub(':BLUPs',' ', levs3)
Data$panel3 <- factor(sub(':fixed','', Data$panel2,), levels = levs3)
Datap$panel3 <- factor(sub(':BLUPs',' ', Data$panel2,), levels = levs3)

      # Fig. 4

{ xyplot(y ~ x | SD_ratio * panel3, Data,
         groups = patient, type = 'n',
         drop.unused.levels = FALSE,
         par.strip.text = list(cex = 0.4),
         scales = list(x = list(alternating = FALSE, cex = 0.7),
                       y = list(alternating =2, cex = 0.5)),
         xlab = list('dosage', cex = 0.7),
         between =
             list(y = c(0,0,0,0.25,0,0,0,.25)),
         ylab = list(
           paste0(
             paste0(rep(' ',18),  collapse = ''),
                  'fixed-effects',
             paste0(rep(' ',52), collapse = ''),
             'BLUPs',
             paste0(rep(' ',32), collapse = ''),
             'data'),
           cex = 0.7)) +
    glayer(panel.ellipse(..., center.pch = 16, center.cex = 0.5,
                         level = 0.5),
           panel.abline(a = 10, b = 1)) +
    xyplot(pred  ~ x | SD_ratio * panel2, Datap, type = 'l',
           groups = patient,
           drop.unused.levels = FALSE)
} -> re.plot
useOuterStrips(re.plot)

      # case-based and cluster-based CV for various data/model combinations

model_lists <- lapply(fits, function(fitlist) do.call(models, fitlist))
cvs_cases <-
  lapply(1:Nv,
         function(i){
           cv(model_lists[[i]], k = 10,
              data = split(Data, ~ SD_ratio)[[i]])
         })
cvs_clusters <-
  lapply(1:Nv,
         function(i){
           cv(model_lists[[i]],  k = 10,
              data = split(Data, ~SD_ratio)[[i]],
              clusterVariables = 'patient')
         })

names(cvs_clusters) <- names(cvs_cases) <- SD_ratio
dsummary <- expand.grid(SD_ratio_i = names(cvs_cases), model = names(cvs_cases[[1]]))
dsummary$cases <-
  sapply(1:nrow(dsummary), function(i){
    with(dsummary[i,], cvs_cases[[SD_ratio_i]][[model]][['CV crit']])
  })
dsummary$clusters <-
  sapply(1:nrow(dsummary), function(i){
    with(dsummary[i,], cvs_clusters[[SD_ratio_i]][[model]][['CV crit']])
  })

      # Fig. 5

xyplot(clusters + cases ~ model|SD_ratio_i, dsummary,
       auto.key = list(space = 'top', columns = 2, cex = 0.7),
       type = 'b', ylab = list('CV criterion (MSE)', cex = .7),
       xlab = '',
       par.strip.text = list(cex = 0.7),
       layout= c(Nv, 1),
       par.settings =
         list(superpose.line=list(lty = c(1, 2), lwd = 1),
              superpose.symbol=list(cex = 0.7)),
       scales = list(y = list(log = TRUE, cex = 0.7),
                     x = list(alternating = FALSE, rot = 45, cex = 0.7)))

# Sec.4: Cross-validating model specification

  # Sec. 4.1: Example: Data transformation and predictor selection for the Auto data

    # Auto data

names(Auto)
xtabs(~ year, data = Auto)
xtabs(~ origin, data = Auto)
xtabs(~ cylinders, data = Auto)

    # data management

Auto$cylinders <- factor(Auto$cylinders,
                         labels = c("3-4", "3-4", "5-6", "5-6", "8"))
Auto$year <- as.factor(Auto$year)
Auto$origin <- factor(Auto$origin,
                      labels = c("AMER", "EUR", "JAP"))
rownames(Auto) <- make.names(Auto$name, unique = TRUE)
Auto$name <- NULL

    # scatterplot matrix of numeric variables

      # Fig. 6

scatterplotMatrix(~ mpg + displacement + horsepower + weight
                  + acceleration,
                  smooth = list(spread = FALSE), data = Auto, pch= ".")

    # working model

m.auto <- lm(mpg ~ ., data = Auto)

      # Fig. 7

crPlots(m.auto, pch = ".", ylab = "C+R", las = 2)

    # transform numeric predictors towards multi-normality

num.predictors <- c("displacement", "horsepower", "weight",
                    "acceleration")
tr.x <- powerTransform(Auto[, num.predictors])
summary(tr.x)

A <- Auto
powers <- tr.x$roundlam
for (pred in num.predictors){
  A[, pred] <- bcPower(A[, pred], lambda = powers[pred])
}

    # refit model to transformed data

m <- update(m.auto, data = A)

    # transform the response

summary(powerTransform(m))
m <- update(m, log(mpg) ~ .)

      # scatterplot matrix for transformed data

scatterplotMatrix(~ log(mpg) + displacement + horsepower + weight
                  + acceleration,
                  smooth=list(spread = FALSE), data = A, pch = ".")

      # C+R plots for transformed data

crPlots(m, pch = ".", ylab = "C+R", las = 2)

    # predictor selection for model fit to transformed data

library("MASS")
m.step <- stepAIC(m, k = log(nrow(A)), trace = FALSE)
brief(m.step)

mse(Auto$mpg, exp(fitted(m.step)))

mse(Auto$mpg, fitted(m.auto))

      # comparing on the log-response scale

mse(Auto$mpg, fitted(m.step))

mse(Auto$mpg, log(fitted(m.auto)))

    # apply CV to the data transformation and model selection process

num.predictors
cvs <- cv(selectTransStepAIC, data = Auto, seed = 76692,
          working.model = m.auto, predictors = num.predictors,
          response = "mpg", AIC = FALSE)
summary(cvs)

compareFolds(cvs)

  # Sec. 4.2: Example: Applying meta CV to polynomial regression for the Auto data

metaCV.auto <- cv(selectModelList, data = Auto,
                  working.model = mlist, save.model = TRUE,
                  seed = 2120)
summary(metaCV.auto)
brief(m.sel <- cvInfo(metaCV.auto, "selected model"))
    # CV for selected model
cv(m.sel, seed = 2120)

  # equivalent, using meta = TRUE

summary(cv(mlist, data = Auto, seed = 2120, meta = TRUE,
           save.model = TRUE))

# Sec. 5: Extending the cv package

  # multinomial logistic regression for the BEPS data

data("BEPS", package = "carData")
summary(BEPS)

library("nnet")
m.beps <- multinom(vote ~ age + gender + economic.cond.national
                        + economic.cond.household + Blair + Hague
                        + Kennedy + Europe * political.knowledge,
                   data = BEPS, trace = FALSE)

    # effect plot for the Europe x political knowledge interaction

    # Fig. 8

plot(effects::Effect(c("Europe", "political.knowledge"), m.beps,
            xlevels = list(Europe = 1:11, political.knowledge = 0:3),
            fixed.predictors = list(given.values = c(gendermale = 0.5))),
     lines = list(col = c("blue", "red", "orange")),
     axes = list(x = list(rug = FALSE), y = list(style = "stacked")))

    # cv::BayesRule() for a binary response

BayesRule


    # defining BayesRuleMulti()

head(BEPS$vote)
yhat <- predict(m.beps, type = "class")
head(yhat)

BayesRuleMulti <- function(y, yhat){
  result <- mean(y != yhat)
  attr(result, "casewise loss") <- "y != yhat"
  result
}

BayesRuleMulti(BEPS$vote, yhat)

xtabs(~ vote, data = BEPS) / nrow(BEPS)

    # try it out (produces error)

try(cv(m.beps, seed = 3465, criterion = BayesRuleMulti))

    # define GetResponse() method for "multinom" objects

GetResponse.multinom <- function(model, ...) {
  insight::get_response(model)
}

head(GetResponse(m.beps))

    # try it out (still fails)

try(cv(m.beps, seed = 3465, criterion = BayesRuleMulti))
# traceback() # first remove try()

    # create cv() method for "multinom" object

cv.multinom <- function (model, data, criterion = BayesRuleMulti,
                         k, reps, seed, ...) {
    model <- update(model, trace = FALSE)
    NextMethod(
      type = "class", criterion = criterion,
      criterion.name = deparse(substitute(criterion))
    )
  }

    # works
summary(cv(m.beps, seed = 3465))

# Sec. 6: Comparing cv to other R software for cross-validation

  # example from rsample package vignette (modified to use LOO CV)

library("rsample")

data("attrition", package = "modeldata")
nrow(attrition)

mod_form <- Attrition ~ JobSatisfaction + Gender + MonthlyIncome

rs_obj <- loo_cv(attrition)

    # user-supplied function

holdout_results <- function(splits, ...) {
  mod <- glm(..., data = analysis(splits), family = binomial)
  holdout <- assessment(splits)
  res <- broom::augment(mod, newdata = holdout)
  lvls <- levels(holdout$Attrition)
  predictions <- factor(ifelse(res$.fitted > 0, lvls[2], lvls[1]),
                        levels=lvls)
  res$correct <- predictions == holdout$Attrition
  res
}

library("purrr")

rs_obj$results <- map(rs_obj$splits, holdout_results, mod_form)
rs_obj$accuracy <- map_dbl(rs_obj$results, function(x) mean(x$correct))
1 - summary(rs_obj$accuracy) # error rate is complement of accuracy


    # script for direct computation of LOO CV for the example

mod.attrition <- glm(mod_form, data = attrition, family = binomial)
n <- nrow(attrition)
yhat <- numeric(n)
for (i in 1:n){
  m <- update(mod.attrition, data = attrition[-i, ])
  yhat[i] <- predict(m, newdata = attrition[i, ], type = "response")
}
mean(((attrition$Attrition == "Yes") - round(yhat)) ^ 2)

    # LOO CV for the example using the caret package

library("ggplot2")
library("caret")
train(x = attrition[, c("JobSatisfaction", "Gender", "MonthlyIncome")],
      y = attrition$Attrition, method = "glm",
      trControl = trainControl(method = "LOOCV"))


    # LOO cv for the example using cv::cv()

(cv.attrition <- cv(mod.attrition, k = "loo", criterion = BayesRule))
1 - mean(rs_obj$accuracy)

    # comparative timings for the various computations

set.seed(53437)
# set times = 100 for greater precision (and much greater patience)
print(microbenchmark::microbenchmark(
  cv.hatvalues = cv(mod.attrition, k = "loo", criterion = BayesRule,
                    method = "hatvalues"),
  cv.wood = cv(mod.attrition, k = "loo", criterion = BayesRule,
               method = "Woodbury"),
  cv.exact = cv(mod.attrition, k = "loo", criterion = BayesRule),
  direct = {
    n <- nrow(attrition)
    yhat <- numeric(n)
    for (i in 1:n){
      m <- update(mod.attrition, data = attrition[-i, ])
      yhat[i] <- predict(m, newdata = attrition[i, ], type = "response")
    }
    mean(((attrition$Attrition == "Yes") - round(yhat)) ^ 2)
  },
  rsample = {
    rs_obj$results <- map(rs_obj$splits, holdout_results, mod_form)
    rs_obj$accuracy <- map_dbl(rs_obj$results,
                               function(x) mean(x$correct))
  },
  caret = train(x = attrition[, c("JobSatisfaction", "Gender",
                                  "MonthlyIncome")],
              y = attrition$Attrition,
              method = "glm",
              trControl = trainControl(method = "LOOCV")),
  times = 10, unit = "relative"), signif = 3)

# session info at end of script:

sessionInfo()

## R version 4.5.0 (2025-04-11)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.4.1
##
## Matrix products: default
## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
##
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
## time zone: America/Toronto
## tzcode source: internal
##
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets
## [7] methods   base
##
## other attached packages:
##  [1] caret_7.0-1         ggplot2_3.5.2       purrr_1.0.4
##  [4] rsample_1.3.0       nnet_7.3-20         MASS_7.3-65
##  [7] latticeExtra_0.6-30 lattice_0.22-7      glmmTMB_1.1.11
## [10] lme4_1.1-37         Matrix_1.7-3        car_3.1-3
## [13] carData_3.0-5       cv_2.0.3            doParallel_1.0.17
## [16] iterators_1.0.14    foreach_1.5.2
##
## loaded via a namespace (and not attached):
##  [1] Rdpack_2.6.4         DBI_1.2.3            pROC_1.18.5
##  [4] deldir_2.0-4         sandwich_3.1-1       rlang_1.1.6
##  [7] magrittr_2.0.3       multcomp_1.4-28      furrr_0.3.1
## [10] e1071_1.7-16         compiler_4.5.0       mgcv_1.9-3
## [13] reshape2_1.4.4       png_0.1-8            vctrs_0.6.5
## [16] stringr_1.5.1        pkgconfig_2.0.3      fastmap_1.2.0
## [19] backports_1.5.0      rmarkdown_2.29       prodlim_2025.04.28
## [22] nloptr_2.2.1         xfun_0.52            cachem_1.1.0
## [25] jsonlite_2.0.0       recipes_1.3.0        jpeg_0.1-11
## [28] broom_1.0.8          R6_2.6.1             stringi_1.8.7
## [31] bslib_0.9.0          RColorBrewer_1.1-3   parallelly_1.43.0
## [34] boot_1.3-31          rpart_4.1.24         lubridate_1.9.4
## [37] jquerylib_0.1.4      numDeriv_2016.8-1.1  Rcpp_1.0.14
## [40] knitr_1.50           future.apply_1.11.3  zoo_1.8-14
## [43] timechange_0.3.0     splines_4.5.0        tidyselect_1.2.1
## [46] rstudioapi_0.17.1    effects_4.2-2        abind_1.4-8
## [49] yaml_2.3.10          timeDate_4041.110    TMB_1.9.17
## [52] codetools_0.2-20     listenv_0.9.1        tibble_3.2.1
## [55] plyr_1.8.9           withr_3.0.2          evaluate_1.0.3
## [58] future_1.40.0        survival_3.8-3       proxy_0.4-27
## [61] survey_4.4-2         pillar_1.10.2        stats4_4.5.0
## [64] rticles_0.27         reformulas_0.4.1     insight_1.2.0
## [67] generics_0.1.3       scales_1.4.0         minqa_1.2.8
## [70] globals_0.17.0       class_7.3-23         glue_1.8.0
## [73] tools_4.5.0          interp_1.1-6         data.table_1.17.0
## [76] ModelMetrics_1.2.2.2 gower_1.0.2          mvtnorm_1.3-3
## [79] grid_4.5.0           tidyr_1.3.1          mitools_2.4
## [82] rbibutils_2.3        ipred_0.9-15         colorspace_2.1-1
## [85] nlme_3.1-168         Formula_1.2-5        cli_3.6.5
## [88] lava_1.8.1           dplyr_1.1.4          gtable_0.3.6
## [91] sass_0.4.10          digest_0.6.37        TH.data_1.1-3
## [94] farver_2.1.2         htmltools_0.5.8.1    lifecycle_1.0.4
## [97] hardhat_1.4.1        microbenchmark_1.5.0
