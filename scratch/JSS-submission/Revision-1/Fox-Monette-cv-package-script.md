

``` r
# Replication Script for Fox and Monette,
# cv: An R Package for Cross-Validating Regression Models

# packages used: cv, car, microbenchmark, lme4, glmmTMB, lattice, latticeExtra,
#    MASS, carData, nnet, effects, rsample, modeldata, purrr, ggplot2, caret

# session info at the start of the script:

sessionInfo()
```

```
## R version 4.4.1 (2024-06-14)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sonoma 14.6.1
## 
## Matrix products: default
## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
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
## other attached packages:
## [1] knitr_1.48
## 
## loaded via a namespace (and not attached):
## [1] compiler_4.4.1    tools_4.4.1       rstudioapi_0.16.0
## [4] xfun_0.46         evaluate_0.24.0
```

``` r
# set some options
options(digits=5)
palette(car::carPalette())

# code for Sec. 1: Introduction

library("cv")
```

```
## Loading required package: doParallel
```

```
## Loading required package: foreach
```

```
## Loading required package: iterators
```

```
## Loading required package: parallel
```

``` r
methods("cv")
```

```
## [1] cv.default*  cv.function* cv.glm*      cv.glmmTMB* 
## [5] cv.lm*       cv.lme*      cv.merMod*   cv.modList* 
## [9] cv.rlm*     
## see '?methods' for accessing help and source code
```

``` r
# code for Sec. 2: Preliminary Example: Polynomial regression

  # the Auto Data

data("Auto", package = "ISLR2")
summary(Auto)
```

```
##       mpg         cylinders     displacement   horsepower   
##  Min.   : 9.0   Min.   :3.00   Min.   : 68   Min.   : 46.0  
##  1st Qu.:17.0   1st Qu.:4.00   1st Qu.:105   1st Qu.: 75.0  
##  Median :22.8   Median :4.00   Median :151   Median : 93.5  
##  Mean   :23.4   Mean   :5.47   Mean   :194   Mean   :104.5  
##  3rd Qu.:29.0   3rd Qu.:8.00   3rd Qu.:276   3rd Qu.:126.0  
##  Max.   :46.6   Max.   :8.00   Max.   :455   Max.   :230.0  
##                                                             
##      weight      acceleration       year        origin    
##  Min.   :1613   Min.   : 8.0   Min.   :70   Min.   :1.00  
##  1st Qu.:2225   1st Qu.:13.8   1st Qu.:73   1st Qu.:1.00  
##  Median :2804   Median :15.5   Median :76   Median :1.00  
##  Mean   :2978   Mean   :15.5   Mean   :76   Mean   :1.58  
##  3rd Qu.:3615   3rd Qu.:17.0   3rd Qu.:79   3rd Qu.:2.00  
##  Max.   :5140   Max.   :24.8   Max.   :82   Max.   :3.00  
##                                                           
##                  name    
##  amc matador       :  5  
##  ford pinto        :  5  
##  toyota corolla    :  5  
##  amc gremlin       :  4  
##  amc hornet        :  4  
##  chevrolet chevette:  4  
##  (Other)           :365
```

``` r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

``` r
# (preferable to polynomial regression:)

plot(mpg ~ horsepower, data=Auto, log="xy")
abline(lm(log10(mpg) ~ log10(horsepower), data=Auto), lwd=2)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-2.png)

``` r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-3.png)

``` r
  # the "lm" methods for cv()
args(cv:::cv.lm)
```

```
## function (model, data = insight::get_data(model), criterion = mse, 
##     k = 10L, reps = 1L, seed = NULL, details = k <= 10L, confint = n >= 
##         400L, level = 0.95, method = c("auto", "hatvalues", "Woodbury", 
##         "naive"), ncores = 1L, ...) 
## NULL
```

``` r
library("car") # for brief() and other functions
```

```
## Loading required package: carData
```

``` r
  # CV for second-degree polynomial

m.auto <- lm(mpg ~ poly(horsepower, 2), data = Auto)
brief(m.auto)
```

```
##            (Intercept) poly(horsepower, 2)1 poly(horsepower, 2)2
## Estimate        23.446              -120.14                44.09
## Std. Error       0.221                 4.37                 4.37
## 
##  Residual SD = 4.37 on 389 df, R-squared = 0.688
```

``` r
(cv.auto <- cv(m.auto, confint = TRUE))
```

```
## R RNG seed set to 548422
```

```
## cross-validation criterion (mse) = 19.348
```

``` r
summary(cv.auto)
```

```
## 10-Fold Cross Validation
## method: Woodbury
## criterion: mse
## cross-validation criterion = 19.348
## bias-adjusted cross-validation criterion = 19.329
## 95% CI for bias-adjusted CV criterion = (15.838, 22.821)
## full-sample criterion = 18.985
```

``` r
summary(cv(m.auto, k = "loo"))
```

```
## n-Fold Cross Validation
## method: hatvalues
## criterion: mse
## cross-validation criterion = 19.248
```

``` r
summary(cv(m.auto, k = "loo", method = "naive", confint = TRUE))
```

```
## n-Fold Cross Validation
## method: naive
## criterion: mse
## cross-validation criterion = 19.248
## bias-adjusted cross-validation criterion = 19.248
## 95% CI for bias-adjusted CV criterion = (15.779, 22.717)
## full-sample criterion = 18.985
```

``` r
  # some relative timings

m.auto.glm <- glm(mpg ~ poly(horsepower, 2), data = Auto)
boot::cv.glm(Auto, m.auto.glm)$delta # MSE, biased-corrected MSE
```

```
## [1] 19.248 19.248
```

``` r
set.seed(19412)
# set times = 100 for greater precision (and patience)
print(microbenchmark::microbenchmark(
  hatvalues = cv(m.auto, k = "loo"),
  Woodbury = cv(m.auto, k = "loo", method = "Woodbury"),
  naive = cv(m.auto, k = "loo", method = "naive"),
  cv.glm = boot::cv.glm(Auto, m.auto.glm),
  times = 10, unit = "relative"), signif = 3)
```

```
## Warning in microbenchmark::microbenchmark(hatvalues = cv(m.auto,
## k = "loo"), : less accurate nanosecond times to avoid potential
## integer overflows
```

```
## Unit: relative
##       expr   min    lq  mean median     uq max neval cld
##  hatvalues   1.0   1.0   1.0   1.00   1.00   1    10 a  
##   Woodbury  11.3  11.3  10.5   9.97   9.34  12    10 a  
##      naive 192.0 192.0 181.0 169.00 163.00 188    10  b 
##     cv.glm 333.0 328.0 307.0 291.00 289.00 303    10   c
```

``` r
  # Sec. 2.1: Comparing competing models

mlist <- vector(10, mode = "list")
for (p in 1:10) mlist[[p]] <- lm(mpg ~ poly(horsepower, p), data = Auto)
names(mlist) <- paste0("m.", 1:10)
mlist[2] # 2nd degree polyomial
```

```
## $m.2
## 
## Call:
## lm(formula = mpg ~ poly(horsepower, p), data = Auto)
## 
## Coefficients:
##          (Intercept)  poly(horsepower, p)1  poly(horsepower, p)2  
##                 23.4                -120.1                  44.1
```

``` r
mlist <- models(mlist)

      # 10-fold CV

cv.auto.10 <- cv(mlist, data = Auto, seed = 2120)
cv.auto.10[2] # 2nd degree polyomial
```

```
## Model m.2:
## cross-validation criterion = 19.346
```

``` r
      # LOO CV

cv.auto.loo <- cv(mlist, data = Auto, k = "loo")
cv.auto.loo[2] # 2nd degree polyomial
```

```
## Model m.2:
## cross-validation criterion = 19.248
```

``` r
      # Fig. 2 (a)

plot(cv.auto.10, main = "Polynomial Regressions, 10-Fold CV",
     axis.args = list(labels = 1:10), xlab = "Degree of Polynomial, p")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-4.png)

``` r
      # Fig. 2 (b)

plot(cv.auto.loo, main = "Polynomial Regressions, LOO CV",
     axis.args = list(labels = 1:10), xlab = "Degree of Polynomial, p")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-5.png)

``` r
# Sec. 3: Cross-validating mixed-effects models

  # Sec. 3.1: High-School and Beyond data

data("MathAchieve", package = "nlme")
dim(MathAchieve)
```

```
## [1] 7185    6
```

``` r
head(MathAchieve, 3)
```

```
## Grouped Data: MathAch ~ SES | School
##   School Minority    Sex    SES MathAch MEANSES
## 1   1224       No Female -1.528   5.876  -0.428
## 2   1224       No Female -0.588  19.708  -0.428
## 3   1224       No   Male -0.528  20.349  -0.428
```

``` r
tail(MathAchieve, 3)
```

```
## Grouped Data: MathAch ~ SES | School
##      School Minority    Sex    SES MathAch MEANSES
## 7183   9586       No Female  1.332  19.641   0.627
## 7184   9586       No Female -0.008  16.241   0.627
## 7185   9586       No Female  0.792  22.733   0.627
```

``` r
data("MathAchSchool", package = "nlme")
dim(MathAchSchool)
```

```
## [1] 160   7
```

``` r
head(MathAchSchool, 2)
```

```
##      School Size Sector PRACAD DISCLIM HIMINTY MEANSES
## 1224   1224  842 Public   0.35   1.597       0  -0.428
## 1288   1288 1855 Public   0.27   0.174       0   0.128
```

``` r
tail(MathAchSchool, 2)
```

```
##      School Size   Sector PRACAD DISCLIM HIMINTY MEANSES
## 9550   9550 1532   Public   0.45   0.791       0   0.059
## 9586   9586  262 Catholic   1.00  -2.416       0   0.627
```

``` r
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
```

```
## Loading required package: Matrix
```

``` r
hsb.lmer <- lmer(mathach ~ mean.ses*cses + sector*cses
                   + (cses | school), data = HSB)
S(hsb.lmer, brief = TRUE)
```

```
## 
## Estimates of Fixed Effects:
##                     Estimate Std. Error z value Pr(>|z|)
## (Intercept)           12.128      0.199   60.86  < 2e-16
## mean.ses               5.333      0.369   14.45  < 2e-16
## cses                   2.945      0.156   18.93  < 2e-16
## sectorCatholic         1.227      0.306    4.00  6.2e-05
## mean.ses:cses          1.039      0.299    3.48  0.00051
## cses:sectorCatholic   -1.643      0.240   -6.85  7.3e-12
## 
## Estimates of Random Effects (Covariance Components):
##  Groups   Name        Std.Dev. Corr
##  school   (Intercept) 1.543        
##           cses        0.318    0.39
##  Residual             6.060        
## 
## Number of obs: 7185, groups:  school, 160
## 
## logLik     df    AIC    BIC 
## -23252     10  46524  46592
```

``` r
    # cluster-based CV

summary(cv(hsb.lmer, k = 10, clusterVariables = "school", seed = 5240))
```

```
## R RNG seed set to 5240
```

```
## 10-Fold Cross Validation based on 160 {school} clusters
## criterion: mse
## cross-validation criterion = 39.157
## bias-adjusted cross-validation criterion = 39.148
## 95% CI for bias-adjusted CV criterion = (38.066, 40.231)
## full-sample criterion = 39.006
```

``` r
    # case-based CV

summary(cv(hsb.lmer, seed = 1575)) # note one convergence failure
```

```
## R RNG seed set to 1575
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl =
## control$checkConv, : Model failed to converge with max|grad| =
## 0.00587207 (tol = 0.002, component 1)
```

```
## boundary (singular) fit: see help('isSingular')
```

```
## 10-Fold Cross Validation
## criterion: mse
## cross-validation criterion = 37.445
## bias-adjusted cross-validation criterion = 37.338
## 95% CI for bias-adjusted CV criterion = (36.288, 38.388)
## full-sample criterion = 36.068
```

``` r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-6.png)

``` r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-7.png)

``` r
      # case-based and cluster-based CV for various data/model combinations

model_lists <- lapply(fits, function(fitlist) do.call(models, fitlist))
cvs_cases <-
  lapply(1:Nv,
         function(i){
           cv(model_lists[[i]], k = 10,
              data = split(Data, ~ SD_ratio)[[i]])
         })
```

```
## R RNG seed set to 982470
```

```
## R RNG seed set to 903914
```

```
## R RNG seed set to 246168
```

```
## R RNG seed set to 948407
```

``` r
cvs_clusters <-
  lapply(1:Nv,
         function(i){
           cv(model_lists[[i]],  k = 10,
              data = split(Data, ~SD_ratio)[[i]],
              clusterVariables = 'patient')
         })
```

```
## R RNG seed set to 382198
```

```
## R RNG seed set to 533366
```

```
## R RNG seed set to 637459
```

```
## R RNG seed set to 542289
```

``` r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-8.png)

``` r
# Sec.4: Cross-validating model specification

  # Sec. 4.1: Example: Data transformation and predictor selection for the Auto data

    # Auto data

names(Auto)
```

```
## [1] "mpg"          "cylinders"    "displacement" "horsepower"  
## [5] "weight"       "acceleration" "year"         "origin"      
## [9] "name"
```

``` r
xtabs(~ year, data = Auto)
```

```
## year
## 70 71 72 73 74 75 76 77 78 79 80 81 82 
## 29 27 28 40 26 30 34 28 36 29 27 28 30
```

``` r
xtabs(~ origin, data = Auto)
```

```
## origin
##   1   2   3 
## 245  68  79
```

``` r
xtabs(~ cylinders, data = Auto)
```

```
## cylinders
##   3   4   5   6   8 
##   4 199   3  83 103
```

``` r
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
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-9.png)

``` r
    # working model

m.auto <- lm(mpg ~ ., data = Auto)

      # Fig. 7

crPlots(m.auto, pch = ".", ylab = "C+R", las = 2)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-10.png)

``` r
    # transform numeric predictors towards multi-normality

num.predictors <- c("displacement", "horsepower", "weight",
                    "acceleration")
tr.x <- powerTransform(Auto[, num.predictors])
summary(tr.x)
```

```
## bcPower Transformations to Multinormality 
##              Est Power Rounded Pwr Wald Lwr Bnd Wald Upr Bnd
## displacement   -0.0509           0      -0.2082       0.1065
## horsepower     -0.1249           0      -0.2693       0.0194
## weight         -0.0870           0      -0.2948       0.1208
## acceleration    0.3061           0      -0.0255       0.6376
## 
## Likelihood ratio test that transformation parameters are equal to 0
##  (all log transformations)
##                                LRT df  pval
## LR test, lambda = (0 0 0 0) 4.8729  4 0.301
## 
## Likelihood ratio test that no transformations are needed
##                                LRT df   pval
## LR test, lambda = (1 1 1 1) 390.08  4 <2e-16
```

``` r
A <- Auto
powers <- tr.x$roundlam
for (pred in num.predictors){
  A[, pred] <- bcPower(A[, pred], lambda = powers[pred])
}

    # refit model to transformed data

m <- update(m.auto, data = A)

    # transform the response

summary(powerTransform(m))
```

```
## bcPower Transformation to Normality 
##    Est Power Rounded Pwr Wald Lwr Bnd Wald Upr Bnd
## Y1    0.0024           0      -0.1607       0.1654
## 
## Likelihood ratio test that transformation parameter is equal to 0
##  (log transformation)
##                              LRT df  pval
## LR test, lambda = (0) 0.00080154  1 0.977
## 
## Likelihood ratio test that no transformation is needed
##                          LRT df   pval
## LR test, lambda = (1) 124.13  1 <2e-16
```

``` r
m <- update(m, log(mpg) ~ .)

      # Fig. 8

scatterplotMatrix(~ log(mpg) + displacement + horsepower + weight
                  + acceleration,
                  smooth=list(spread = FALSE), data = A, pch = ".")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-11.png)

``` r
      # Fig. 9

crPlots(m.auto, pch = ".", ylab = "C+R", las = 2)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-12.png)

``` r
    # predictor selection for model fit to transformed data

library("MASS")
m.step <- stepAIC(m, k = log(nrow(A)), trace = FALSE)
brief(m.step)
```

```
##            (Intercept) horsepower weight acceleration year71
## Estimate         9.435    -0.2763 -0.609      -0.1314 0.0280
## Std. Error       0.262     0.0561  0.056       0.0532 0.0289
##              year72  year73 year74 year75 year76 year77 year78
## Estimate   -0.00711 -0.0395 0.0528 0.0532 0.0743 0.1379 0.1459
## Std. Error  0.02845  0.0260 0.0300 0.0293 0.0282 0.0289 0.0275
##            year79 year80 year81 year82 originEUR originJAP
## Estimate   0.2360 0.3353 0.2629 0.3234    0.0558    0.0436
## Std. Error 0.0291 0.0311 0.0306 0.0296    0.0168    0.0175
## 
##  Residual SD = 0.105 on 374 df, R-squared = 0.909
```

``` r
mse(Auto$mpg, exp(fitted(m.step)))
```

```
## [1] 6.5121
## attr(,"casewise loss")
## [1] "(y - yhat)^2"
```

``` r
mse(Auto$mpg, fitted(m.auto))
```

```
## [1] 8.0932
## attr(,"casewise loss")
## [1] "(y - yhat)^2"
```

``` r
    # apply CV to the data transformation and model selection process

num.predictors
```

```
## [1] "displacement" "horsepower"   "weight"       "acceleration"
```

``` r
cvs <- cv(selectTransStepAIC, data = Auto, seed = 76692,
          working.model = m.auto, predictors = num.predictors,
          response = "mpg", AIC = FALSE)
```

```
## R RNG seed set to 76692
```

``` r
summary(cvs)
```

```
## 10-Fold Cross Validation
## criterion: mse
## cross-validation criterion = 7.4856
## bias-adjusted cross-validation criterion = 7.3435
## full-sample criterion = 6.5121
```

``` r
compareFolds(cvs)
```

```
## CV criterion by folds:
##  fold.1  fold.2  fold.3  fold.4  fold.5  fold.6  fold.7  fold.8 
##  6.0006 10.5093  7.4444 10.6179  5.0492  7.1832  4.0777 10.5019 
##  fold.9 fold.10 
##  9.1812  4.2507 
## 
## Coefficients by folds:
##         (Intercept) horsepower lam.acceleration lam.displacement
## Fold 1      9.71384   -0.17408          0.50000          0.00000
## Fold 2      9.21713   -0.31480          0.00000          0.00000
## Fold 3      9.61824   -0.19248          0.00000          0.00000
## Fold 4      8.69910   -0.25523          0.50000          0.00000
## Fold 5      9.14403   -0.14934          0.00000          0.00000
## Fold 6      9.63481   -0.16739          0.50000          0.00000
## Fold 7      9.98933   -0.36847          0.00000          0.00000
## Fold 8      9.06301   -0.29721          0.00000          0.00000
## Fold 9      8.88315   -0.22684          0.00000          0.00000
## Fold 10     9.61727   -0.17086          0.00000          0.00000
##         lam.horsepower lam.weight   lambda   weight   year71
## Fold 1         0.00000    0.00000  0.00000 -0.74636  0.03764
## Fold 2         0.00000    0.00000  0.00000 -0.47728  0.02173
## Fold 3         0.00000    0.00000  0.00000 -0.72085  0.01128
## Fold 4         0.00000    0.00000  0.00000 -0.53846  0.02153
## Fold 5         0.00000    0.00000  0.00000 -0.69081  0.02531
## Fold 6         0.00000    0.00000  0.00000 -0.74049  0.02456
## Fold 7        -0.15447    0.00000  0.00000 -0.72843  0.02532
## Fold 8         0.00000    0.00000  0.00000 -0.46392  0.02702
## Fold 9         0.00000    0.00000  0.00000 -0.47136  0.00860
## Fold 10        0.00000    0.00000  0.00000 -0.73550  0.02937
##           year72   year73   year74   year75   year76   year77
## Fold 1  -0.00327 -0.02477  0.05606  0.07080  0.07250  0.14420
## Fold 2  -0.01488 -0.03770  0.04312  0.04031  0.06718  0.13094
## Fold 3  -0.02569 -0.03872  0.05187  0.03837  0.06399  0.11593
## Fold 4  -0.02922 -0.05181  0.04136  0.04072  0.05537  0.12292
## Fold 5  -0.01062 -0.04625  0.05039  0.05596  0.07044  0.13356
## Fold 6   0.00759 -0.03412  0.06266  0.06940  0.07769  0.14211
## Fold 7  -0.01271 -0.04144  0.04568  0.03614  0.07385  0.12976
## Fold 8  -0.02041 -0.05605  0.04437  0.06573  0.08135  0.13158
## Fold 9  -0.03620 -0.04835  0.01906  0.03018  0.05846  0.10536
## Fold 10 -0.00899 -0.03814  0.05408  0.04881  0.07862  0.14101
##           year78   year79   year80   year81   year82 acceleration
## Fold 1   0.14281  0.23266  0.35127  0.25635  0.30546             
## Fold 2   0.14917  0.21871  0.33192  0.26196  0.30943     -0.18909
## Fold 3   0.12601  0.20499  0.32821  0.24478  0.29204             
## Fold 4   0.14083  0.22878  0.32947  0.25140  0.27248     -0.03484
## Fold 5   0.14724  0.24675  0.33331  0.26938  0.32594             
## Fold 6   0.14647  0.23532  0.34761  0.26737  0.33062             
## Fold 7   0.14040  0.23976  0.33998  0.27652  0.30659             
## Fold 8   0.13987  0.23011  0.32880  0.25886  0.30538     -0.17676
## Fold 9   0.11722  0.20665  0.31533  0.23352  0.29375     -0.14514
## Fold 10  0.14313  0.23258  0.35649  0.26214  0.32421             
##         displacement cylinders5-6 cylinders8 originEUR originJAP
## Fold 1                                                          
## Fold 2      -0.09197                                            
## Fold 3                                                          
## Fold 4                   -0.09080   -0.10909                    
## Fold 5                                         0.06261      0.04
## Fold 6                                                          
## Fold 7                                                          
## Fold 8      -0.10542                                            
## Fold 9      -0.13452                                            
## Fold 10
```

``` r
  # Sec. 4.2: Example: Applying recursive CV to polynomial regression for the Auto data

recursiveCV.auto <- cv(selectModelList, data = Auto,
                       working.model = mlist, save.model = TRUE,
                       seed = 2120)
```

```
## R RNG seed set to 2120
```

``` r
summary(recursiveCV.auto)
```

```
## 10-Fold Cross Validation
## criterion: mse
## cross-validation criterion = 20.012
## bias-adjusted cross-validation criterion = 20.619
## full-sample criterion = 18.746
```

``` r
brief(m.sel <- cvInfo(recursiveCV.auto, "selected model"))
```

```
##            (Intercept) poly(horsepower, p)1 poly(horsepower, p)2
## Estimate        23.446               -120.1                 44.1
## Std. Error       0.217                  4.3                  4.3
##            poly(horsepower, p)3 poly(horsepower, p)4
## Estimate                  -3.95                -5.19
## Std. Error                 4.30                 4.30
##            poly(horsepower, p)5 poly(horsepower, p)6
## Estimate                   13.3                -8.55
## Std. Error                  4.3                 4.30
##            poly(horsepower, p)7
## Estimate                   7.98
## Std. Error                 4.30
## 
##  Residual SD = 4.3 on 384 df, R-squared = 0.702
```

``` r
# CV for selected model
cv(m.sel, seed = 2120)
```

```
## R RNG seed set to 2120
```

```
## cross-validation criterion (mse) = 18.898
```

``` r
  # equivalent, using recursive = TRUE

summary(cv(mlist, data = Auto, seed = 2120, recursive = TRUE,
           save.model = TRUE))
```

```
## R RNG seed set to 2120
```

```
## 10-Fold Cross Validation
## cross-validation criterion = 20.012
## bias-adjusted cross-validation criterion = 20.619
## full-sample criterion = 18.746
```

``` r
# Sec. 5: Extending the cv package

  # multinomial logistic regression for the BEPS data

data("BEPS", package = "carData")
summary(BEPS)
```

```
##                vote          age       economic.cond.national
##  Conservative    :462   Min.   :24.0   Min.   :1.00          
##  Labour          :720   1st Qu.:41.0   1st Qu.:3.00          
##  Liberal Democrat:343   Median :53.0   Median :3.00          
##                         Mean   :54.2   Mean   :3.25          
##                         3rd Qu.:67.0   3rd Qu.:4.00          
##                         Max.   :93.0   Max.   :5.00          
##  economic.cond.household     Blair          Hague     
##  Min.   :1.00            Min.   :1.00   Min.   :1.00  
##  1st Qu.:3.00            1st Qu.:2.00   1st Qu.:2.00  
##  Median :3.00            Median :4.00   Median :2.00  
##  Mean   :3.14            Mean   :3.33   Mean   :2.75  
##  3rd Qu.:4.00            3rd Qu.:4.00   3rd Qu.:4.00  
##  Max.   :5.00            Max.   :5.00   Max.   :5.00  
##     Kennedy         Europe      political.knowledge    gender   
##  Min.   :1.00   Min.   : 1.00   Min.   :0.00        female:812  
##  1st Qu.:2.00   1st Qu.: 4.00   1st Qu.:0.00        male  :713  
##  Median :3.00   Median : 6.00   Median :2.00                    
##  Mean   :3.14   Mean   : 6.73   Mean   :1.54                    
##  3rd Qu.:4.00   3rd Qu.:10.00   3rd Qu.:2.00                    
##  Max.   :5.00   Max.   :11.00   Max.   :3.00
```

``` r
library("nnet")
m.beps <- multinom(vote ~ age + gender + economic.cond.national
                        + economic.cond.household + Blair + Hague
                        + Kennedy + Europe * political.knowledge,
                   data = BEPS, trace = FALSE)

    # effect plot for the Europe x political knowledge interaction

    # Fig. 10

plot(effects::Effect(c("Europe", "political.knowledge"), m.beps,
            xlevels = list(Europe = 1:11, political.knowledge = 0:3),
            fixed.predictors = list(given.values = c(gendermale = 0.5))),
     lines = list(col = c("blue", "red", "orange")),
     axes = list(x = list(rug = FALSE), y = list(style = "stacked")))
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-13.png)

``` r
    # cv::BayesRule() for a binary response

BayesRule
```

```
## function (y, yhat) 
## {
##     if (!all(y %in% c(0, 1))) 
##         stop("response values not all 0 or 1")
##     if (any(yhat < 0) || any(yhat > 1)) 
##         stop("fitted values outside of interval [0, 1]")
##     yhat <- round(yhat)
##     result <- mean(y != yhat)
##     attr(result, "casewise loss") <- "y != round(yhat)"
##     result
## }
## <bytecode: 0x12621a928>
## <environment: namespace:cv>
```

``` r
    # defining BayesRuleMulti()

head(BEPS$vote)
```

```
## [1] Liberal Democrat Labour           Labour          
## [4] Labour           Labour           Labour          
## Levels: Conservative Labour Liberal Democrat
```

``` r
yhat <- predict(m.beps, type = "class")
head(yhat)
```

```
## [1] Labour           Labour           Labour          
## [4] Labour           Liberal Democrat Labour          
## Levels: Conservative Labour Liberal Democrat
```

``` r
BayesRuleMulti <- function(y, yhat){
  result <- mean(y != yhat)
  attr(result, "casewise loss") <- "y != yhat"
  result
}

BayesRuleMulti(BEPS$vote, yhat)
```

```
## [1] 0.31869
## attr(,"casewise loss")
## [1] "y != yhat"
```

``` r
xtabs(~ vote, data = BEPS) / nrow(BEPS)
```

```
## vote
##     Conservative           Labour Liberal Democrat 
##          0.30295          0.47213          0.22492
```

``` r
    # try it out (produces error)

try(cv(m.beps, seed = 3465, criterion = BayesRuleMulti))
```

```
## Error in GetResponse.default(model) : non-vector response
```

``` r
    # define GetResponse() method for "multinom" objects

GetResponse.multinom <- function(model, ...) {
  insight::get_response(model)
}

head(GetResponse(m.beps))
```

```
## [1] Liberal Democrat Labour           Labour          
## [4] Labour           Labour           Labour          
## Levels: Conservative Labour Liberal Democrat
```

``` r
    # try it out (still fails)

try(cv(m.beps, seed = 3465, criterion = BayesRuleMulti))
```

```
## R RNG seed set to 3465
```

```
## Error in match.arg(type) : 'arg' should be one of "class", "probs"
```

``` r
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
```

```
## R RNG seed set to 3465
```

```
## 10-Fold Cross Validation
## criterion: BayesRuleMulti
## cross-validation criterion = 0.32459
## bias-adjusted cross-validation criterion = 0.32368
## 95% CI for bias-adjusted CV criterion = (0.30017, 0.34718)
## full-sample criterion = 0.31869
```

``` r
# Sec. 6: Comparing cv to other R software for cross-validation

  # example from rsample package vignette (modified to use LOO CV)

library("rsample")

data("attrition", package = "modeldata")
nrow(attrition)
```

```
## [1] 1470
```

``` r
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
```

```
## 
## Attaching package: 'purrr'
```

```
## The following object is masked from 'package:car':
## 
##     some
```

```
## The following objects are masked from 'package:foreach':
## 
##     accumulate, when
```

``` r
rs_obj$results <- map(rs_obj$splits, holdout_results, mod_form)
rs_obj$accuracy <- map_dbl(rs_obj$results, function(x) mean(x$correct))
1 - summary(rs_obj$accuracy) # error rate is complement of accuracy
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   0.000   0.000   0.161   0.000   0.000
```

``` r
    # script for direct computation of LOO CV for the example

mod.attrition <- glm(mod_form, data = attrition, family = binomial)
n <- nrow(attrition)
yhat <- numeric(n)
for (i in 1:n){
  m <- update(mod.attrition, data = attrition[-i, ])
  yhat[i] <- predict(m, newdata = attrition[i, ], type = "response")
}
mean(((attrition$Attrition == "Yes") - round(yhat)) ^ 2)
```

```
## [1] 0.16122
```

``` r
    # LOO CV for the example using the caret package

library("ggplot2")
```

```
## 
## Attaching package: 'ggplot2'
```

```
## The following objects are masked _by_ '.GlobalEnv':
## 
##     layer, mpg
```

```
## The following object is masked from 'package:latticeExtra':
## 
##     layer
```

``` r
library("caret")
```

```
## 
## Attaching package: 'caret'
```

```
## The following object is masked from 'package:purrr':
## 
##     lift
```

``` r
train(x = attrition[, c("JobSatisfaction", "Gender", "MonthlyIncome")],
      y = attrition$Attrition, method = "glm",
      trControl = trainControl(method = "LOOCV"))
```

```
## Generalized Linear Model 
## 
## 1470 samples
##    3 predictor
##    2 classes: 'No', 'Yes' 
## 
## No pre-processing
## Resampling: Leave-One-Out Cross-Validation 
## Summary of sample sizes: 1469, 1469, 1469, 1469, 1469, 1469, ... 
## Resampling results:
## 
##   Accuracy  Kappa
##   0.83878   0
```

``` r
    # LOO cv for the example using cv::cv()

(cv.attrition <- cv(mod.attrition, k = "loo", criterion = BayesRule))
```

```
## cross-validation criterion (BayesRule) = 0.16122
```

``` r
1 - mean(rs_obj$accuracy)
```

```
## [1] 0.16122
```

``` r
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
```

```
## Unit: relative
##          expr    min     lq   mean median     uq  max neval  cld
##  cv.hatvalues    1.0    1.0    1.0    1.0    1.0    1    10 a   
##       cv.wood   78.4   69.6   73.5   67.2   59.7  125    10 a   
##      cv.exact 3390.0 3100.0 2950.0 3040.0 2710.0 2640    10  b  
##        direct 3960.0 3530.0 3380.0 3460.0 3030.0 3080    10  bc 
##       rsample 4380.0 3890.0 3730.0 3800.0 3420.0 3470    10   c 
##         caret 5040.0 4520.0 4570.0 4370.0 3840.0 6090    10    d
```

``` r
# session info at end of script:

sessionInfo()
```

```
## R version 4.4.1 (2024-06-14)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sonoma 14.6.1
## 
## Matrix products: default
## BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
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
##  [1] caret_6.0-94        ggplot2_3.5.1       purrr_1.0.2        
##  [4] rsample_1.2.1       nnet_7.3-19         MASS_7.3-61        
##  [7] latticeExtra_0.6-30 lattice_0.22-6      glmmTMB_1.1.9      
## [10] lme4_1.1-35.5       Matrix_1.7-0        car_3.1-3          
## [13] carData_3.0-5       cv_2.0.3            doParallel_1.0.17  
## [16] iterators_1.0.14    foreach_1.5.2       knitr_1.48         
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.2.1      timeDate_4032.109    
##  [3] dplyr_1.1.4           TH.data_1.1-2        
##  [5] pROC_1.18.5           rpart_4.1.23         
##  [7] digest_0.6.36         timechange_0.3.0     
##  [9] lifecycle_1.0.4       survival_3.7-0       
## [11] magrittr_2.0.3        compiler_4.4.1       
## [13] rlang_1.1.4           tools_4.4.1          
## [15] utf8_1.2.4            data.table_1.15.4    
## [17] interp_1.1-6          plyr_1.8.9           
## [19] RColorBrewer_1.1-3    multcomp_1.4-26      
## [21] abind_1.4-5           withr_3.0.1          
## [23] numDeriv_2016.8-1.1   effects_4.2-3        
## [25] stats4_4.4.1          grid_4.4.1           
## [27] fansi_1.0.6           e1071_1.7-14         
## [29] colorspace_2.1-1      future_1.34.0        
## [31] globals_0.16.3        scales_1.3.0         
## [33] insight_0.20.2        cli_3.6.3            
## [35] mvtnorm_1.2-5         survey_4.4-2         
## [37] generics_0.1.3        future.apply_1.11.2  
## [39] rstudioapi_0.16.0     reshape2_1.4.4       
## [41] proxy_0.4-27          minqa_1.2.8          
## [43] DBI_1.2.3             stringr_1.5.1        
## [45] splines_4.4.1         mitools_2.4          
## [47] vctrs_0.6.5           hardhat_1.4.0        
## [49] boot_1.3-30           sandwich_3.1-0       
## [51] Formula_1.2-5         listenv_0.9.1        
## [53] jpeg_0.1-10           gower_1.0.1          
## [55] tidyr_1.3.1           recipes_1.1.0        
## [57] glue_1.7.0            parallelly_1.38.0    
## [59] nloptr_2.1.1          codetools_0.2-20     
## [61] stringi_1.8.4         lubridate_1.9.3      
## [63] gtable_0.3.5          deldir_2.0-4         
## [65] munsell_0.5.1         tibble_3.2.1         
## [67] furrr_0.3.1           pillar_1.9.0         
## [69] ipred_0.9-15          lava_1.8.0           
## [71] R6_2.5.1              TMB_1.9.14           
## [73] microbenchmark_1.4.10 evaluate_0.24.0      
## [75] png_0.1-8             backports_1.5.0      
## [77] broom_1.0.6           class_7.3-22         
## [79] Rcpp_1.0.13           prodlim_2024.06.25   
## [81] nlme_3.1-165          mgcv_1.9-1           
## [83] xfun_0.46             ModelMetrics_1.2.2.2 
## [85] zoo_1.8-12            pkgconfig_2.0.3
```

