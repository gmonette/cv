# Cross-Validate a Model-Selection Procedure

The [`cv()`](https://gmonette.github.io/cv/reference/cv.md) `"function"`
method is a general function to cross-validate a model-selection
procedure, such as the following: `selectStepAIC()` is a procedure that
applies the [`stepAIC()`](https://rdrr.io/pkg/MASS/man/stepAIC.html)
model-selection function in the MASS package; `selectTrans()` is a
procedure for selecting predictor and response transformations in
regression, which uses the
[`powerTransform()`](https://rdrr.io/pkg/car/man/powerTransform.html)
function in the car package; `selectTransAndStepAIC()` combines
predictor and response transformations with predictor selection; and
`selectModelList()` uses cross-validation to select a model from a list
of models created by
[`models()`](https://gmonette.github.io/cv/reference/cv.modList.md) and
employs (meta) cross-validation to assess the predictive accuracy of
this procedure.

## Usage

``` r
# S3 method for class '`function`'
cv(
  model,
  data,
  criterion = mse,
  k = 10L,
  reps = 1L,
  seed = NULL,
  working.model = NULL,
  y.expression = NULL,
  confint = n >= 400L,
  level = 0.95,
  details = k <= 10L,
  save.model = FALSE,
  ncores = 1L,
  ...
)

selectStepAIC(
  data,
  indices,
  model,
  criterion = mse,
  AIC = TRUE,
  details = TRUE,
  save.model = FALSE,
  ...
)

selectTrans(
  data,
  indices,
  details = TRUE,
  save.model = FALSE,
  model,
  criterion = mse,
  predictors,
  response,
  family = c("bcPower", "bcnPower", "yjPower", "basicPower"),
  family.y = c("bcPower", "bcnPower", "yjPower", "basicPower"),
  rounded = TRUE,
  ...
)

selectTransStepAIC(
  data,
  indices,
  details = TRUE,
  save.model = FALSE,
  model,
  criterion = mse,
  predictors,
  response,
  family = c("bcPower", "bcnPower", "yjPower", "basicPower"),
  family.y = c("bcPower", "bcnPower", "yjPower", "basicPower"),
  rounded = TRUE,
  AIC = TRUE,
  ...
)

selectModelList(
  data,
  indices,
  model,
  criterion = mse,
  k = 10L,
  k.meta = k,
  details = k <= 10L,
  save.model = FALSE,
  seed = FALSE,
  quietly = TRUE,
  ...
)

compareFolds(object, digits = 3, ...)

# S3 method for class 'cvSelect'
coef(object, average, NAs = 0, ...)

# S3 method for class 'cvSelect'
cvInfo(
  object,
  what = c("CV criterion", "adjusted CV criterion", "full CV criterion", "confint", "SE",
    "k", "seed", "method", "criterion name", "selected model"),
  ...
)
```

## Arguments

- model:

  a regression model object fit to data, or for the
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) `"function"`
  method, a model-selection procedure function (see Details).

- data:

  full data frame for model selection.

- criterion:

  a CV criterion ("cost" or lack-of-fit) function.

- k:

  perform k-fold cross-validation (default is 10); `k` may be a number
  or `"loo"` or `"n"` for n-fold (leave-one-out) cross-validation.

- reps:

  number of times to replicate k-fold CV (default is `1`)

- seed:

  for R's random number generator; not used for n-fold cross-validation.
  If not explicitly set, a seed is randomly generated and saved to make
  the results reproducible. In some cases, for internal use only, `seed`
  is set to `FALSE` to suppress automatically setting the seed.

- working.model:

  a regression model object fit to data, typically to begin a
  model-selection process; for use with `selectModelList()`, a list of
  competing models created by
  [`models()`](https://gmonette.github.io/cv/reference/cv.modList.md).

- y.expression:

  normally the response variable is found from the `model` or
  `working.model` argument; but if, for a particular selection
  procedure, the `model` or `working.model` argument is absent, or if
  the response can't be inferred from the model, the response can be
  specified by an expression, such as `expression(log(income))`, to be
  evaluated within the data set provided by the `data` argument.

- confint:

  if `TRUE` (the default if the number of cases is 400 or greater),
  compute a confidence interval for the bias-corrected CV criterion, if
  the criterion is the average of casewise components.

- level:

  confidence level (default `0.95`).

- details:

  if `TRUE`, save detailed information about the value of the CV
  criterion for the cases in each fold and the regression coefficients
  (and possibly other information) with that fold deleted; default is
  `TRUE` if `k` is 10 or smaller, `FALSE` otherwise.

- save.model:

  save the model that's selected using the *full* data set (default,
  `FALSE`).

- ncores:

  number of cores to use for parallel computations (default is `1`,
  i.e., computations aren't done in parallel)

- ...:

  for
  [`cvSelect()`](https://gmonette.github.io/cv/reference/cvCompute.md)
  and the [`cv()`](https://gmonette.github.io/cv/reference/cv.md)
  `"function"` method, arguments to be passed to `procedure()`; for
  `selectStepAIC()` and `selectTransStepAIC()`, arguments to be passed
  to `stepAIC()`.

- indices:

  indices of cases in data defining the current fold.

- AIC:

  if `TRUE` (the default) use the AIC as the model-selection criterion;
  if `FALSE`, use the BIC. The `k` argument to
  [`stepAIC()`](https://rdrr.io/pkg/MASS/man/stepAIC.html) is set
  accordingly (note that this is distinct from the number of folds `k`).

- predictors:

  character vector of names of the predictors in the model to transform;
  if missing, no predictors will be transformed.

- response:

  name of the response variable; if missing, the response won't be
  transformed.

- family:

  transformation family for the predictors, one of
  `"bcPower", "bcnPower", "yjPower", "basicPower"`, with `"bcPower"` as
  the default. These are the names of transformation functions in the
  car package; see
  [`bcPower()`](https://rdrr.io/pkg/car/man/bcPower.html).

- family.y:

  transformation family for the response, with `"bcPower"` as the
  default.

- rounded:

  if `TRUE` (the default) use nicely rounded versions of the estimated
  transformation parameters (see
  [`bcPower()`](https://rdrr.io/pkg/car/man/bcPower.html)).

- k.meta:

  the number of folds for meta CV; defaults to the value of `k`; may be
  specified as `"loo"` or `"n"` as well as an integer.

- quietly:

  if `TRUE` (the default), simple messages (for example about the value
  to which the random-number generator seed is set), but not warnings or
  errors, are suppressed.

- object:

  an object of class `"cvSelect"`.

- digits:

  significant digits for printing coefficients (default `3`).

- average:

  if supplied, a function, such as `mean` or `median`, to use us in
  averaging estimates across folds; if missing, the estimates for each
  fold are returned.

- NAs:

  values to substitute for `NA`s in calculating averaged estimates; the
  default, `0`, is appropriate, e.g., for regression coefficients; the
  value `1` might be appropriate for power-transformation estimates.

- what:

  the information to extract from a `"cvSelect"` object, one of:
  `"CV criterion"`, `"adjusted CV criterion"`, `"full CV criterion"`
  (the CV criterion applied to the model fit to the full data set),
  `"SE"` (the standard error of the adjusted CV criterion), `"confint"`
  (confidence interval for the adjusted CV criterion), `"k"`, (the
  number of folds), `"seed"` (the seed employed for R's random-number
  generator), `"method"` (the computational method employed, e.g., for a
  `"lm"` model object), `"criterion name"` (the CV criterion employed),
  or `"selected model"` (the model object for the model that was
  selected); not all of these elements may be present, in which case
  [`cvInfo()`](https://gmonette.github.io/cv/reference/cv.md) would
  return `NULL`.

## Value

An object of class `"cvSelect"`, inheriting from class `"cv"`, with the
CV criterion (`"CV crit"`), the bias-adjusted CV criterion
(`"adj CV crit"`), the criterion for the model applied to the full data
(`"full crit"`), the confidence interval and level for the bias-adjusted
CV criterion (`"confint"`), the number of folds (`"k"`), the seed for
R's random-number generator (`"seed"`), and (optionally) a list of
coefficients (or, in the case of `selectTrans()`, estimated
transformation parameters, and in the case of `selectTransAndStepAIC()`,
both regression coefficients and transformation parameters) for the
selected models for each fold (`"coefficients"`). If `reps` \> `1`, then
an object of class `c("cvSelectList", "cvList")` is returned, which is
literally a list of `c("cvSelect", "cv")` objects.

## Details

The model-selection function supplied as the `procedure` (for
[`cvSelect()`](https://gmonette.github.io/cv/reference/cvCompute.md)) or
`model` (for [`cv()`](https://gmonette.github.io/cv/reference/cv.md))
argument should accept the following arguments:

- `data`:

  set to the `data` argument to
  [`cvSelect()`](https://gmonette.github.io/cv/reference/cvCompute.md)
  or [`cv()`](https://gmonette.github.io/cv/reference/cv.md).

- `indices`:

  the indices of the rows of `data` defining the current fold; if
  missing, the model-selection procedure is applied to the full `data`.

- other arguments:

  to be passed via `...` from
  [`cvSelect()`](https://gmonette.github.io/cv/reference/cvCompute.md)
  or [`cv()`](https://gmonette.github.io/cv/reference/cv.md).

`procedure()` or `model()` should return a list with the following named
elements: `fit.i`, the vector of predicted values for the cases in the
current fold computed from the model omitting these cases; `crit.all.i`,
the CV criterion computed for all of the cases using the model omitting
the current fold; and (optionally) `coefficients`, parameter estimates
from the model computed omitting the current fold.

When the `indices` argument is missing, `procedure()` returns the
cross-validation criterion for all of the cases based on the model fit
to all of the cases.

For examples of model-selection functions for the `procedure` argument,
see the code for `selectStepAIC()`, `selectTrans()`, and
`selectTransAndStepAIC()`.

For additional information, see the "Cross-validating model selection"
vignette (`vignette("cv-select", package="cv")`) and the "Extending the
cv package" vignette
([`vignette("cv-extend", package="cv")`](https://gmonette.github.io/cv/articles/cv-extend.md)).

## Functions

- `` cv(`function`) ``:
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) method for
  applying a model model-selection (or specification) procedure.

- `selectStepAIC()`: select a regression model using the
  [`stepAIC()`](https://rdrr.io/pkg/MASS/man/stepAIC.html) function in
  the MASS package.

- `selectTrans()`: select transformations of the predictors and response
  using
  [`powerTransform()`](https://rdrr.io/pkg/car/man/powerTransform.html)
  in the car package.

- `selectTransStepAIC()`: select transformations of the predictors and
  response, and then select predictors.

- `selectModelList()`: select a model using (meta) CV.

- `compareFolds()`: print the coefficients from the selected models for
  the several folds.

- `coef(cvSelect)`: extract the coefficients from the selected models
  for the several folds and possibly average them.

## See also

[`stepAIC`](https://rdrr.io/pkg/MASS/man/stepAIC.html),
[`bcPower`](https://rdrr.io/pkg/car/man/bcPower.html),
[`powerTransform`](https://rdrr.io/pkg/car/man/powerTransform.html),
[`cv`](https://gmonette.github.io/cv/reference/cv.md).

## Examples

``` r
  if (requireNamespace("ISLR2", quietly=TRUE)){
    withAutoprint({
      data("Auto", package="ISLR2")
      m.auto <- lm(mpg ~ . - name - origin - year, data=Auto)
      cv(selectStepAIC, Auto, seed=123, working.model=m.auto,
         AIC=FALSE, k=5, reps=2) # via BIC
      # The user is invited to increase reps and to try using AIC with the
      # following line:
      # cv(selectStepAIC, Auto, seed=123, working.model=m.auto) # via AIC
    })
  } else {
    cat("\n install the 'ISLR2' package to run these examples\n")
  }
#> > data("Auto", package = "ISLR2")
#> > m.auto <- lm(mpg ~ . - name - origin - year, data = Auto)
#> > cv(selectStepAIC, Auto, seed = 123, working.model = m.auto, AIC = FALSE, 
#> +     k = 5, reps = 2)
#> R RNG seed set to 123
#> R RNG seed set to 319729
#> cross-validation criterion (mse)
#> Replicate 1: 18.04787 
#> Replicate 2: 18.40478 
#> Average: 18.16684 
if (requireNamespace("carData", quietly=TRUE)){
  withAutoprint({
    data("Prestige", package="carData")
    m.pres <- lm(prestige ~ income + education + women,
               data=Prestige)
    cvt <- cv(selectTrans, data=Prestige, working.model=m.pres, seed=123,
            predictors=c("income", "education", "women"),
            response="prestige", family="yjPower")
    cvt
    compareFolds(cvt)
    coef(cvt, average=median, NAs=1) # NAs not really needed here
    cv(m.pres, seed=123)
  })
  } else {
  cat("install the 'carData' package to run these examples\n")
}
#> > data("Prestige", package = "carData")
#> > m.pres <- lm(prestige ~ income + education + women, data = Prestige)
#> > cvt <- cv(selectTrans, data = Prestige, working.model = m.pres, seed = 123, 
#> +     predictors = c("income", "education", "women"), response = "prestige", family = "yjPower")
#> R RNG seed set to 123
#> > cvt
#> cross-validation criterion (mse) = 58.68193 
#> > compareFolds(cvt)
#> CV criterion by folds:
#>    fold.1    fold.2    fold.3    fold.4    fold.5    fold.6    fold.7    fold.8 
#> 100.01805  27.79099  66.84054  60.97286  30.75663  17.98406 110.27446  66.35537 
#>    fold.9   fold.10 
#>  37.32475  67.45709 
#> 
#> Coefficients by folds:
#>         lam.education lam.income lam.women lambda
#> Fold 1          1.000      0.330     0.330      1
#> Fold 2          1.000      0.330     0.330      1
#> Fold 3          1.000      0.330     0.000      1
#> Fold 4          1.000      0.330     0.330      1
#> Fold 5          1.000      0.330     0.000      1
#> Fold 6          1.000      0.330     0.000      1
#> Fold 7          1.000      0.000     0.000      1
#> Fold 8          1.000      0.330     0.000      1
#> Fold 9          1.000      0.500     0.330      1
#> Fold 10         1.000      0.330     0.159      1
#> > coef(cvt, average = median, NAs = 1)
#> lam.education    lam.income     lam.women        lambda 
#>    1.00000000    0.33000000    0.07948918    1.00000000 
#> > cv(m.pres, seed = 123)
#> R RNG seed set to 123
#> cross-validation criterion (mse) = 66.91929 
if (requireNamespace("ISLR2", quietly=TRUE)){
withAutoprint({
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
          response="mpg", AIC=FALSE)
cvs
compareFolds(cvs)
})
}
#> > Auto$year <- as.factor(Auto$year)
#> > Auto$origin <- factor(Auto$origin, labels = c("America", "Europe", "Japan"))
#> > rownames(Auto) <- make.names(Auto$name, unique = TRUE)
#> > Auto$name <- NULL
#> > m.auto <- lm(mpg ~ ., data = Auto)
#> > cvs <- cv(selectTransStepAIC, data = Auto, seed = 76692, working.model = m.auto, 
#> +     criterion = medAbsErr, predictors = c("cylinders", "displacement", "horsepower", 
#> +         "weight", "acceleration"), response = "mpg", AIC = FALSE)
#> R RNG seed set to 76692
#> > cvs
#> cross-validation criterion (medAbsErr) = 1.476272 
#> > compareFolds(cvs)
#> CV criterion by folds:
#>   fold.1   fold.2   fold.3   fold.4   fold.5   fold.6   fold.7   fold.8 
#> 1.563930 1.562867 1.369843 1.282772 1.246066 1.582635 1.340277 1.183112 
#>   fold.9  fold.10 
#> 1.138892 1.585400 
#> 
#> Coefficients by folds:
#>         (Intercept) horsepower lam.acceleration lam.cylinders lam.displacement
#> Fold 1      9.71384   -0.17408          0.50000      -0.50000          0.18555
#> Fold 2      9.09355   -0.31285          0.00000      -0.50000          0.19452
#> Fold 3      9.61824   -0.19248          0.00000      -0.50000          0.15461
#> Fold 4      9.49410   -0.25380          0.00000      -0.50000          0.12800
#> Fold 5      9.14403   -0.14934          0.00000      -0.50000          0.18388
#> Fold 6      9.63481   -0.16739          0.00000      -0.50000          0.17122
#> Fold 7      9.60487   -0.18119          0.00000       0.00000          0.19063
#> Fold 8      8.92286   -0.29237          0.00000      -0.50000          0.33000
#> Fold 9      8.71492   -0.22484          0.00000      -0.50000          0.16881
#> Fold 10     9.61727   -0.17086          0.00000       0.00000          0.16638
#>         lam.horsepower lam.weight   lambda   weight   year71   year72   year73
#> Fold 1         0.00000    0.00000  0.00000 -0.74636  0.03764 -0.00327 -0.02477
#> Fold 2         0.00000    0.00000  0.00000 -0.48573  0.02141 -0.01416 -0.03720
#> Fold 3         0.00000    0.00000  0.00000 -0.72085  0.01128 -0.02569 -0.03872
#> Fold 4         0.00000    0.00000  0.00000 -0.57844  0.02269 -0.02130 -0.05069
#> Fold 5         0.00000    0.00000  0.00000 -0.69081  0.02531 -0.01062 -0.04625
#> Fold 6         0.00000    0.00000  0.00000 -0.74049  0.02456  0.00759 -0.03412
#> Fold 7         0.00000    0.00000  0.00000 -0.72741  0.02341 -0.01438 -0.04241
#> Fold 8         0.00000    0.00000  0.00000 -0.48913  0.02720 -0.01844 -0.05477
#> Fold 9         0.00000    0.00000  0.00000 -0.48188  0.00896 -0.03530 -0.04815
#> Fold 10        0.00000    0.00000  0.00000 -0.73550  0.02937 -0.00899 -0.03814
#>           year74   year75   year76   year77   year78   year79   year80   year81
#> Fold 1   0.05606  0.07080  0.07250  0.14420  0.14281  0.23266  0.35127  0.25635
#> Fold 2   0.04374  0.04010  0.06713  0.13116  0.14857  0.21825  0.33253  0.26173
#> Fold 3   0.05187  0.03837  0.06399  0.11593  0.12601  0.20499  0.32821  0.24478
#> Fold 4   0.04488  0.04238  0.05907  0.12717  0.13902  0.22986  0.33253  0.25407
#> Fold 5   0.05039  0.05596  0.07044  0.13356  0.14724  0.24675  0.33331  0.26938
#> Fold 6   0.06266  0.06940  0.07769  0.14211  0.14647  0.23532  0.34761  0.26737
#> Fold 7   0.04368  0.03327  0.07175  0.12777  0.13816  0.23722  0.33822  0.27453
#> Fold 8   0.04670  0.06595  0.08192  0.13289  0.13899  0.23003  0.33022  0.25909
#> Fold 9   0.01986  0.02981  0.05843  0.10584  0.11625  0.20625  0.31601  0.23350
#> Fold 10  0.05408  0.04881  0.07862  0.14101  0.14313  0.23258  0.35649  0.26214
#>           year82 acceleration displacement cylinders originEurope originJapan
#> Fold 1   0.30546                                                             
#> Fold 2   0.30876     -0.19017     -0.03233                                   
#> Fold 3   0.29204                                                             
#> Fold 4   0.28048     -0.14790               -0.30064                         
#> Fold 5   0.32594                                          0.06261        0.04
#> Fold 6   0.33062                                                             
#> Fold 7   0.30436                                                             
#> Fold 8   0.30500     -0.17392     -0.01683                                   
#> Fold 9   0.29284     -0.14623     -0.05409                                   
#> Fold 10  0.32421                                                             
# \donttest{
data("Duncan", package="carData")
m1 <- lm(prestige ~ income + education, data=Duncan)
m2 <- lm(prestige ~ income + education + type, data=Duncan)
m3 <- lm(prestige ~ (income + education)*type, data=Duncan)
summary(cv.sel <- cv(selectModelList, data=Duncan, seed=5963,
                     working.model=models(m1, m2, m3),
                     save.model=TRUE)) # meta CV
#> R RNG seed set to 5963
#> 10-Fold Cross Validation
#> criterion: mse
#> cross-validation criterion = 112.0172
#> bias-adjusted cross-validation criterion = 139.8472
#> full-sample criterion = 113.7905
cvInfo(cv.sel, "selected model")
#> 
#> Call:
#> lm(formula = prestige ~ income + education + type, data = Duncan)
#> 
#> Coefficients:
#> (Intercept)       income    education     typeprof       typewc  
#>     -0.1850       0.5975       0.3453      16.6575     -14.6611  
#> 
# }
```
