# Utility Functions for the cv Package

These functions are primarily useful for writing methods for the
[`cv()`](https://gmonette.github.io/cv/reference/cv.md) generic
function. They are used internally in the package and can also be used
for extensions (see the vignette "Extending the cv package,
[`vignette("cv-extend", package="cv")`](https://gmonette.github.io/cv/articles/cv-extend.md)).

## Usage

``` r
cvCompute(
  model,
  data = insight::get_data(model),
  criterion = mse,
  criterion.name,
  k = 10L,
  reps = 1L,
  seed,
  details = k <= 10L,
  confint,
  level = 0.95,
  method = NULL,
  ncores = 1L,
  type = "response",
  start = FALSE,
  f,
  fPara = f,
  locals = list(),
  model.function = NULL,
  model.function.name = NULL,
  ...
)

cvMixed(
  model,
  package,
  data = insight::get_data(model),
  criterion = mse,
  criterion.name,
  k,
  reps = 1L,
  confint,
  level = 0.95,
  seed,
  details,
  ncores = 1L,
  clusterVariables,
  predict.clusters.args = list(object = model, newdata = data),
  predict.cases.args = list(object = model, newdata = data),
  fixed.effects,
  ...
)

cvSelect(
  procedure,
  data,
  criterion = mse,
  criterion.name,
  model,
  y.expression,
  k = 10L,
  confint = n >= 400,
  level = 0.95,
  reps = 1L,
  save.coef,
  details = k <= 10L,
  save.model = FALSE,
  seed,
  ncores = 1L,
  ...
)

folds(n, k)

fold(folds, i_, ...)

# S3 method for class 'folds'
fold(folds, i_, ...)

# S3 method for class 'folds'
print(x, ...)

GetResponse(model, ...)

# Default S3 method
GetResponse(model, ...)

# S3 method for class 'merMod'
GetResponse(model, ...)

# S3 method for class 'lme'
GetResponse(model, ...)

# S3 method for class 'glmmTMB'
GetResponse(model, ...)

# S3 method for class 'modList'
GetResponse(model, ...)

checkFormula(model, data.names)
```

## Arguments

- model:

  a regression model object.

- data:

  data frame to which the model was fit (not usually necessary, except
  for `cvSelect()`).

- criterion:

  cross-validation criterion ("cost" or lack-of-fit) function of form
  `f(y, yhat)` where `y` is the observed values of the response and
  `yhat` the predicted values; the default is
  [`mse`](https://gmonette.github.io/cv/reference/cost-functions.md)
  (the mean-squared error).

- criterion.name:

  a character string giving the name of the CV criterion function in the
  returned `"cv"` object).

- k:

  perform k-fold cross-validation (default is `10`); `k` may be a number
  or `"loo"` or `"n"` for n-fold (leave-one-out) cross-validation; for
  `folds()`, `k` must be a number.

- reps:

  number of times to replicate k-fold CV (default is `1`).

- seed:

  for R's random number generator; optional, if not supplied a random
  seed will be selected and saved; not needed for n-fold
  cross-validation.

- details:

  if `TRUE` (the default if the number of folds `k <= 10`), save
  detailed information about the value of the CV criterion for the cases
  in each fold and the regression coefficients with that fold deleted.

- confint:

  if `TRUE` (the default if the number of cases is 400 or greater),
  compute a confidence interval for the bias-corrected CV criterion, if
  the criterion is the average of casewise components.

- level:

  confidence level (default `0.95`).

- method:

  computational method to apply; use by some
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods.

- ncores:

  number of cores to use for parallel computations (default is `1`,
  i.e., computations aren't done in parallel).

- type:

  used by some [`cv()`](https://gmonette.github.io/cv/reference/cv.md)
  methods, such as the default method, where `type` is passed to the
  `type` argument of
  [`predict()`](https://rdrr.io/r/stats/predict.html); the default is
  `type="response"`, which is appropriate, e.g., for a `"glm"` model and
  may be recognized or ignored by
  [`predict()`](https://rdrr.io/r/stats/predict.html) methods for other
  model classes.

- start:

  used by some [`cv()`](https://gmonette.github.io/cv/reference/cv.md)
  methods; if `TRUE` (the default is `FALSE`), the `start` argument, set
  to the vector of regression coefficients for the model fit to the full
  data, is passed to [`update()`](https://rdrr.io/r/stats/update.html),
  possibly making the CV updates faster, e.g. for a GLM.

- f:

  function to be called by `cvCompute()` for each fold.

- fPara:

  function to be called by `cvCompute()` for each fold using parallel
  computation.

- locals:

  a named list of objects that are required in the local environment of
  `cvCompute()` for `f()` or `fPara()`.

- model.function:

  a regression function, typically for a new
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) method,
  residing in a package that's not a declared dependency of the cv
  package, e.g.,
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html).

- model.function.name:

  the quoted name of the regression function, e.g., `"multinom"`.

- ...:

  to match generic; passed to
  [`predict()`](https://rdrr.io/r/stats/predict.html) for the default
  method, and to `fPara()` (for parallel computations) in `cvCompute()`.

- package:

  the name of the package in which mixed-modeling function (or
  functions) employed resides; used to get the namespace of the package.

- clusterVariables:

  a character vector of names of the variables defining clusters for a
  mixed model with nested or crossed random effects; if missing,
  cross-validation is performed for individual cases rather than for
  clusters

- predict.clusters.args:

  a list of arguments to be used to predict the whole data set from a
  mixed model when performing CV on clusters; the first two elements
  should be `model` and `newdata`; see the "Extending the cv package"
  vignette
  ([`vignette("cv-extend", package="cv")`](https://gmonette.github.io/cv/articles/cv-extend.md)).

- predict.cases.args:

  a list of arguments to be used to predict the whole data set from a
  mixed model when performing CV on cases; the first two elements should
  be `model` and `newdata`; see the "Extending the cv package" vignette
  ([`vignette("cv-extend", package="cv")`](https://gmonette.github.io/cv/articles/cv-extend.md)).

- fixed.effects:

  a function to be used to compute fixed-effect coefficients for
  cluster-based CV when `details = TRUE`.

- procedure:

  a model-selection procedure function (see Details).

- y.expression:

  normally the response variable is found from the `model` argument; but
  if, for a particular selection procedure, the `model` argument is
  absent, or if the response can't be inferred from the model, the
  response can be specified by an expression, such as
  `expression(log(income))`, to be evaluated within the data set
  provided by the `data` argument.

- save.coef:

  save the coefficients from the selected models? Deprecated in favor of
  the `details` argument; if specified, `details` is set is set to the
  value of `save.coef`.

- save.model:

  save the model that's selected using the *full* data set.

- n:

  number of cases, for constructed folds.

- folds:

  an object of class `"folds"`.

- i\_:

  a fold number for an object of class `"folds"`.

- x:

  a `"cv"`, `"cvList"`, or `"folds"` object to be printed

- data.names:

  names of variables in the data set to which the model was fit; if
  missing, an attempt will be made to extract the data from the model.

## Value

The utility functions return various kinds of objects:

- `cvCompute()` returns an object of class `"cv"`, with the CV criterion
  (`"CV crit"`), the bias-adjusted CV criterion (`"adj CV crit"`), the
  criterion for the model applied to the full data (`"full crit"`), the
  confidence interval and level for the bias-adjusted CV criterion
  (`"confint"`), the number of folds (`"k"`), and the seed for R's
  random-number generator (`"seed"`). If `details=TRUE`, then the
  returned object will also include a `"details"` component, which is a
  list of two elements: `"criterion"`, containing the CV criterion
  computed for the cases in each fold; and `"coefficients"`, regression
  coefficients computed for the model with each fold deleted. Some
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods
  calling `cvCompute()` may return a subset of these components and may
  add additional information. If `reps` \> `1`, then an object of class
  `"cvList"` is returned, which is literally a list of `"cv"` objects.

- `cvMixed()` also returns an object of class `"cv"` or `"cvList"`.

- `cvSelect` returns an object of class `"cvSelect"` inheriting from
  `"cv"`, or an object of class `"cvSelectList"` inheriting from
  `"cvList"`.

- `folds()` returns an object of class folds, for which there are
  `fold()` and [`print()`](https://rdrr.io/r/base/print.html) methods.

- `GetResponse()` returns the (numeric) response variable from the
  model.

  The supplied `default` method returns the `model$y` component of the
  model object, or, if `model` is an S4 object, the result returned by
  the
  [`get_response()`](https://easystats.github.io/insight/reference/get_response.html)
  function in the insight package. If this result is `NULL`, the result
  of `model.response(model.frame(model))` is returned, checking in any
  case whether the result is a numeric vector.

  There are also `"lme"`, `"merMod"` and `"glmmTMB"` methods that
  convert factor responses to numeric 0/1 responses, as would be
  appropriate for a generalized linear mixed model with a binary
  response.

- `checkFormula()` returns `TRUE` if all variables in the model formula
  are also in the data to which the model is fit; `FALSE` is this is not
  the case (and q warning is printed); or `NA` if the function couldn't
  extract a model formula.

## Functions

- `cvCompute()`: used internally by
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods (not
  for direct use); exported to support new
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods.

- `cvMixed()`: used internally by
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods for
  mixed-effect models (not for direct use); exported to support new
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods.

- `cvSelect()`: used internally by
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods for
  cross-validating a model-selection procedure; may also be called
  directly for this purpose, but use via
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) is preferred.
  `cvSelect()` is exported primarily to support new model-selection
  procedures.

- `folds()`: used internally by
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods (not
  for direct use).

- `fold()`: to extract a fold from a `"folds"` object.

- `fold(folds)`: `fold()` method for `"folds"` objects.

- `print(folds)`: [`print()`](https://rdrr.io/r/base/print.html) method
  for `"folds"` objects.

- `GetResponse()`: function to return the response variable from a
  regression model.

- `GetResponse(default)`: default method.

- `GetResponse(merMod)`: `"merMod"` method.

- `GetResponse(lme)`: `"lme"` method.

- `GetResponse(glmmTMB)`: `"glmmTMB"` method.

- `GetResponse(modList)`: `"modList"` method.

- `checkFormula()`: check a model formula to determine whether it
  include variables not in the data to which the model was fit; prints a
  warning if this is not the case.

## See also

[`cv`](https://gmonette.github.io/cv/reference/cv.md),
[`cv.merMod`](https://gmonette.github.io/cv/reference/cv.merMod.md),
[`cv.function`](https://gmonette.github.io/cv/reference/cv.function.md).

## Examples

``` r
fit <- lm(mpg ~ gear, mtcars)
GetResponse(fit)
#>           Mazda RX4       Mazda RX4 Wag          Datsun 710      Hornet 4 Drive 
#>                21.0                21.0                22.8                21.4 
#>   Hornet Sportabout             Valiant          Duster 360           Merc 240D 
#>                18.7                18.1                14.3                24.4 
#>            Merc 230            Merc 280           Merc 280C          Merc 450SE 
#>                22.8                19.2                17.8                16.4 
#>          Merc 450SL         Merc 450SLC  Cadillac Fleetwood Lincoln Continental 
#>                17.3                15.2                10.4                10.4 
#>   Chrysler Imperial            Fiat 128         Honda Civic      Toyota Corolla 
#>                14.7                32.4                30.4                33.9 
#>       Toyota Corona    Dodge Challenger         AMC Javelin          Camaro Z28 
#>                21.5                15.5                15.2                13.3 
#>    Pontiac Firebird           Fiat X1-9       Porsche 914-2        Lotus Europa 
#>                19.2                27.3                26.0                30.4 
#>      Ford Pantera L        Ferrari Dino       Maserati Bora          Volvo 142E 
#>                15.8                19.7                15.0                21.4 

set.seed(123)
(ffs <- folds(n=22, k=5))
#> 5 folds of approximately 4 cases each
#>  fold 1: 15 19 14 3 10
#>  fold 2: 11 5 4 20 6
#>  fold 3: 9 18 16 21
#>  fold 4: 12 1 22 7
#>  fold 5: 17 13 8 2
fold(ffs, 2)
#> [1] 11  5  4 20  6
```
