# Cross-Validate Several Models Fit to the Same Data

A [`cv()`](https://gmonette.github.io/cv/reference/cv.md) method for an
object of class `"modlist"`, created by the `models()` function. This
[`cv()`](https://gmonette.github.io/cv/reference/cv.md) method
simplifies the process of cross-validating several models on the same
set of CV folds and may also be used for meta CV, where CV is used to
select one from among several models. `models()` performs some "sanity"
checks, warning if the models are of different classes, and reporting an
error if they are fit to apparently different data sets or different
response variables.

## Usage

``` r
# S3 method for class 'modList'
cv(
  model,
  data,
  criterion = mse,
  k,
  reps = 1L,
  seed,
  quietly = TRUE,
  meta = FALSE,
  ...
)

models(...)

# S3 method for class 'cvModList'
print(x, ...)

# S3 method for class 'cvModList'
summary(object, ...)

# S3 method for class 'cvModList'
plot(
  x,
  y,
  spread = c("range", "sd"),
  confint = TRUE,
  xlab = "",
  ylab,
  main,
  axis.args = list(labels = names(x), las = 3L),
  col = palette()[2L],
  lwd = 2L,
  grid = TRUE,
  ...
)

# S3 method for class 'cvModList'
as.data.frame(x, row.names = NULL, optional = TRUE, ...)
```

## Arguments

- model:

  a list of regression model objects, created by `models()`.

- data:

  (required) the data set to which the models were fit.

- criterion:

  the CV criterion ("cost" or lack-of-fit) function, defaults to
  [`mse`](https://gmonette.github.io/cv/reference/cost-functions.md).

- k:

  the number of CV folds; may be omitted, in which case the value will
  depend on the default for the
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) method invoked
  for the individual models.

- reps:

  number of replications of CV for each model (default is 1).

- seed:

  (optional) seed for R's pseudo-random-number generator, to be used to
  create the same set of CV folds for all of the models; if omitted, a
  seed will be randomly generated and saved. Not used for leave-one-out
  CV.

- quietly:

  if `TRUE` (the default), simple messages (for example about the value
  to which the random-number generator seed is set), but not warnings or
  errors, are suppressed.

- meta:

  if `TRUE` (the default is `FALSE`), cross-validation is performed
  recursively to select a "best" model deleting each fold in turn by
  calculating the CV estimate of the criterion for the remaining folds;
  this is equivalent to employing the
  [`selectModelList()`](https://gmonette.github.io/cv/reference/cv.function.md)
  model-selection procedure.

- ...:

  for `cv.modList()`, additional arguments to be passed to the
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) method applied
  to each model.

  For `models()`, two or more competing models fit to the the same data;
  the several models may be named. It is also possible to specify a
  single argument, which should then be list of models (which has the
  effect of turning a list of models into a `"modList"` object).

  For the [`print()`](https://rdrr.io/r/base/print.html) method,
  arguments to be passed to the
  [`print()`](https://rdrr.io/r/base/print.html) method for the
  individual model cross-validations.

  For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method, arguments to be passed to the base
  [`plot()`](https://rdrr.io/r/base/plot.html) function.

- x:

  an object of class `"cvModList"` to be printed or plotted.

- object:

  an object to summarize.

- y:

  the name of the element in each `"cv"` object to be plotted; defaults
  to `"adj CV crit"`, if it exists, or to `"CV crit"`.

- spread:

  if `"range"`, the default, show the range of CV criteria for each
  model along with their average; if `"sd"`, show the average plus or
  minus 1 standard deviation.

- confint:

  if `TRUE` (the default) and if confidence intervals are in any of the
  `"cv"` objects, then plot the confidence intervals around the CV
  criteria.

- xlab:

  label for the x-axis (defaults to blank).

- ylab:

  label for the y-axis (if missing, a label is constructed).

- main:

  main title for the graph (if missing, a label is constructed).

- axis.args:

  a list of arguments for the
  [`axis()`](https://rdrr.io/r/graphics/axis.html) function, used to
  draw the horizontal axis. In addition to the axis arguments given
  explicitly, `side=1` (the horizontal axis) and `at=seq(along=x)`
  (i.e., 1 to the number of models) are used and can't be modified.

- col:

  color for the line and points, defaults to the second element of the
  color palette; see
  [`palette()`](https://rdrr.io/r/grDevices/palette.html).

- lwd:

  line width for the line (defaults to 2).

- grid:

  if `TRUE` (the default), include grid lines on the graph.

- row.names:

  optional row names for the result, defaults to `NULL`.

- optional:

  to match the
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) generic
  function; if `FALSE` (the default is `TRUE`), then the names of the
  columns of the returned data frame, including the names of
  coefficients, are coerced to syntactically correct names.

## Value

`models()` returns a `"modList"` object, the
[`cv()`](https://gmonette.github.io/cv/reference/cv.md) method for which
returns a `"cvModList"` object, or, when `meta=TRUE`, an object of class
`c("cvSelect", "cv")`.

## Functions

- `cv(modList)`: [`cv()`](https://gmonette.github.io/cv/reference/cv.md)
  method for `"modList"` objects.

- `models()`: create a list of models.

- `print(cvModList)`: [`print()`](https://rdrr.io/r/base/print.html)
  method for `"cvModList"` objects.

- `summary(cvModList)`:
  [`summary()`](https://rdrr.io/r/base/summary.html) method for
  `"cvModList"` objects.

- `plot(cvModList)`:
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
  `"cvModList"` objects.

- `as.data.frame(cvModList)`:
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) method
  for `"cvModList"` objects.

## See also

[`cv`](https://gmonette.github.io/cv/reference/cv.md),
[`cv.merMod`](https://gmonette.github.io/cv/reference/cv.merMod.md),
[`selectModelList`](https://gmonette.github.io/cv/reference/cv.function.md).

## Examples

``` r
if(interactive()){
  if (requireNamespace("carData", quietly=TRUE)){
  withAutoprint({
  data("Duncan", package="carData")
  m1 <- lm(prestige ~ income + education, data=Duncan)
  m2 <- lm(prestige ~ income + education + type, data=Duncan)
  m3 <- lm(prestige ~ (income + education)*type, data=Duncan)
  (cv.models <- cv(models(m1=m1, m2=m2, m3=m3),
                   data=Duncan, seed=7949, reps=5))
  D.cv.models <- as.data.frame(cv.models)
  head(D.cv.models)
  summary(D.cv.models, criterion ~ model + rep, include="folds")
  plot(cv.models)
  (cv.models.ci <- cv(models(m1=m1, m2=m2, m3=m3),
                      data=Duncan, seed=5963, confint=TRUE, level=0.50))
                   # nb: n too small for accurate CIs
  plot(cv.models.ci)
  (cv.models.meta <- cv(models(m1=m1, m2=m2, m3=m3),
                        data=Duncan, seed=5963,
                        meta=TRUE, save.model=TRUE))
  cvInfo(cv.models.meta, "selected model")
  })
  } else {
  cat("install the 'carData' package to run these examples\n")
  }
}
```
