# Cost Functions for Fitted Regression Models

Compute cost functions (cross-validation criteria) for fitted regression
models.

## Usage

``` r
mse(y, yhat)

rmse(y, yhat)

medAbsErr(y, yhat)

BayesRule(y, yhat)

BayesRule2(y, yhat)
```

## Arguments

- y:

  response

- yhat:

  fitted value

## Value

In general, cost functions should return a single numeric value
measuring lack-of-fit. `mse()` returns the mean-squared error; `rmse()`
returns the root-mean-squared error; `medAbsErr()` returns the median
absolute error; and `BayesRule()` and `BayesRule2()` return the
proportion of misclassified cases.

## Details

Cost functions (cross-validation criteria) are meant to measure
lack-of-fit. Several cost functions are provided:

1.  `mse()` returns the mean-squared error of prediction for a numeric
    response variable `y` and predictions `yhat`; and `rmse()` returns
    the root-mean-squared error and is just the square-root of `mse()`.

2.  `medAbsErr()` returns the median absolute error of prediction for a
    numeric response `y` and predictions `yhat`.

3.  `BayesRule()` and `BayesRule2()` report the proportion of incorrect
    predictions for a dichotomous response variable `y`, assumed coded
    (or coercible to) `0` and `1`. The `yhat` values are predicted
    probabilities and are rounded to 0 or 1. The distinction between
    `BayesRule()` and `BayesRule2()` is that the former checks that the
    `y` values are all either `0` or `1` and that the `yhat` values are
    all between 0 and 1, while the latter doesn't and is therefore
    faster.

## Functions

- `mse()`: Mean-square error.

- `rmse()`: Root-mean-square error.

- `medAbsErr()`: Median absolute error.

- `BayesRule()`: Bayes Rule for a binary response.

- `BayesRule2()`: Bayes rule for a binary response (without bounds
  checking).

## See also

[`cv`](https://gmonette.github.io/cv/reference/cv.md),
[`cv.merMod`](https://gmonette.github.io/cv/reference/cv.merMod.md),
[`cv.function`](https://gmonette.github.io/cv/reference/cv.function.md).

## Examples

``` r
if (requireNamespace("carData", quietly=TRUE)){
withAutoprint({
data("Duncan", package="carData")
m.lm <- lm(prestige ~ income + education, data=Duncan)
mse(Duncan$prestige, fitted(m.lm))

data("Mroz", package="carData")
m.glm <- glm(lfp ~ ., data=Mroz, family=binomial)
BayesRule(Mroz$lfp == "yes", fitted(m.glm))
})
} else {
cat("\n install 'carData' package to run these examples\n")
}
#> > data("Duncan", package = "carData")
#> > m.lm <- lm(prestige ~ income + education, data = Duncan)
#> > mse(Duncan$prestige, fitted(m.lm))
#> [1] 166.8155
#> attr(,"casewise loss")
#> [1] "(y - yhat)^2"
#> > data("Mroz", package = "carData")
#> > m.glm <- glm(lfp ~ ., data = Mroz, family = binomial)
#> > BayesRule(Mroz$lfp == "yes", fitted(m.glm))
#> [1] 0.3067729
#> attr(,"casewise loss")
#> [1] "y != round(yhat)"
```
