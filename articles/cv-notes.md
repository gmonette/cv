# Computational and technical notes on cross-validating regression models

## Efficient computations for linear and generalized linear models

The most straightforward way to implement cross-validation in R for
statistical modeling functions that are written in the canonical manner
is to use [`update()`](https://rdrr.io/r/stats/update.html) to refit the
model with each fold removed. This is the approach taken in the default
method for [`cv()`](https://gmonette.github.io/cv/reference/cv.md), and
it is appropriate if the cases are independently sampled. Refitting the
model in this manner for each fold is generally feasible when the number
of folds in modest, but can be prohibitively costly for leave-one-out
cross-validation when the number of cases is large.

The `"lm"` and `"glm"` methods for
[`cv()`](https://gmonette.github.io/cv/reference/cv.md) take advantage
of computational efficiencies by avoiding refitting the model with each
fold removed. Consider, in particular, the weighted linear model
$`\mathbf{y}_{n \times 1} = \mathbf{X}_{n \times p}\boldsymbol{\beta}_{p \times 1} + \boldsymbol{\varepsilon}_{n \times 1}`$,
where
$`\boldsymbol{\varepsilon} \sim \mathbf{N}_n \left(\mathbf{0}, \sigma^2 \mathbf{W}^{-1}_{n \times n}\right)`$.
Here, $`\mathbf{y}`$ is the response vector, $`\mathbf{X}`$ the model
matrix, and $`\boldsymbol{\varepsilon}`$ the error vector, each for
$`n`$ cases, and $`\boldsymbol{\beta}`$ is the vector of $`p`$
population regression coefficients. The errors are assumed to be
multivariately normally distributed with 0 means and covariance matrix
$`\sigma^2 \mathbf{W}^{-1}`$, where $`\mathbf{W} = \mathrm{diag}(w_i)`$
is a diagonal matrix of inverse-variance weights. For the linear model
with constant error variance, the weight matrix is taken to be
$`\mathbf{W} = \mathbf{I}_n`$, the order-$`n`$ identity matrix.

The weighted-least-squares (WLS) estimator of $`\boldsymbol{\beta}`$ is
(see, e.g., Fox, 2016, sec. 12.2.2) [^1]
``` math
\mathbf{b}_{\mathrm{WLS}} = \left( \mathbf{X}^T \mathbf{W} \mathbf{X} \right)^{-1} 
  \mathbf{X}^T \mathbf{W} \mathbf{y}
```

Fitted values are then
$`\widehat{\mathbf{y}} = \mathbf{X}\mathbf{b}_{\mathrm{WLS}}`$.

The LOO fitted value for the $`i`$th case can be efficiently computed by
$`\widehat{y}_{-i} = y_i - e_i/(1 - h_i)`$ where
$`h_i = \mathbf{x}^T_i \left( \mathbf{X}^T \mathbf{W} \mathbf{X} \right)^{-1} \mathbf{x}_i`$
(the so-called “hatvalue”). Here, $`\mathbf{x}^T_i`$ is the $`i`$th row
of $`\mathbf{X}`$, and $`\mathbf{x}_i`$ is the $`i`$th row written as a
column vector. This approach can break down when one or more hatvalues
are equal to 1, in which case the formula for $`\widehat{y}_{-i}`$
requires division by 0.

To compute cross-validated fitted values when the folds contain more
than one case, we make use of the Woodbury matrix identity (Hager,
1989),
``` math
\left(\mathbf{A}_{m \times m} + \mathbf{U}_{m \times k} 
\mathbf{C}_{k \times k} \mathbf{V}_{k \times m} \right)^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U} \left(\mathbf{C}^{-1} + 
\mathbf{VA}^{-1}\mathbf{U} \right)^{-1} \mathbf{VA}^{-1}
```
where $`\mathbf{A}`$ is a nonsingular order-$`n`$ matrix. We apply this
result by letting
``` math
\begin{align*}
    \mathbf{A} &= \mathbf{X}^T \mathbf{W} \mathbf{X} \\
    \mathbf{U} &= \mathbf{X}_\mathbf{j}^T \\
    \mathbf{V} &= - \mathbf{X}_\mathbf{j} \\
    \mathbf{C} &= \mathbf{W}_\mathbf{j} \\
\end{align*}
```
where the subscript $`\mathbf{j} = (i_{j1}, \ldots, i_{jm})^T`$
represents the vector of indices for the cases in the $`j`$th fold,
$`j = 1, \ldots, k`$. The negative sign in
$`\mathbf{V} = - \mathbf{X}_\mathbf{j}`$ reflects the *removal*, rather
than addition, of the cases in $`\mathbf{j}`$.

Applying the Woodbury identity isn’t quite as fast as using the
hatvalues, but it is generally much faster than refitting the model. A
disadvantage of the Woodbury identity, however, is that it entails
explicit matrix inversion and thus may be numerically unstable. The
inverse of $`\mathbf{A} = \mathbf{X}^T \mathbf{W} \mathbf{X}`$ is
available directly in the `"lm"` object, but the second term on the
right-hand side of the Woodbury identity requires a matrix inversion
with each fold deleted. (In contrast, the inverse of each
$`\mathbf{C} = \mathbf{W}_\mathbf{j}`$ is straightforward because
$`\mathbf{W}`$ is diagonal.)

The Woodbury identity also requires that the model matrix be of full
rank. We impose that restriction in our code by removing redundant
regressors from the model matrix for all of the cases, but that doesn’t
preclude rank deficiency from surfacing when a fold is removed. Rank
deficiency of $`\mathbf{X}`$ doesn’t disqualify cross-validation because
all we need are fitted values under the estimated model.

[`glm()`](https://rdrr.io/r/stats/glm.html) computes the
maximum-likelihood estimates for a generalized linear model by iterated
weighted least squares (see, e.g., Fox & Weisberg, 2019, sec. 6.12). The
last iteration is therefore just a WLS fit of the “working response” on
the model matrix using “working weights.” Both the working weights and
the working response at convergence are available from the information
in the object returned by [`glm()`](https://rdrr.io/r/stats/glm.html).

We then treat re-estimation of the model with a case or cases deleted as
a WLS problem, using the hatvalues or the Woodbury matrix identity. The
resulting fitted values for the deleted fold aren’t exact—that is,
except for the Gaussian family, the result isn’t identical to what we
would obtain by literally refitting the model—but in our (limited)
experience, the approximation is very good, especially for LOO CV, which
is when we would be most tempted to use it. Nevertheless, because these
results are approximate, the default for the `"glm"`
[`cv()`](https://gmonette.github.io/cv/reference/cv.md) method is to
perform the exact computation, which entails refitting the model with
each fold omitted.

Let’s compare the efficiency of the various computational methods for
linear and generalized linear models. Consider, for example,
leave-one-out cross-validation for the quadratic regression of `mpg` on
`horsepower` in the `Auto` data, from the introductory “Cross-validating
regression models” vignette, repeated here:

``` r

data("Auto", package="ISLR2")
m.auto <- lm(mpg ~ poly(horsepower, 2), data = Auto)
summary(m.auto)
#> 
#> Call:
#> lm(formula = mpg ~ poly(horsepower, 2), data = Auto)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -14.714  -2.594  -0.086   2.287  15.896 
#> 
#> Coefficients:
#>                      Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)            23.446      0.221   106.1   <2e-16 ***
#> poly(horsepower, 2)1 -120.138      4.374   -27.5   <2e-16 ***
#> poly(horsepower, 2)2   44.090      4.374    10.1   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 4.37 on 389 degrees of freedom
#> Multiple R-squared:  0.688,  Adjusted R-squared:  0.686 
#> F-statistic:  428 on 2 and 389 DF,  p-value: <2e-16

library("cv")
#> Loading required package: doParallel
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
summary(cv(m.auto, k = "loo") ) # default method = "hatvalues"
#> n-Fold Cross Validation
#> method: hatvalues
#> criterion: mse
#> cross-validation criterion = 19.248
summary(cv(m.auto, k = "loo", method = "naive"))
#> n-Fold Cross Validation
#> method: naive
#> criterion: mse
#> cross-validation criterion = 19.248
#> bias-adjusted cross-validation criterion = 19.248
#> full-sample criterion = 18.985
summary(cv(m.auto, k = "loo", method = "Woodbury"))
#> n-Fold Cross Validation
#> method: Woodbury
#> criterion: mse
#> cross-validation criterion = 19.248
#> bias-adjusted cross-validation criterion = 19.248
#> full-sample criterion = 18.985
```

This is a small regression problem and all three computational
approaches are essentially instantaneous, but it is still of interest to
investigate their relative speed. In this comparison, we include the
[`cv.glm()`](https://gmonette.github.io/cv/reference/cv.md) function
from the **boot** package (Canty & Ripley, 2022; Davison & Hinkley,
1997), which takes the naive approach, and for which we have to fit the
linear model as an equivalent Gaussian GLM. We use the
`microbenchmark()` function from the package of the same name for the
timings (Mersmann, 2023):

``` r

m.auto.glm <- glm(mpg ~ poly(horsepower, 2), data = Auto)
boot::cv.glm(Auto, m.auto.glm)$delta
#> [1] 19.248 19.248

microbenchmark::microbenchmark(
  hatvalues = cv(m.auto, k = "loo"),
  Woodbury = cv(m.auto, k = "loo", method = "Woodbury"),
  naive = cv(m.auto, k = "loo", method = "naive"),
  cv.glm = boot::cv.glm(Auto, m.auto.glm),
  times = 10
)
#> Warning in microbenchmark::microbenchmark(hatvalues = cv(m.auto, k = "loo"), :
#> less accurate nanosecond times to avoid potential integer overflows
#> Unit: milliseconds
#>       expr      min       lq     mean   median       uq      max neval
#>  hatvalues   2.0507   2.2768   2.6806   2.3058   2.7121   5.4399    10
#>   Woodbury  13.5814  14.8696  19.0034  16.4796  24.1375  27.7238    10
#>      naive 234.6956 258.4918 278.0291 289.1408 297.0773 301.7305    10
#>     cv.glm 450.2114 462.6246 537.8314 510.7763 607.9489 745.1678    10
```

On our computer, using the hatvalues is about an order of magnitude
faster than employing Woodbury matrix updates, and more than two orders
of magnitude faster than refitting the model.[^2]

Similarly, let’s return to the logistic-regression model fit to Mroz’s
data on women’s labor-force participation, also employed as an example
in the introductory vignette:

``` r

data("Mroz", package="carData")
m.mroz <- glm(lfp ~ ., data = Mroz, family = binomial)
summary(m.mroz)
#> 
#> Call:
#> glm(formula = lfp ~ ., family = binomial, data = Mroz)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.18214    0.64438    4.94  7.9e-07 ***
#> k5          -1.46291    0.19700   -7.43  1.1e-13 ***
#> k618        -0.06457    0.06800   -0.95  0.34234    
#> age         -0.06287    0.01278   -4.92  8.7e-07 ***
#> wcyes        0.80727    0.22998    3.51  0.00045 ***
#> hcyes        0.11173    0.20604    0.54  0.58762    
#> lwg          0.60469    0.15082    4.01  6.1e-05 ***
#> inc         -0.03445    0.00821   -4.20  2.7e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 1029.75  on 752  degrees of freedom
#> Residual deviance:  905.27  on 745  degrees of freedom
#> AIC: 921.3
#> 
#> Number of Fisher Scoring iterations: 4

summary(cv(m.mroz, # default method = "exact"
   k = "loo", 
   criterion = BayesRule))
#> n-Fold Cross Validation
#> method: exact
#> criterion: BayesRule
#> cross-validation criterion = 0.32005
#> bias-adjusted cross-validation criterion = 0.3183
#> 95% CI for bias-adjusted CV criterion = (0.28496, 0.35164)
#> full-sample criterion = 0.30677
summary(cv(m.mroz,
   k = "loo",
   criterion = BayesRule,
   method = "Woodbury"))
#> n-Fold Cross Validation
#> method: Woodbury
#> criterion: BayesRule
#> cross-validation criterion = 0.32005
#> bias-adjusted cross-validation criterion = 0.3183
#> 95% CI for bias-adjusted CV criterion = (0.28496, 0.35164)
#> full-sample criterion = 0.30677
summary(cv(m.mroz,
   k = "loo",
   criterion = BayesRule,
   method = "hatvalues"))
#> n-Fold Cross Validation
#> method: hatvalues
#> criterion: BayesRule
#> cross-validation criterion = 0.32005
```

As for linear models, we report some timings for the various
[`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods of
computation in LOO CV as well as for the
[`cv.glm()`](https://gmonette.github.io/cv/reference/cv.md) function
from the **boot** package (which, recall, refits the model with each
case removed, and thus is comparable to
[`cv()`](https://gmonette.github.io/cv/reference/cv.md) with
`method="exact"`):

``` r

microbenchmark::microbenchmark(
  hatvalues = cv(
    m.mroz,
    k = "loo",
    criterion = BayesRule,
    method = "hatvalues"
  ),
  Woodbury = cv(
    m.mroz,
    k = "loo",
    criterion = BayesRule,
    method = "Woodbury"
  ),
  exact = cv(m.mroz, k = "loo", criterion = BayesRule),
  cv.glm = boot::cv.glm(Mroz, m.mroz,
                        cost = BayesRule),
  times = 10
)
#> Unit: milliseconds
#>       expr       min        lq      mean    median        uq       max neval
#>  hatvalues    2.1878    2.4001    2.6145    2.5735    2.8389    3.1786    10
#>   Woodbury   48.3798   51.2279   52.9838   52.8029   54.7773   58.1968    10
#>      exact 2184.4771 2260.2463 2378.7911 2367.4649 2456.2159 2662.5837    10
#>     cv.glm 2442.4449 2523.9813 2693.9351 2736.7670 2772.3282 2949.3461    10
```

There is a substantial time penalty associated with exact computations.

## Computation of the bias-corrected CV criterion and confidence intervals

Let $`\mathrm{CV}(\mathbf{y}, \widehat{\mathbf{y}})`$ represent a
cross-validation cost criterion, such as mean-squared error, computed
for all of the $`n`$ values of the response $`\mathbf{y}`$ based on
fitted values $`\widehat{\mathbf{y}}`$ from the model fit to all of the
data. We require that $`\mathrm{CV}(\mathbf{y}, \widehat{\mathbf{y}})`$
is the mean of casewise components, that is,
$`\mathrm{CV}(\mathbf{y}, \widehat{\mathbf{y}}) = \frac{1}{n}\sum_{i=1}^n\mathrm{cv}(y_i, \widehat{y}_i)`$.[^3]
For example,
$`\mathrm{MSE}(\mathbf{y}, \widehat{\mathbf{y}}) = \frac{1}{n}\sum_{i=1}^n (y_i - \widehat{y}_i)^2`$.

We divide the $`n`$ cases into $`k`$ folds of approximately
$`n_j \approx n/k`$ cases each, where $`n = \sum n_j`$. As above, let
$`\mathbf{j}`$ denote the indices of the cases in the $`j`$th fold.

Now define
$`\mathrm{CV}_j = \mathrm{CV}(\mathbf{y}, \widehat{\mathbf{y}}^{(j)})`$.
The superscript $`(j)`$ on $`\widehat{\mathbf{y}}^{(j)}`$ represents
fitted values computed for all of the cases from the model with fold
$`j`$ omitted. Let $`\widehat{\mathbf{y}}^{(-i)}`$ represent the vector
of fitted values for all $`n`$ cases where the fitted value for the
$`i`$th case is computed from the model fit with the fold including the
$`i`$th case omitted (i.e., fold $`j`$ for which $`i \in \mathbf{j}`$).

Then the cross-validation criterion is just
$`\mathrm{CV} = \mathrm{CV}(\mathbf{y}, \widehat{\mathbf{y}}^{(-i)})`$.
Following Davison & Hinkley (1997, pp. 293–295), the bias-adjusted
cross-validation criterion is
``` math
\mathrm{CV}_{\mathrm{adj}} = \mathrm{CV} + \mathrm{CV}(\mathbf{y}, \widehat{\mathbf{y}}) - \frac{1}{n} \sum_{j=1}^{k} n_j \mathrm{CV}_j
```

We compute the standard error of CV as
``` math
\mathrm{SE}(\mathrm{CV}) = \frac{1}{\sqrt n} \sqrt{ \frac{\sum_{i=1}^n \left[ \mathrm{cv}(y_i, \widehat{y}_i^{(-i)} ) - \mathrm{CV} \right]^2 }{n - 1} }
```
that is, as the standard deviation of the casewise components of CV
divided by the square-root of the number of cases.

We then use $`\mathrm{SE}(\mathrm{CV})`$ to construct a
$`100 \times (1 - \alpha)`$% confidence interval around the *adjusted*
CV estimate of error:
``` math
\left[ \mathrm{CV}_{\mathrm{adj}} - z_{1 - \alpha/2}\mathrm{SE}(\mathrm{CV}), \mathrm{CV}_{\mathrm{adj}} + z_{1 - \alpha/2}\mathrm{SE}(\mathrm{CV})  \right]
```
where $`z_{1 - \alpha/2}`$ is the $`1 - \alpha/2`$ quantile of the
standard-normal distribution (e.g, $`z \approx 1.96`$ for a 95%
confidence interval, for which $`1 - \alpha/2 = .975`$).

Bates, Hastie, & Tibshirani (2023) show that the coverage of this
confidence interval is poor for small samples, and they suggest a much
more computationally intensive procedure, called *nested
cross-validation*, to compute better estimates of error and confidence
intervals with better coverage for small samples. We may implement Bates
et al.’s approach in a later release of the **cv** package. At present
we use the confidence interval above for sufficiently large $`n`$,
which, based on Bates et al.’s results, we take by default to be
$`n \ge 400`$.

## Why the complement of AUC isn’t a casewise CV criterion

Consider calculating AUC for folds in which a validation set contains
$`n_v`$ observations. To calculate AUC in the validation set, we need
the vector of prediction criteria,
$`\widehat{\mathbf{y}}_{v_{(n_v \times 1)}} = (\widehat{y}_1, ..., \widehat{y}_{n_v})^T`$,
and the vector of observed responses in the validation set,
$`\mathbf{y}_{v_{(n_v \times 1)}} = (y_1, \ldots, y_{n_v})^T`$ with
$`y_i \in \{0,1\}, \; i = 1, \ldots, n_v`$.

To construct the ROC curve, only the ordering of the values in
$`\mathbf{\widehat{y}}_v`$ is relevant. Thus, assuming that there are no
ties, and reordering observations if necessary, we can set
$`\mathbf{\widehat{y}}_v = (1, 2, \ldots, n_v)^T`$.

If the AUC can be expressed as the casewise mean or sum of a function
$`\mathrm{cv}(\widehat{y}_i,y_i)`$, where
$`\mathrm{cv}: \{1,2,...,n_v\}\times\{0,1\} \rightarrow [0,1]`$, then
``` math
\begin{equation}
\label{eq:cw}
\tag{1}
\sum_{i=1}^{n_v} \mathrm{cv}(\widehat{y}_i,y_i) = \mathrm{AUC}(\mathbf{\widehat{y}}_v,\mathbf{y}_v)
\end{equation}
```
must hold for all $`2^{n_v}`$ possible values of
$`\mathbf{y}_v = (y_1,...,y_{n_v})^T`$. If all $`y\mathrm{s}`$ have the
same value, either 1 or 0, then the definition of AUC is ambiguous. AUC
could be considered undefined, or it could be set to 0 if all $`y`$s are
0 and to 1 if all $`y`$s are 1. If AUC is considered to be undefined in
these cases, we have $`2^{n_v} - 2`$ admissible values for
$`\mathbf{y}_v`$.

Thus, equation () produces either $`2^{n_v}`$ or $`2^{n_v}-2`$
constraints. Although there are only $`2n_v`$ possible values for the
$`\mathrm{cv(\cdot)}`$ function, equation () could, nevertheless, have
consistent solutions. We therefore need to determine whether there is a
value of $`n_v`$ for which () has no consistent solution for all
admissible values of $`\mathbf{y}_v`$. In that eventuality, we will have
shown that AUC cannot, in general, be expressed through a casewise sum.

If $`n_v=3`$, we show below that () has no consistent solution if we
include all possibilities for $`\mathbf{y}_v`$, but does if we exclude
cases where all $`y`$s have the same value. If $`n_v=4`$, we show that
there are no consistent solutions in either case.

The following R function computes AUC from $`\mathbf{\widehat{y}}_v`$
and $`\mathbf{y}_v`$, accommodating the cases where $`\mathbf{y}_v`$ is
all 0s or all 1s:

``` r

AUC <- function(y, yhat = seq_along(y)) {
  s <- sum(y)
  if (s == 0)
    return(0)
  if (s == length(y))
    return(1)
  Metrics::auc(y, yhat)
}
```

We then define a function to generate all possible $`\mathbf{y}_v`$s of
length $`n_v`$ as rows of the matrix
$`\mathbf{Y}_{(2^{n_v} \times n_v)}`$:

``` r

Ymat <- function(n_v, exclude_identical = FALSE) {
  stopifnot(n_v > 0 &&
              round(n_v) == n_v)    # n_v must be a positive integer
  ret <- sapply(0:(2 ^ n_v - 1),
                function(x)
                  as.integer(intToBits(x)))[1:n_v,]
  ret <- if (is.matrix(ret))
    t(ret)
  else
    matrix(ret)
  colnames(ret) <- paste0("y", 1:ncol(ret))
  if (exclude_identical)
    ret[-c(1, nrow(ret)),]
  else
    ret
}
```

For $`n_v=3`$,

``` r

Ymat(3)
#>      y1 y2 y3
#> [1,]  0  0  0
#> [2,]  1  0  0
#> [3,]  0  1  0
#> [4,]  1  1  0
#> [5,]  0  0  1
#> [6,]  1  0  1
#> [7,]  0  1  1
#> [8,]  1  1  1
```

If we exclude $`\mathbf{y}_v`$s with identical values, then

``` r

Ymat(3, exclude_identical = TRUE)
#>      y1 y2 y3
#> [1,]  1  0  0
#> [2,]  0  1  0
#> [3,]  1  1  0
#> [4,]  0  0  1
#> [5,]  1  0  1
#> [6,]  0  1  1
```

Here is $`\mathbf{Y}`$ with corresponding values of AUC:

``` r

cbind(Ymat(3), AUC = apply(Ymat(3), 1, AUC))
#>      y1 y2 y3 AUC
#> [1,]  0  0  0 0.0
#> [2,]  1  0  0 0.0
#> [3,]  0  1  0 0.5
#> [4,]  1  1  0 0.0
#> [5,]  0  0  1 1.0
#> [6,]  1  0  1 0.5
#> [7,]  0  1  1 1.0
#> [8,]  1  1  1 1.0
```

The values of $`\mathrm{cv}(\widehat{y}_i, y_i)`$ that express AUC as a
sum of casewise values are solutions of equation (), which can be
written as solutions of the following system of $`2^{n_v}`$ linear
simultaneous equations in $`2n_v`$ unknowns:
``` math
\begin{equation}
\label{eq:lin}
\tag{2}
(\mathbf{U} -\mathbf{Y}) \mathbf{c}_0 + \mathbf{Y} \mathbf{c}_1
=
[\mathbf{U} -\mathbf{Y}, \mathbf{Y}]
\begin{bmatrix}
\mathbf{c}_0 \\ \mathbf{c}_1
\end{bmatrix}
= \mathrm{AUC}(\mathbf{\widehat{Y}},\mathbf{Y})
\end{equation}
```
where $`\mathbf{U}_{(2^{n_v} \times n_v)}`$ is a matrix of 1s
conformable with $`\mathbf{Y}`$;
$`\mathbf{c}_0 = [\mathrm{cv}(1,0), c(2,0), ..., \mathrm{cv}(n_v,0)]^T`$;
$`\mathbf{c}_1 = [\mathrm{cv}(1,1), c(2,1), ..., \mathrm{cv}(n_v,1)]^T`$;
$`[\mathbf{U} -\mathbf{Y}, \mathbf{Y}]_{(2^{n_v} \times 2n_v)}`$ and
$`\begin{bmatrix}\begin{aligned}
\mathbf{c}_0 \\ \mathbf{c}_1
\end{aligned}
\end{bmatrix}_{(2n_v \times 1)}`$ are partitioned matrices; and
$`\mathbf{\widehat{Y}}_{(2^{n_v} \times n_v)}`$ is a matrix each of
whose rows consists of the integers 1 to $`n_v`$.

We can test whether equation () has a solution for any given $`n_v`$ by
trying to solve it as a least-squares problem, considering whether the
residuals of the associated linear model are all 0, using the “design
matrix” $`[\mathbf{U} -\mathbf{Y}, \mathbf{Y}]`$ to predict the
“outcome”
$`\mathrm{AUC}(\mathbf{\widehat{Y}},\mathbf{Y})_{(2^{n_v} \times 1)}`$:

``` r

resids <- function(n_v,
                   exclude_identical = FALSE,
                   tol = sqrt(.Machine$double.eps)) {
  Y <- Ymat(n_v, exclude_identical = exclude_identical)
  AUC <- apply(Y, 1, AUC)
  X <- cbind(1 - Y, Y)
  opts <- options(warn = -1)
  on.exit(options(opts))
  fit <- lsfit(X, AUC, intercept = FALSE)
  ret <- max(abs(residuals(fit)))
  if (ret < tol) {
    ret <- 0
    solution <- coef(fit)
    names(solution) <- paste0("c(", c(1:n_v, 1:n_v), ",",
                              rep(0:1, each = n_v), ")")
    attr(ret, "solution") <- zapsmall(solution)
  }
  ret
}
```

The case $`n_v=3`$, excluding identical $`y`$s, has a solution:

``` r

resids(3, exclude_identical = TRUE)
#> [1] 0
#> attr(,"solution")
#> c(1,0) c(2,0) c(3,0) c(1,1) c(2,1) c(3,1) 
#>    1.0    0.0   -0.5    0.5    0.0    0.0
```

But, if identical $`y`$s are included, the equation is not consistent:

``` r

resids(3, exclude_identical = FALSE)
#> [1] 0.125
```

For $`n_v=4`$, there are no solutions in either case:

``` r

resids(4, exclude_identical = TRUE)
#> [1] 0.083333
resids(4, exclude_identical = FALSE)
#> [1] 0.25
```

Consequently, the widely employed AUC measure of fit for binary
regression cannot in general be used for a casewise cross-validation
criterion.

## References

Arlot, S., & Celisse, A. (2010). A survey of cross-validation procedures
for model selection. *Statistics Surveys*, *4*, 40–79. Retrieved from
<https://doi.org/10.1214/09-SS054>

Bates, S., Hastie, T., & Tibshirani, R. (2023). Cross-validation: What
does it estimate and how well does it do it? *Journal of the American
Statistical Association*, *in press*. Retrieved from
<https://doi.org/10.1080/01621459.2023.2197686>

Canty, A., & Ripley, B. D. (2022). *Boot: Bootstrap R (S-plus)
functions*.

Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap methods and their
applications*. Cambridge: Cambridge University Press.

Fox, J. (2016). *Applied regression analysis and generalized linear
models* (Second edition). Thousand Oaks CA: Sage.

Fox, J., & Weisberg, S. (2019). *An R companion to applied regression*
(Third edition). Thousand Oaks CA: Sage.

Hager, W. W. (1989). Updating the inverse of a matrix. *SIAM Review*,
*31*(2), 221–239.

Mersmann, O. (2023). *Microbenchmark: Accurate timing functions*.
Retrieved from <https://CRAN.R-project.org/package=microbenchmark>

[^1]: This is a definitional formula, which assumes that the model
    matrix $`\mathbf{X}`$ is of full column rank, and which can be
    subject to numerical instability when $`\mathbf{X}`$ is
    ill-conditioned. [`lm()`](https://rdrr.io/r/stats/lm.html) uses the
    singular-value decomposition of the model matrix to obtain
    computationally more stable results.

[^2]: Out of impatience, we asked `microbenchmark()` to execute each
    command only 10 times rather than the default 100. With the
    exception of the last columns, the output is self-explanatory. The
    last column shows which methods have average timings that are
    statistically distinguishable. Because of the small number of
    repetitions (i.e., 10), the `"hatvalues"` and `"Woodbury"` methods
    aren’t distinguishable, but the difference between these methods
    persists when we perform more repetitions—we invite the reader to
    redo this computation with the default `times=100` repetitions.

[^3]: Arlot & Celisse (2010) term the casewise loss,
    $`\mathrm{cv}(y_i, \widehat{y}_i)`$, the “contrast function.”
