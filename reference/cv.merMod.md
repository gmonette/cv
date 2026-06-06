# Cross-Validate Mixed-Effects Model

[`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods for
mixed-effect models of class `"merMod"`, fit by the
[`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) and
[`glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) functions in the
lme4 package; for models of class `"lme"` fit by the
[`lme()`](https://rdrr.io/pkg/nlme/man/lme.html) function in the nlme
package; and for models of class `"glmmTMB"` fit by the
[`glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html) function in
the glmmTMB package.

## Usage

``` r
# S3 method for class 'merMod'
cv(
  model,
  data = insight::get_data(model),
  criterion = mse,
  k = NULL,
  reps = 1L,
  seed,
  details = NULL,
  ncores = 1L,
  clusterVariables,
  ...
)

# S3 method for class 'lme'
cv(
  model,
  data = insight::get_data(model),
  criterion = mse,
  k = NULL,
  reps = 1L,
  seed,
  details = NULL,
  ncores = 1L,
  clusterVariables,
  ...
)

# S3 method for class 'glmmTMB'
cv(
  model,
  data = insight::get_data(model),
  criterion = mse,
  k = NULL,
  reps = 1L,
  seed,
  details = NULL,
  ncores = 1L,
  clusterVariables,
  ...
)
```

## Arguments

- model:

  a mixed-effects model object for which a
  [`cv()`](https://gmonette.github.io/cv/reference/cv.md) method is
  available.

- data:

  data frame to which the model was fit (not usually necessary).

- criterion:

  cross-validation ("cost" or lack-of-fit) criterion function of form
  `f(y, yhat)` where `y` is the observed values of the response and
  `yhat` the predicted values; the default is
  [`mse`](https://gmonette.github.io/cv/reference/cost-functions.md)
  (the mean-squared error).

- k:

  perform k-fold cross-validation; `k` may be a number or `"loo"` or
  `"n"` for n-fold (leave-one-out) cross-validation; the default is `10`
  if cross-validating individual cases and `"loo"` if cross-validating
  clusters.

- reps:

  number of times to replicate k-fold CV (default is `1`), or greater),
  compute a confidence interval for the bias-corrected CV criterion, if
  the criterion is the average of casewise components.

- seed:

  for R's random number generator; optional, if not supplied a random
  seed will be selected and saved; not needed for n-fold
  cross-validation.

- details:

  if `TRUE` (the default if the number of folds `k <= 10`), save
  detailed information about the value of the CV criterion for the cases
  in each fold and the regression coefficients with that fold deleted.

- ncores:

  number of cores to use for parallel computations (default is `1`,
  i.e., computations aren't done in parallel).

- clusterVariables:

  a character vector of names of the variables defining clusters for a
  mixed model with nested or crossed random effects; if missing,
  cross-validation is performed for individual cases rather than for
  clusters.

- ...:

  for [`cv()`](https://gmonette.github.io/cv/reference/cv.md) methods,
  to match generic, and for
  [`cvMixed()`](https://gmonette.github.io/cv/reference/cvCompute.md),
  arguments to be passed to
  [`update()`](https://rdrr.io/r/stats/update.html).

## Value

The methods `cv.merMod()`, `cv.lme()`, and `cv.glmmTMB()`, return
objects of class `"cv"`, or, if `reps > 1`, of class `"cvList"` (see
[`cv()`](https://gmonette.github.io/cv/reference/cv.md)).

## Details

For mixed-effects models, cross-validation can be done by "clusters" or
by individual observations. If the former, predictions are based only on
fixed effects; if the latter, predictions include the random effects
(i.e., are the best linear unbiased predictors or "BLUPS").

The model supplied must work properly with
[`update()`](https://rdrr.io/r/stats/update.html), and in particular the
formula for the model should not include variables that are not in the
data set to which the model was fit. See the last (faulty) example in
the help for [`cv()`](https://gmonette.github.io/cv/reference/cv.md).

## Functions

- `cv(merMod)`: [`cv()`](https://gmonette.github.io/cv/reference/cv.md)
  method for [`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) and
  [`glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) models from the
  lme4 package.

- `cv(lme)`: [`cv()`](https://gmonette.github.io/cv/reference/cv.md)
  method for [`lme()`](https://rdrr.io/pkg/nlme/man/lme.html) models
  from the nlme package.

- `cv(glmmTMB)`: [`cv()`](https://gmonette.github.io/cv/reference/cv.md)
  method for [`glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)
  models from the glmmTMB package.

## See also

[`cv`](https://gmonette.github.io/cv/reference/cv.md),
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html),
[`glmer`](https://rdrr.io/pkg/lme4/man/glmer.html),
[`lme`](https://rdrr.io/pkg/nlme/man/lme.html),
[`glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)

## Examples

``` r
library("lme4")
# from ?lmer:
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: Reaction ~ Days + (Days | Subject)
#>    Data: sleepstudy
#> REML criterion at convergence: 1743.628
#> Random effects:
#>  Groups   Name        Std.Dev. Corr 
#>  Subject  (Intercept) 24.741        
#>           Days         5.922   0.07 
#>  Residual             25.592        
#> Number of obs: 180, groups:  Subject, 18
#> Fixed Effects:
#> (Intercept)         Days  
#>      251.41        10.47  
summary(cv(fm1, clusterVariables="Subject")) # LOO CV of clusters
#> n-Fold Cross Validation based on 18 {Subject} clusters
#> criterion: mse
#> cross-validation criterion = 2460.604
#> bias-adjusted cross-validation criterion = 2454.627
#> full-sample criterion = 2251.398
summary(cv(fm1, seed=447)) # 10-fold CV of cases
#> R RNG seed set to 447
#> 10-Fold Cross Validation
#> criterion: mse
#> cross-validation criterion = 869.533
#> bias-adjusted cross-validation criterion = 847.5586
#> full-sample criterion = 549.342
summary(cv(fm1, clusterVariables="Subject", k=5,
   seed=834, reps=3)) # 5-fold CV of clusters, repeated 3 times
#> R RNG seed set to 834
#> R RNG seed set to 448690
#> R RNG seed set to 916534
#> 
#> Replicate 1:
#> 5-Fold Cross Validation based on 18 {Subject} clusters
#> criterion: mse
#> cross-validation criterion = 2458.192
#> bias-adjusted cross-validation criterion = 2434.117
#> full-sample criterion = 2251.398
#> 
#> Replicate 2:
#> 5-Fold Cross Validation based on 18 {Subject} clusters
#> criterion: mse
#> cross-validation criterion = 2511.651
#> bias-adjusted cross-validation criterion = 2479.989
#> full-sample criterion = 2251.398
#> 
#> Replicate 3:
#> 5-Fold Cross Validation based on 18 {Subject} clusters
#> criterion: mse
#> cross-validation criterion = 2403.755
#> bias-adjusted cross-validation criterion = 2387.049
#> full-sample criterion = 2251.398
#> 
#> Average:
#> 5-Fold Cross Validation based on 18 {Subject} clusters
#> criterion: mse
#> cross-validation criterion = 2457.947 (44.04918)
#> bias-adjusted cross-validation criterion = 2433.818 (37.94423)
#> full-sample criterion = 2251.398

library("nlme")
#> 
#> Attaching package: ‘nlme’
#> The following object is masked from ‘package:lme4’:
#> 
#>     lmList
# from ?lme
(fm2 <- lme(distance ~ age + Sex, data = Orthodont,
            random = ~ 1))
#> Linear mixed-effects model fit by REML
#>   Data: Orthodont 
#>   Log-restricted-likelihood: -218.7563
#>   Fixed: distance ~ age + Sex 
#> (Intercept)         age   SexFemale 
#>  17.7067130   0.6601852  -2.3210227 
#> 
#> Random effects:
#>  Formula: ~1 | Subject
#>         (Intercept) Residual
#> StdDev:    1.807425 1.431592
#> 
#> Number of Observations: 108
#> Number of Groups: 27 
summary(cv(fm2)) # LOO CV of cases
#> R RNG seed set to 765199
#> 10-Fold Cross Validation
#> criterion: mse
#> cross-validation criterion = 2.666103
#> bias-adjusted cross-validation criterion = 2.589541
#> full-sample criterion = 1.582435
summary(cv(fm2, clusterVariables="Subject",
        k=5, seed=321)) # 5-fold CV of clusters
#> R RNG seed set to 321
#> 5-Fold Cross Validation based on 27 {Subject} clusters
#> criterion: mse
#> cross-validation criterion = 5.875411
#> bias-adjusted cross-validation criterion = 5.780695
#> full-sample criterion = 5.017326

library("glmmTMB")
# from ?glmmTMB
(m1 <- glmmTMB(count ~ mined + (1|site),
               zi=~mined,
               family=poisson, data=Salamanders))
#> Formula:          count ~ mined + (1 | site)
#> Zero inflation:         ~mined
#> Data: Salamanders
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#> 1908.4695 1930.8080 -949.2348 1898.4695       639 
#> Random-effects (co)variances:
#> 
#> Conditional model:
#>  Groups Name        Std.Dev.
#>  site   (Intercept) 0.28    
#> 
#> Number of obs: 644 / Conditional model: site, 23
#> 
#> Fixed Effects:
#> 
#> Conditional model:
#> (Intercept)      minedno  
#>      0.0879       1.1419  
#> 
#> Zero-inflation model:
#> (Intercept)      minedno  
#>       1.139       -1.736  
summary(cv(m1, seed=97816, k=5,
          clusterVariables="site")) # 5-fold CV of clusters
#> R RNG seed set to 97816
#> 5-Fold Cross Validation based on 23 {site} clusters
#> criterion: mse
#> cross-validation criterion = 6.006117
#> bias-adjusted cross-validation criterion = 6.002191
#> 95% CI for bias-adjusted CV criterion = (2.385452, 9.61893)
#> full-sample criterion = 5.970489
summary(cv(m1, seed=34506, k=5)) # 5-fold CV of cases
#> R RNG seed set to 34506
#> 5-Fold Cross Validation
#> criterion: mse
#> cross-validation criterion = 6.058988
#> bias-adjusted cross-validation criterion = 6.032337
#> 95% CI for bias-adjusted CV criterion = (2.483579, 9.581094)
#> full-sample criterion = 5.79783
```
