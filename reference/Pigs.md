# Body Weights of 48 Pigs in 9 Successive Weeks

This data set appears in Table 3.1 of Diggle, Liang, and Zeger (1994).

## Usage

``` r
data("Pigs", package = "cv")
```

## Format

A data frame with 432 rows and 3 columns.

- id:

  Pig id number, 1–48.

- week:

  Week number, 1–9.

- weight:

  Weight in kg.

## Source

P. J. Diggle, K.-Y. Liang, and S. L. Zeger, *Analysis of Longitudinal
Data* (Oxford, 1994).

## Examples

``` r
library("lme4")
#> Loading required package: Matrix
m.p <- lmer(weight ~ week + (1 | id) + (1 | week),
            data=Pigs, REML=FALSE,
            control=lmerControl(optimizer="bobyqa"))
summary(m.p)
#> Linear mixed model fit by maximum likelihood  ['lmerMod']
#> Formula: weight ~ week + (1 | id) + (1 | week)
#>    Data: Pigs
#> Control: lmerControl(optimizer = "bobyqa")
#> 
#>       AIC       BIC    logLik -2*log(L)  df.resid 
#>    2037.6    2058.0   -1013.8    2027.6       427 
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.7750 -0.5418  0.0054  0.4762  3.9816 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev.
#>  id       (Intercept) 14.83622 3.8518  
#>  week     (Intercept)  0.08499 0.2915  
#>  Residual              4.29733 2.0730  
#> Number of obs: 432, groups:  id, 48; week, 9
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept) 19.35561    0.63340   30.56
#> week         6.20990    0.05393  115.14
#> 
#> Correlation of Fixed Effects:
#>      (Intr)
#> week -0.426
cv(m.p, clusterVariables=c("id", "week"), k=10, seed=8469)
#> R RNG seed set to 8469
#> cross-validation criterion (mse) = 19.2352 
```
