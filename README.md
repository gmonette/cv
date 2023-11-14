<!-- badges: start -->

[![R-CMD-check](https://github.com/gmonette/cv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gmonette/cv/actions/workflows/R-CMD-check.yaml) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable) [![Last Commit](https://img.shields.io/github/last-commit/gmonette/cv)](https://github.com/gmonette/cv) [![CRAN](https://www.r-pkg.org/badges/version/cv)](https://cran.r-project.org/package=cv)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/cv)](https://cran.r-project.org/package=cv)
[![](https://img.shields.io/badge/pkgdown%20site-brightgreen)](https://gmonette.github.io/cv/) 

<!-- badges: end -->

# cv package for R: Various Functions for Cross-Validation of Regression Models <img src="man/figures/cv-hex.png" style="float:right; height:200px;"/>

Some of the functions supplied by the package:

-   `cv()` is a generic function with a default method and computationally efficient `"lm"` and `"glm"` methods, along with a method for a list of competing models. There are also experimental `"merMod"`, `"lme"`, and `"glmmTMB"` methods for mixed-effects models. `cv()` supports parallel computations.

-   `mse()` (mean-squared error), `medAbsErr()` (median absolute error), and `BayesRule()` are cross-validation criteria ("cost functions"), suitable for use with `cv()`.

-   `cvSelect()` cross-validates a selection procedure for a regression model. `cvSelect()` also supports parallel computations.

-   `selectStepAIC()` is a model-selection procedure, suitable for use with `cvSelect()`, based on the `stepAIC()` function in the **MASS** package.

-   `selectTrans()` is a procedure for selecting predictor and response transformations in regression, suitable for use with `cvSelect()`, based on the `powerTransform()` function in the **car** package.

-   `selectTransStepAIC()` is a procedure also suitable for use with `cvSelect()`,  that first selects predictor and response transformations and then selects predictors.

For additional information on using the **cv** package, see the "[Cross-validation of regression models](https://gmonette.github.io/cv/articles/cv.html)" vignette (`vignette("cv", package="cv")`). The **cv** package is designed to be extensible to other classes of regression models and other model-selection procedures; for details, see the "[Extending the cv package](https://gmonette.github.io/cv/articles/cv-extend.html)" vignette (`vignette("cv-extend", package="cv")`).

To install the current version of the **cv** package from CRAN:

```
install.packages("cv")
```

To install the development version of the **cv** package from GitHub:

```
if (!require(remotes)) install.packages("remotes")
remotes::install_github("gmonette/cv", build_vignettes=TRUE,
  dependencies=TRUE)
```
