# cv 1.1.0

- cv() et al. now work properly with "non-casewise average" CV criteria such as the new rmse() and medAbsErr(), not just with "casewise-average" fit criteria such as mse() and BayesRule().

- Bias adjustment and confidence intervals (which are new) are computed only for casewise-average CV criteria.

- Generally suppress spurious messages about setting the seed in cv.modList() for LOO CV.

- Fix bugs in selectTrans() that caused errors when one of response and predictors arguments not specified.

- Fix bug in cvMixed() that prevented parallel computations (reported by Craig See).

- Fix small bug in cvSelect(), returning properly named "coefficients" element when save.coef is TRUE.

- Fix bug in cv.lm() and cv.glm() with method="hatvalues" for cost criteria other than mse().

- Add selectTransAndStepAIC() procedure for use with cvSelect().

- Add medAbsErr() and rmse() cost criteria; make rmse() the default, rather than mse().

- Add coef.cvSelect() method.

- Add cv.rlm() method.

- plot.cvModList() can show averages +/- SDs as well as (the default) averages and ranges.

- Add Pigs data set.

- Improvements to documentation and expanded cv.Rmd vignette.

- Mixed-models methods no longer flagged as "experimental."

- Mixed-models CV functions no longer limited to nested random effects.

# cv 1.0.1

- Initial CRAN version.

# cv 0.1.0

- Initial version.
