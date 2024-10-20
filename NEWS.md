# cv 2.1.0

- New Arima() function, provides formula interface to standard-R arima() function.

- New cvOrdered() to support timeseries models, used by new cv.ARIMA() method.

- Corresponding changes to folds() and  fold() to support timeseries data.

# cv 2.0.3

- New plot.cv() and plot.cvList() methods.

- New cvInfo() accessor function with "cv", "cvList", "cvModList", and "cvSelect" methods.

- Differentiate print() and summary() methods for "cv", "cvList", and "cvModlist" objects.

- Fixes to computing per-fold details for mixed-models cv() methods.

- Rename "recursive CV" as "meta CV" and edit functions, arguments, examples, etc., to reflect this change.

- Small fixes.

# cv 2.0.2

- changed cache options for cv-mixed vignette

# cv 2.0.1

- New examples for cross-validation with mixed models.

- Updated GetResponse.glmmTMB() method.

- Small fix to docs.

- Small improvements.

# cv 2.0.0

- New cv.function() method meant to replace cvSelect(), direct use of which is now discouraged.

- New selectModelList() to be used with cv.function() (or with cvSelect()). selectModelList() implements recursive cross-validation, where the fit of a model selected by CV is assessed by CV. The same procedure is also available by setting recursive=TRUE in a call to cv.modList().

- cv.default() and other cv() methods acquire a details argument, which if TRUE includes information about the folds in the returned object.

- New as.data.frame.cv() and related methods for turning the detailed results returned by cv() methods into a data frame, with new print() and summary() methods for the objects produced.

- Improvements to code, introducing folds(), fold(), and related functions.

- Refactoring of code; cv() methods now all call cvCompute() (which is new), cvMixed(), or cvSelect().

- Reorganization of package file structure and of documentation.

- Make the cv.default() method more robust, particularly for parallel computations.

- Reorganize package vignettes (of which there are now 5).

- Other small improvements.

# cv 1.1.0

- cv() et al. now work properly with "non-casewise average" CV criteria such as the new rmse() and medAbsErr(), not just with "casewise-average" fit criteria such as mse() and BayesRule().

- Bias adjustment and confidence intervals (which are new) are computed only for casewise-average CV criteria. Demonstrate that 1 - AUC isn't a casewise-average criterion.

- Generally suppress spurious messages about setting the seed in cv.modList() for LOO CV.

- Fix bugs in selectTrans() that caused errors when one of response and predictors arguments not specified.

- Fix bug in cvMixed() that prevented parallel computations (reported by Craig See).

- Fix small bug in cvSelect(), returning properly named "coefficients" element when save.coef is TRUE.

- Fix bug in cv.lm() and cv.glm() with method="hatvalues" for cost criteria other than mse().

- Add selectTransStepAIC() procedure for use with cvSelect().

- Add medAbsErr() and rmse() cost criteria.

- Add coef.cvSelect() method.

- Add cv.rlm() method.

- plot.cvModList() can show averages +/- SDs, and averages and CIs, as well as averages and ranges.

- Add Pigs data set.

- change getResponse() and methods to GetResponse() to avoid name clash with nlme.

- Improvements and updates to documentation, and expanded cv.Rmd vignette.

- Mixed-models methods no longer flagged as "experimental."

- Mixed-models CV functions no longer limited to nested random effects.

# cv 1.0.1

- Initial CRAN version.

# cv 0.1.0

- Initial version.
