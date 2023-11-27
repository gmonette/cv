# cv 1.1.0

- Generally suppress spurious messages about setting the seed in cv.modList() for LOO CV.

- Fix bugs in selectTrans() that caused errors when one of response and predictors arguments not specified.

- Fix bug in cvMixed() that prevented parallel computations (reported by Craig See).

- Fix small bug in cvSelect(), returning properly named "coefficients" element when save.coef is TRUE.

- Add selectTransAndStepAIC() procedure for use with cvSelect().

- Add medAbsErr() and rmse() cost criteria.

- Add coef.cvSelect() method.

- Add cv.rlm() method.

- Add Pigs data set.

- Improvements to documentation and expanded cv.Rmd vignette.

- Mixed-models methods no longer flagged as "experimental."

- Mixed-models CV functions no longer limited to nested random effects.

# cv 1.0.1

- Initial CRAN version.

# cv 0.1.0

- Initial version.
