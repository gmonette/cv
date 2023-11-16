# cv 1.0.2

- Generally suppress spurious messages about setting the seed in cv.modList() for LOO CV.

- Fix bugs in selectTrans() that caused errors when one of response and predictors arguments not specified.

- Fix bug in cvMixed() that prevented parallel computations (reported by Craig See).

- Fix small bug in cvSelect(), returning properly named "coefficients" element when save.coef is TRUE.

- Add selectTransAndStepAIC() procedure for use with cvSelect().

- Add medAbsErr() cost criterion.

- Add coef.cvSelect() method.

- Improvements to documentation and expanded cv.Rmd vignette.

# cv 1.0.1

- Initial CRAN version.

# cv 0.1.0

- Initial version.
