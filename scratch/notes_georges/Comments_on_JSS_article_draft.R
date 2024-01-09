#'
#' Pages with reference to Jan 6, JSS-article-reduced.pdf
#'
#'
#'
#' - Edit a number of references to 'this vignette'
#' Typos
#' - \code{xm] to \code{xm}
#' - \code{women] to \code{women}
#' - 'a advantage' to 'an advantage'
#' - 'here we specify family="yjPower' (missing final ")
#'
#' p. 8
#'
#' Can keep related objects in a list
#' and avoid typing out 'm.1, m.2, ...'
#'
#' Also, using
#' 'do.call(models, mlist)'
#' uses the list names for the models.
#'
mlist <- list()
for (p in 1:10) mlist[[p]] <- lm(mpg ~ poly(horsepower, p), data = Auto)
names(mlist) <- paste0("m.", 1:10)

summary(mlist[['m.2']])  # for example, the quadratic fit


mlist <- do.call(models, mlist)

# 10-fold CV
cv.auto.10 <- cv(mlist, data = Auto, seed = 2120)

cv.auto.10[1:2]  # note that the names created for the list are used for the models

#'
#'
#
#'
#' p.12:
#'
#' We can avoid the complexity of summarizing with dplyr by using base::ave as on
#' p. 17.  The process seems easier and more transparent using base tools
#' than dplyr.
#'
library(cv)
library(ISLR2)
library(nlme)

HSB <- MathAchieve
names(HSB) <- tolower(names(HSB))

HSB <-
  within(HSB,
    {
      mean.ses <- ave(ses, school)
      cses <- ses - mean.ses
    }
  )

#'
#' p. 52:
#'
#' "This approach can break down when one more hatvalues are equal to 1, ...
#' requires division by 0."
#'
#' could add:
#'
#' ---
#' In this case the training set omitting the observation with hatvalue = 1
#' is rank-deficient and the predictors for the case left out
#' are outside the linear span of the predictors in the training set.
#' ---
#'
#'
#'
#'

