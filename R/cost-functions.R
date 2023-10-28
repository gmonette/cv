#' Cost Functions for Fitted Regression Models
#'
#' Compute cost functions (cross-validation criteria) for fitted
#' regression models.
#'
#' @aliases costFunctions
#'
#' @param y response
#' @param yhat fitted value
#'
#' @details
#' Two cost functions (cross-validation criteria) are provided: (1)
#' \code{mse()} returns the mean-squared error of prediction for
#' a numeric response variable \code{y} and predictions \code{yhat}.
#' (2) \code{BayesRule()} and \code{BayesRule2()} report the proportion
#' of correct predictions for a dichotomous response variable \code{y}, assumed
#' coded \code{0} and \code{1}. The \code{yhat} values are
#' predicted probabilities and are rounded to 0 or 1. The distinction
#' between \code{BayesRule()} and \code{BayesRule2()} is that the former
#' checks that the \code{y} values are all either \code{0} or \code{1}
#' and that the \code{yhat} values are all between 0 and 1, while
#' the latter doesn't and is therefore faster.
#'
#' @returns In general, cost functions should return a single numeric
#' value measuring lack-of-fit. `mse()` returns the mean-squared error; `BayesRule()` and
#' `BayesRule2()` return the proportion of misclassified cases.

#' @describeIn cost-functions Mean-square error
#' @export
mse <- function(y, yhat){
  mean((y - yhat)^2)
}

#' @describeIn cost-functions Bayes Rule for a binary response
#' @export
BayesRule <- function(y, yhat){
  if (!all(y %in% c(0, 1))) stop("response values not all 0 or 1")
  if (any(yhat < 0) || any(yhat > 1)) stop("fitted values outside of interval [0, 1]")
  yhat <- round(yhat)
  mean(y != yhat) # proportion in error
}

#' @describeIn cost-functions Bayes rule for a binary response (without bounds checking)
#' @export
BayesRule2 <- function(y, yhat){
  yhat <- round(yhat)
  mean(y != yhat) # proportion in error
}
