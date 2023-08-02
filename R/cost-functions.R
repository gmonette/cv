#' Cost Functions for Fitted Regression Models
#'
#' Compute cost functions (cross-validation criteria) for fitted
#' regression models.
#'
#' @param y response
#' @param yhat fitted value
#'
#' @describeIn cost-functions Mean-Square Error
#' @export
mse <- function(y, yhat){
  mean((y - yhat)^2)
}

#' @describeIn cost-functions Bayes for a binary response
#' @export
BayesRule <- function(y, yhat){
  if (!all(y %in% c(0, 1))) stop("response values not all 0 or 1")
  if (any(yhat < 0) || any(yhat > 1)) stop("fitted values outside of interval [0, 1]")
  yhat <- round(yhat)
  mean(y != yhat) # proportion in error
}

#' @describeIn cost-functions Bayes for a binary response (without bounds checking)
#' @export
BayesRule2 <- function(y, yhat){
  yhat <- round(yhat)
  mean(y != yhat) # proportion in error
}
