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
#' @seealso \code{\link{cv}}, \code{\link{cv.merMod}},
#' \code{\link{cv.function}}.
#'
#' @details
#' Cost functions (cross-validation criteria) are meant to measure lack-of-fit. Several cost functions are provided:
#' 1. \code{mse()} returns the mean-squared error of prediction for
#' a numeric response variable \code{y} and predictions \code{yhat}; and
#' \code{rmse()} returns the root-mean-squared error and is just the
#' square-root of \code{mse()}.
#' 2. \code{medAbsErr()} returns the median absolute error of prediction for a numeric
#' response \code{y} and predictions \code{yhat}.
#' 3. \code{BayesRule()} and \code{BayesRule2()} report the proportion
#' of incorrect predictions for a dichotomous response variable \code{y}, assumed
#' coded (or coercible to) \code{0} and \code{1}. The \code{yhat} values are
#' predicted probabilities and are rounded to 0 or 1. The distinction
#' between \code{BayesRule()} and \code{BayesRule2()} is that the former
#' checks that the \code{y} values are all either \code{0} or \code{1}
#' and that the \code{yhat} values are all between 0 and 1, while
#' the latter doesn't and is therefore faster.
#'
#' @returns In general, cost functions should return a single numeric
#' value measuring lack-of-fit. \code{mse()} returns the mean-squared error;
#' \code{rmse()} returns the root-mean-squared error;
#' \code{medAbsErr()} returns the median absolute error;
#' and \code{BayesRule()} and
#' \code{BayesRule2()} return the proportion of misclassified cases.
#' @examples
#' data("Duncan", package="carData")
#' m.lm <- lm(prestige ~ income + education, data=Duncan)
#' mse(Duncan$prestige, fitted(m.lm))
#'
#' data("Mroz", package="carData")
#' m.glm <- glm(lfp ~ ., data=Mroz, family=binomial)
#' BayesRule(Mroz$lfp == "yes", fitted(m.glm))

#' @describeIn cost-functions Mean-square error.
#' @export
mse <- function(y, yhat) {
  result <- mean((y - yhat) ^ 2)
  attr(result, "casewise loss") <- "(y - yhat)^2"
  result
}

#' @describeIn cost-functions Root-mean-square error.
#' @export
rmse <- function(y, yhat) {
  sqrt(mean((y - yhat) ^ 2))
}

#' @describeIn cost-functions Median absolute error.
#' @importFrom stats median
#' @export
medAbsErr <- function(y, yhat) {
  median(abs(y - yhat))
}

#' @describeIn cost-functions Bayes Rule for a binary response.
#' @export
BayesRule <- function(y, yhat) {
  if (!all(y %in% c(0, 1)))
    stop("response values not all 0 or 1")
  if (any(yhat < 0) ||
      any(yhat > 1))
    stop("fitted values outside of interval [0, 1]")
  yhat <- round(yhat)
  result <- mean(y != yhat) # proportion in error
  attr(result, "casewise loss") <- "y != round(yhat)"
  result
}

#' @describeIn cost-functions Bayes rule for a binary response (without bounds checking).
#' @export
BayesRule2 <- function(y, yhat) {
  yhat <- round(yhat)
  result <- mean(y != yhat) # proportion in error
  attr(result, "casewise loss") <- "y != round(yhat)"
  result
}
