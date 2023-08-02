#' Extract Response Variable(s)
#'
#' Generic function to extract the response variable(s) from a fitted model.
#'
#' @param model a fitted model
#' @param ... additional parameters for specific methods
#'
#' @returns a vector or matrix containing the values of response variable(s)
#'
#' @examples
#'     fit <- lm(cbind(hp, mpg) ~ gear, mtcars)
#'     getResponse(fit)
#' @export
getResponse <- function(model, ...){
  UseMethod("getResponse")
}

#' @describeIn getResponse default method
#' @export
getResponse.default <- function(model, ...){
  y <- model$y
  if (is.null(y)) y <- model.response(model.frame(model))
  if (!is.numeric(y)) stop("non-numeric response")
  y
}

