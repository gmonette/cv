#' Extract Response Variable(s)
#'
#' Generic function to extract the response variable(s) from a fitted model.
#'
#' @param model a fitted model
#' @param ... additional parameters for specific methods
#'
#' @returns a vector or matrix containing the values of response variable(s)
#'
#' @details
#' The supplied \code{default} method returns the \code{model$y} component
#' of the model object if it exists and otherwise the result of
#' \code{model.response(model.frame(model))}, checking in either case whether
#' the result is numeric.
#'
#' @examples
#'     fit <- lm(cbind(hp, mpg) ~ gear, mtcars)
#'     getResponse(fit)
#' @export
getResponse <- function(model, ...){
  UseMethod("getResponse")
}

#' @describeIn getResponse \code{default} method
#' @export
getResponse.default <- function(model, ...){
  y <- model$y
  if (is.null(y)) y <- model.response(model.frame(model))
  if (!is.numeric(y)) stop("non-numeric response")
  y
}

