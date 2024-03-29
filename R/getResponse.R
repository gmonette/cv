#' Extract Response Variable
#'
#' Generic function to extract the response variable from a fitted model.
#'
#' @param model a fitted model
#' @param ... additional parameters for specific methods
#'
#' @returns a numeric vector containing the values of the response variable.
#'
#' @details
#' The supplied \code{default} method returns the \code{model$y} component
#' of the model object, or, if \code{model} is an S4 object, the result
#' returned by the \code{\link[insight]{get_response}()} function in
#' the \pkg{insight} package. If this result is \code{NULL}, the result of
#' \code{model.response(model.frame(model))} is returned, checking in any case whether
#' the result is a numeric vector.
#'
#' There is also an \code{"lme"} method, and \code{"merMod"}
#' and \code{"glmmTMB"} methods that convert factor
#' responses to numeric 0/1 responses, as would be appropriate
#' for a generalized linear mixed model with a binary response.
#'
#' @examples
#'     fit <- lm(mpg ~ gear, mtcars)
#'     GetResponse(fit)
#' @export
GetResponse <- function(model, ...){
  UseMethod("GetResponse")
}

#' @describeIn GetResponse \code{default} method
#' @export
GetResponse.default <- function(model, ...){
  y <- if (!isS4(model)) model$y else insight::get_response(model)
  if (is.null(y)) y <- model.response(model.frame(model))
  if (!is.vector(y)) stop("non-vector response")
  if (!is.numeric(y)) stop("non-numeric response")
  y
}

#' @describeIn GetResponse \code{merMod} method
#' @export
GetResponse.merMod <- function(model, ...){
  y <- insight::get_response(model)
  if (is.factor(y)) {
    levels <- levels(y)
    failure <- levels[1]
    if (length(levels) > 2){
      message("Note: the response has more than 2 levels.\n",
              " The first level ('", failure,
              "') denotes failure (0),\n",
              " the others success (1)")
    }
    y <- as.numeric(y != failure)
  }
  if (!is.vector(y)) stop("non-vector response")
  if (!is.numeric(y)) stop("non-numeric response")
  y
}

#' @describeIn GetResponse \code{merMod} method
#' @export
GetResponse.lme <- function(model, ...) insight::get_response(model)

#' @describeIn GetResponse \code{glmmTMB} method
#' @export
GetResponse.glmmTMB <- GetResponse.merMod

#' @describeIn GetResponse \code{modList} method
#' @export
GetResponse.modList <- function(model, ...) GetResponse(model[[1]])
