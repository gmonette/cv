#' Cross-Validate Several Models Fit to the Same Data
#'
#' A \code{\link{cv}()} method for an object of class  \code{"modlist"},
#' created by the \code{models()} function. This \code{cv()} method simplifies
#' the process of cross-validating several models on the same set of folds.
#' \code{models()} performs some
#' "sanity" checks, warning if the models are of different classes, and
#' reporting an error if they are fit to apparently different data sets or
#' different response variables.
#' @param model a list of regression model objects,
#' created by \code{models()}.
#' @param data (required) the data set to which the models were fit.
#' @param criterion the CV criterion (cost) function, defaults to
#' \code{\link{mse}}.
#' @param k the number of CV folds; may be omitted, in which case the value
#' will depend on the default for the \code{cv()} method invoked for the
#' individual models.
#' @param reps ignored (present only to match the generic function).
#' @param seed (optional) seed for R's pseudo-random-number generator,
#' to be used to create the same set of CV folds for all of the models;
#' if omitted, a seed will be randomly generated and saved.
#' @param quietly If \code{TRUE} (the default), simple messages (for example about the
#' value to which the random-number generator seed is set), but not warnings or
#' errors, are suppressed.
#' @param ... for \code{cv()}, additional arguments to be passed to the \code{cv()} method
#' applied to each model. For \code{models()}, two or more competing models fit to the
#' the same data; the several models may be named. For the \code{print()}
#' method, arguments to be passed to the \code{print()} method for
#' the individual model cross-validations.
#' @param x an object of class \code{"cvModList"} to be printed.
#' @return \code{models()} returns a \code{"modList"} object, the
#' \code{cv()} method for which returns a \code{"cvModList"} object.
#' @examples
#' data("Duncan", package="carData")
#' m1 <- lm(prestige ~ income + education, data=Duncan)
#' m2 <- lm(prestige ~ income + education + type, data=Duncan)
#' m3 <- lm(prestige ~ (income + education)*type, data=Duncan)
#' cv(models(m1=m1, m2=m2, m3=m3), data=Duncan, seed=7949)

#' @describeIn models create a list of models
#' @export
models <- function(...){
  models <- list(...)
  if (length(models) < 2L) stop("fewer than 2 models to be compared")
  classes <- sapply(models, function(m) class(m)[1L])
  n <- sapply(models, function(m) nrow(insight::get_data(m)))
  if (!all(n[1L] == n[-1L])) {
    stop("models are fit to data sets of differing numbers of cases")
  }
  response <- getResponse(models[[1L]])
  for (i in 2L:length(models)){
    if (!isTRUE(all.equal(response, getResponse(models[[i]]),
                          check.attributes=FALSE))){
      stop("models are not all fit to the same response variable")
    }
  }
  if (length(unique(classes)) > 1L)
    warning("models are not all of the same primary class")
  nms <- names(models)
  if (is.null(nms)) {
    names(models) <- paste0("model.", seq_along(models))
  } else {
    unnamed <- which(nms == "")
    names(models)[unnamed] <- paste0("model.", seq_along(unnamed))
  }
  class(models) <- "modList"
  models
}

#' @describeIn models \code{cv()} method for \code{"modList"} objects
#' @exportS3Method
cv.modList <- function(model, data, criterion=mse, k, reps, seed, quietly=TRUE, ...){
  n.models <- length(model)
  if (missing(seed)) seed <- sample(1e6, 1L)
  result <- vector(n.models, mode="list")
  names(result) <- names(model)
  class(result) <- "cvModList"
  for (i in 1L:n.models){
    result[[i]] <- if (missing(k)){
      if (quietly){
        suppressMessages(cv(model[[i]], data=data, seed=seed, ...))
      } else {
        cv(model[[i]], data=data, seed=seed, ...)
      }
    } else {
      if (quietly){
        suppressMessages(cv(model[[i]], data=data, k=k, seed=seed, ...))
      } else {
        cv(model[[i]], data=data, k=k, seed=seed, ...)
      }
    }
  }
  result
}

#' @describeIn models \code{print()} method for \code{"cvModList"} objects
#' @exportS3Method
print.cvModList <- function(x, ...){
  nms <- names(x)
  for (i in seq_along(x)){
    cat(paste0("\nModel ", nms[i], ":\n"))
    print(x[[i]], ...)
  }
  return(invisible(x))
}

#' @export
`[.cvModList` <- function(x, ...){
  result <- NextMethod()
  class(result) <- "cvModList"
  result
}
