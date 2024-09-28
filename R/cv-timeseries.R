#' Cross-Validate Mixed-Effects Model
#'
#' \code{\link{cv}()} methods for time-series models of class \code{"gls"}, fit
#' by the \code{\link[nlme]{gls}()} in the \pkg{lme4} package, and
#' by \code{Arima()}, which provides a formula interface to the
#' \code{link{arima}()}) function.
#'
#' @param model an object of class \code{"gls"} produced by the
#' \code{\link[nlme]{gls}()} function, or of class \code{"ARIMA"}
#' produced by the \code{Arima()} function.
#' @param data data frame with the data to which the model was fit;
#' can usually be inferred from the model; for \code{Arima()}, a data
#' frame with data to which the model is to be fit.
#' @param criterion function to compute the CV cost criterion
#' (default \code{\link{mse}}).
#' @param k number of folds, an integer or "n" or "loo" for
#' leave-one-out CV; the default is \code{10}.
#' @param reps ignored (to match \code{\link{cv}()} generic function.)
#' @param seed ignored (to match \code{\link{cv}()} generic function.)
#' @param i.only if \code{TRUE} (the default is \code{FALSE}), predict
#' the response for cases in the i-th fold from the model fit to data
#' in the preceding fold only rather than fit to data from \emph{all}
#' preceding folds.
#' @param criterion.name name of the CV criterion; can usually be
#' inferred from \code{criterion}.
#' @param details return fold-wise statistics, including the CV criterion
#' and parameter estimates for cases in each fold after the first;
#' the default is \code{TRUE} for \eqn{k \le 10}.
#' @param confint if \code{TRUE} (the default if \eqn{n \ge 400}), report
#' a confidence interval for the adjusted CV criterion.
#' @param level level for the confidence interval (default \code{0.95}).
#' @param ncores if \code{ncores} \eqn{> 1}, the computation is parallelized.
#'
#' @describeIn cv.gls \code{cv()} method for \code{\link[nlme]{gls}()} models from the \pkg{nlme} package.
#' @export
cv.gls <- function(model,
                   data = insight::get_data(model),
                   criterion = mse,
                   k = 10L,
                   reps,
                   seed,
                   criterion.name = deparse(substitute(criterion)),
                   i.only = FALSE,
                   details = k <= 10L,
                   confint = n >= 400L,
                   level = 0.95,
                   ncores = 1L,
                   ...) {

  f <- function(i, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i, i.only = i.only)
    model.i <- update(model, data = data[indices.i, ])
    fit.all.i <- predict(model.i, newdata = data, ...)
    fit.i <- fit.all.i[fold(folds, i, predict = TRUE)]
    list(
      fit.i = fit.i,
      crit.all.i = criterion(y, fit.all.i),
      coef.i = coef(model.i)
    )
  }

  fPara <- function(i,
                    model.function = model.function,
                    model.function.name,
                    cor.function = cor.function,
                    cor.function.name = cor.function.name,
                    ...) {
    # helper function to compute cv criterion for each fold
    #  with parallel computations
    indices.i <- fold(folds, i, i.only = i.only)
    assign(model.function.name, model.function)
    assign(cor.function.name, cor.function)
    # the following deals with a scoping issue that can
    #   occur with args passed via ... (which is saved in dots)
    predict.args <- c(list(
      object = update(model, data = data[indices.i, ]),
      newdata = data
    ), dots)
    fit.all.i <- do.call(predict, predict.args)
    fit.i <- fit.all.i[fold(folds, i, predict = TRUE)]
    list(
      fit.i = fit.i,
      crit.all.i = criterion(y, fit.all.i),
      coef.i = coef(predict.args$object)
    )
  }

  call <- getCall(model)

  model.function <- call[[1L]]
  model.function.name <- as.character(model.function)
  model.function <- eval(model.function)
  cor.function <- getCall(model)[[which(names(call) == "correlation")]][[1]]
  cor.function.name <- as.character(cor.function)
  cor.function <- eval(cor.function)

  n <- nrow(data)

  cvOrdered(
    model = model,
    data = data,
    criterion = criterion,
    criterion.name = criterion.name,
    k = k,
    i.only = i.only,
    details = details,
    confint = confint,
    level = level,
    ncores = ncores,
    f = f,
    fPara = fPara,
    model.function = model.function,
    model.function.name = model.function.name,
    cor.function = cor.function,
    cor.function.name = cor.function.name,
    ...
  )
}

#' @param formula either a one-sided formula giving the response variable
#' in an ARIMA model with no predictors, or a two-sided formula with the
#' response on the left and terms for the predictors on the right (as in
#' a typical R regression model).
#' @param subset subsetting expression.
#' @param na.action function to process missing data; the default,
#' \code{na.pass}, will pass missing data to the \code{\link{arima}()}
#' function.
#' @param order the \eqn{p, d, q} specification of the ARIMA model;
#' see \code{\link{arima}()} for details.
#' @param ... further arguments to be passed to \code{\link{arima}()}
#' or \code{Arima()}.
#'
#' @describeIn cv.gls model-formula wrapper for the \code{\link{arima}()} function.
#' @export
Arima <- function(formula, data, subset=NULL, na.action=na.pass,
                  order = c(1L, 0L, 0L), ...){
  dots <- list(...)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)
  which.int <- which("(Intercept)" == colnames(x))
  if (length(which.int > 0)) x <- x[, -which.int, drop=FALSE]
  result <- list(formula=formula, data=data, subset=subset,
                 na.action=na.action, order=order, call=cl, dots=dots)
  if (length(formula) == 2){
    if (!(ncol(x) == 1) && is.numeric(x[, 1]))
      stop("formula must specify a single response")
    result$arima <- stats::arima(x[, 1], order=order, ...)
    result$response <- x[, 1]
  } else {
    y <- model.response(mf, "numeric")
    if (!is.vector(y) && is.numeric(y))
      stop("formula must specify a single response")
    result$arima <- stats::arima(y, order=order, xreg=x, ...)
    result$response <- y
    result$model.matrix <- x
  }
  class(result) <- "ARIMA"
  result
}

#' @param x an object of class \code{"ARIMA"}.
#' @describeIn cv.gls \code{print()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
print.ARIMA <- function(x, ...){
  xx <- x$arima
  xx$call <- getCall(x)
  print(xx)
  invisible(x)
}

#' @param object an object of class \code{"ARIMA"}.
#' @describeIn cv.gls \code{update()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
#'
update.ARIMA <- function(object, ...){
  cl0 <- object$call
  cl <- match.call()
  args <- object$dots
  args[c("formula", "data", "subset", "na.action", "order")] <-
    object[c("formula", "data", "subset", "na.action", "order")]
  dots <- list(...)
  names <- names(dots)
  for (name in names){
    args[[name]]  <- dots[[name]]
  }
  result <- do.call(Arima, args)
  names <- names(cl)
  for (name in names){
    if (name == "" || name == "object") next
    cl0[[name]] <- cl[[name]]
  }
  result$call <- cl0
  result
}

#' @describeIn cv.gls \code{coef()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
coef.ARIMA <- function(object, ...) coef(object$arima)

#' @param n.ahead number of future cases to predict; the default is \code{1}.
#' @param newdata data frame with \code{n.ahead} rows containing the
#' predictors (if any) for the predicted future cases.
#' @param se.fit if \code{TRUE} (the default is \code{FALSE}), compute
#' the standard errors of the predictions.
#' @describeIn cv.gls \code{predict()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
predict.ARIMA <- function(object, n.ahead = 1L, newdata = NULL,
                          se.fit = FALSE, ...){
  if (n.ahead == 0) {
    return(rep(NA, length(object$response)))
  }
  x <- object$data[ , colnames(newdata), drop=FALSE]
  predict(object$arima, n.ahead=n.ahead, newxreg=newdata,
          se.fit=se.fit, ...)

}

#' @describeIn cv.gls \code{cv()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
cv.ARIMA <- function(model,
                   data = model$data,
                   criterion = mse,
                   k = 10L,
                   reps,
                   seed,
                   i.only = FALSE,
                   criterion.name = deparse(substitute(criterion)),
                   details = k <= 10L,
                   ...) {

  f <- function(i, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i, i.only = i.only)
    indices.j <- fold(folds, i, predict=TRUE)
    model.i <- update(model, data = data[indices.i, , drop=FALSE])
    fit.i <- predict(model.i, n.ahead=length(indices.j),
                     newdata=data[indices.j, x.names, drop=FALSE])
    list(
      fit.i = fit.i,
      coef.i = coef(model.i)
    )
  }

  x.names <- colnames(model$model.matrix)

  cvOrdered(
    model = model,
    data = data,
    criterion = criterion,
    criterion.name = criterion.name,
    k = k,
    i.only = i.only,
    details = details,
    f = f,
    locals = list(x.names = x.names),
    n.ahead=0, # passed to predict() for the full sample
               # to produce na for CV criterion
    ...
  )
}
