#' Cross-Validate Mixed-Effects Model
#'
#' \code{\link{cv}()} methods for time-series models of class \code{"gls"}, fit
#' by \code{\link[nlme]{gls}()} in the \pkg{lme4} package, and
#' for models fit by \code{Arima()}, which provides a formula interface to the
#' \code{\link{arima}()}) function.
#'
#' @param model an object of class \code{"gls"} produced by the
#' \code{\link[nlme]{gls}()} function, or of class \code{"ARIMA"}
#' produced by the \code{Arima()} function.
#' @param data data frame with the data to which the model was fit;
#' can usually be inferred from the model; for \code{Arima()}, a data
#' frame with data to which the model is to be fit.
#' @param criterion function to compute the CV cost criterion
#' (default \code{\link{mse}}).
#' @param k number of folds, an integer \eqn{\gt 2}; the default is \code{10}.
#' @param reps ignored (to match \code{\link{cv}()} generic function).
#' @param seed ignored (to match \code{\link{cv}()} generic function).
#' @param fold.type if \code{"cumulate"} (the default), predict
#' the response for cases in the i-th fold from the model fit to data
#' all preceding folds; if \code{"preceding"}, predict using cases in
#' the immediately preceding fold only; if \code{"all"}, predict using all
#' other folds, including those in the future (available for \code{"gls"}
#' models but not for \code{"ARIMA"} models).
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
#' @examples
#'
#' # model from help("gls", package="nlme")
#' if (require("nlme", quietly=TRUE)){
#' withAutoprint({
#' fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
#'            correlation = corAR1(form = ~ 1 | Mare))
#' fm1
#' summary(cv(fm1, k=5, confint = TRUE))
#' summary(cv(fm1, k=5, confint = TRUE, fold.type="preceding"))
#' summary(cv(fm1, k=5, confint = TRUE, fold.type="all"))
#' })
#' }
#'
#' if (require("stats", quietly=TRUE) &&
#'     require("datasets", quietly=TRUE)){
#' withAutoprint({
#' # model adapted from help("arima")
#' LH <- data.frame(lh = lh)
#' lh.arima <- Arima(~lh, data=LH)
#' lh.arima
#' summary(cv(lh.arima, k=5))
#' summary(cv(lh.arima, k=5, fold.type="preceding"))
#'
#' # model adapted from help("arima")
#' Lake <- data.frame(level=LakeHuron, year=time(LakeHuron))
#' lake.arima <- Arima(level ~ I(year - 1920), data=Lake,
#'                   order=c(2, 0, 0))
#' lake.arima
#' summary(cv(lake.arima, k=5))
#' summary(cv(lake.arima, k=5, fold.type="preceding"))
#' })
#' }

#' @describeIn cv.gls \code{cv()} method for \code{\link[nlme]{gls}()} models from the \pkg{nlme} package.
#' @export
cv.gls <- function(model,
                   data = insight::get_data(model),
                   criterion = mse,
                   k = 10L,
                   reps,
                   seed,
                   criterion.name = deparse(substitute(criterion)),
                   fold.type=c("cumulative", "preceding", "all"),
                   details = k <= 10L,
                   confint = n >= 400L,
                   level = 0.95,
                   ncores = 1L,
                   ...) {

  f <- function(i, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i, fold.type = fold.type)
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
    indices.i <- fold(folds, i, fold.type = fold.type)
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

  fold.type <- match.arg(fold.type)

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
    fold.type = fold.type,
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
#' \code{\link{na.pass}}, will pass missing data to the \code{\link{arima}()}
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

#' @describeIn cv.gls \code{model.matrix()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
model.matrix.ARIMA <- function(object, ...) object$model.matrix

#' @param n.ahead number of future cases to predict.
#' @param newdata data frame with rows containing the
#' predictors (if any) for predicted future cases.
#' @param se.fit if \code{TRUE} (the default is \code{FALSE}), compute
#' the standard errors of the predictions.
#'
#' @returns the \code{\link{cv}()} methods return objects of class \code{"cv"}.
#'   \code{Arima()} returns an object of class \code{"Arima"} with the
#'   following components: \code{formula}, the model formula; \code{data},
#'   the data set to which the model was fit; \code{subset}, the subset
#'   expression (if specified); \code{na.action}, see \code{\link{na.pass}}; \code{order}, the
#'   order of the ARIMA model; \code{call}, the function call;
#'   \code{dots}, any other arguments specified; \code{arima},
#'   the object returned by the \code{\link{arima}()} function,
#'   which \code{Arima()} calls;
#'   \code{response}, the response variable; \code{model.matrix},
#'   the model matrix, if there are predictors in the model.
#'
#' @describeIn cv.gls \code{predict()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
predict.ARIMA <- function(object, n.ahead, newdata = NULL,
                          se.fit = FALSE, ...){
  if (missing(n.ahead) && is.null(newdata)) return(NULL)
  x <- model.matrix(object)
  new.x <- if (!is.null(newdata) && !is.null(model.matrix(object))) {
    n.ahead <- nrow(newdata)
    model.frame(object$formula[-2], data=newdata)
  } else {
    NULL
  }
  predict(object$arima, n.ahead=n.ahead, newxreg=new.x,
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
                   fold.type = c("cumulative", "preceding"),
                   criterion.name = deparse(substitute(criterion)),
                   details = k <= 10L,
                   ncores = 1L,
                   ...) {

  f <- function(i, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i, fold.type=fold.type)
    indices.j <- fold(folds, i, predict=TRUE)
    model.i <- update(model, data = data[indices.i, , drop=FALSE])
    fit.i <- predict(model.i, n.ahead=length(indices.j),
                     newdata=data[indices.j, , drop=FALSE])
    list(
      fit.i = fit.i,
      coef.i = coef(model.i)
    )
  }

  fPara <- function(i,
                    model.function = model.function,
                    model.function.name,
                    ...) {
    # helper function to compute cv criterion for each fold
    #  with parallel computations
    indices.i <- fold(folds, i, fold.type=fold.type)
    indices.j <- fold(folds, i, predict=TRUE)
    assign(model.function.name, model.function)
    # the following deals with a scoping issue that can
    #   occur with args passed via ... (which is saved in dots)
    predict.args <- c(list(
      object = update(model, data = data[indices.i, , drop=FALSE]),
      newdata = data[indices.j, , drop=FALSE],
      n.ahead = length(indices.j)
    ), dots)

    fit.i <- do.call(predict, predict.args)
    list(
      fit.i = fit.i,
      coef.i = coef(predict.args$object)
    )
  }

  fold.type <- match.arg(fold.type)

  x.names <- colnames(model$model.matrix)

  call <- getCall(model)
  model.function <- call[[1L]]
  model.function.name <- as.character(model.function)
  model.function <- eval(model.function)

  result <- cvOrdered(
    model = model,
    data = data,
    criterion = criterion,
    criterion.name = criterion.name,
    k = k,
    fold.type=fold.type,
    confint=FALSE,
    details = details,
    ncores = ncores,
    f = f,
    fPara = fPara,
    model.function = model.function,
    model.function.name = model.function.name,
    locals = list(x.names = x.names),
    ...
  )
  result[["full crit"]] <- NULL
  result
}
