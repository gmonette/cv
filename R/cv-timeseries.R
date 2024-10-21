#' Cross-Validate Mixed-Effects Model
#'
#' \code{\link{cv}()} methods for time-series models
#' fit by \code{Arima()}, which provides a formula interface to the
#' \code{\link{arima}()}) function.
#'
#' @param model an object of class \code{"ARIMA"}
#' produced by the \code{Arima()} function.
#' @param data data frame with the data to which the model was fit;
#' can usually be inferred from the model; for \code{Arima()}, a data
#' frame with data to which the model is to be fit.
#' @param criterion function to compute the CV cost criterion
#' (default \code{\link{mse}}).
#' @param k number of folds, an integer \eqn{\gt 2}; the default is \code{"n"}, in which case
#' the first fold is determined by \code{begin.with} and subsequent folds
#' each contain a single case; that makes sense if \code{fold.type = "cumulative"}.
#' @param reps ignored (to match \code{\link{cv}()} generic function).
#' @param seed ignored (to match \code{\link{cv}()} generic function).
#' @param fold.type if \code{"cumulative"} (the default), predict
#' the response for 1 or more cases after the i-th fold from the model fit to data
#' in the \code{i}th and all preceding folds; if \code{"preceding"}, predict using cases in
#' the \code{i}th fold only; if \code{"window"}, folds comprise a moving window
#' of \code{begin.with} cases.
#' @param lead how far ahead to predict (can be a vector of positive integers);
#' the default is \code{1}.
#' @param criterion.name name of the CV criterion; can usually be
#' inferred from \code{criterion}.
#' @param details return fold-wise parameter estimates for cases in each fold after the first;
#' the default is \code{TRUE} for \eqn{k \le 10,000}.
#' @param ncores if \code{ncores} \eqn{> 1}, the computation is parallelized.
#'
#' @examples
#' if (require("stats", quietly=TRUE) &&
#'     require("datasets", quietly=TRUE)){
#' withAutoprint({
#' # model adapted from help("arima")
#' LH <- data.frame(lh = lh)
#' lh.arima <- Arima(~lh, data=LH)
#' lh.arima
#' plot(lh.arima)
#' summary(cv.lh <- cv(lh.arima, lead=1:5))
#' plot(cv.lh)
#' summary(cv(lh.arima, lead=1:5, fold.type="window"))
#' # too few folds (5), just to illustrate fold.type="preceding":
#' summary(cv(lh.arima, k=5, fold.type="preceding"))
#'
#' # model adapted from help("arima")
#' Lake <- data.frame(level=LakeHuron, year=time(LakeHuron))
#' lake.arima <- Arima(level ~ I(year - 1920), data=Lake,
#'                   order=c(2, 0, 0))
#' lake.arima
#' plot(lake.arima)
#' summary(cv.lake <- cv(lake.arima, lead=1:5))
#' plot(cv.lake)
#' })
#' }

#' @param formula either a one-sided formula giving the response variable
#' in an ARIMA model with no predictors, or a two-sided formula with the
#' response on the left and terms for the predictors on the right (as in
#' a typical R regression model); if the timeseries is differenced, the intercept
#' is removed from the model even if the formula implies an intercept.
#' @param subset subsetting expression.
#' @param na.action function to process missing data; the default,
#' \code{\link{na.pass}}, will pass missing data to the \code{\link{arima}()}
#' function.
#' @param order the \eqn{p, d, q} specification of the ARIMA model;
#' see \code{\link{arima}()} for details; the default is \eqn{p = 1, d = 0, q = 0},
#' an AR(1) model.
#' @param seasonal specification of the seasonal part of the ARIMA model;
#' see \code{\link{arima}()} for details; the default is not to include
#' a seasonal part of the model.
#' @param ... further arguments to be passed to \code{\link{arima}()}
#' or \code{Arima()}.
#'
#' @describeIn Arima model-formula wrapper for the \code{\link{arima}()} function.
#' @export
Arima <- function(formula, data, subset=NULL, na.action=na.pass,
                  order = c(1L, 0L, 0L),
                  seasonal = list(order = c(0L, 0L, 0L), period = NA),
                  ...){
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
                 na.action=na.action, order=order, seasonal=seasonal,
                 call=cl, dots=dots)
  if (length(formula) == 2){
    if (!(ncol(x) == 1) && is.numeric(x[, 1]))
      stop("formula must specify a single response")
    result$arima <- stats::arima(x[, 1], order=order, seasonal=seasonal, ...)
    result$response <- x[, 1]
  } else {
    y <- model.response(mf, "numeric")
    if (!is.vector(y) && is.numeric(y))
      stop("formula must specify a single response")
    result$arima <- stats::arima(y, order=order, seasonal=seasonal, xreg=x, ...)
    result$response <- y
    result$model.matrix <- x
  }
  class(result) <- "ARIMA"
  result
}

#' @param x an object of class \code{"ARIMA"}.
#' @describeIn Arima \code{print()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
print.ARIMA <- function(x, ...){
  xx <- x$arima
  xx$call <- getCall(x)
  print(xx)
  invisible(x)
}

#' @param y ignored, to match \code{plot()} generic.
#' @param xlab label for horizontal ("time") axis; defaults to
#' \code{"Time"}.
#' @param main title for diagnostic plots.
#' @param col color for points and lines.
#' @describeIn Arima \code{plot()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
plot.ARIMA <- function(x, y, xlab="time",
                       main="Diagnosic Plots",
                       col=car::carPalette()[2], ...){
  residuals <- residuals(x)
  save <- par(mfrow=c(2, 2), oma=c(0, 0, 1, 0))
  on.exit(par(save))
  plot(fitted(x), xlab=xlab, ylab=expression(hat(y)),
       main = "Fitted Values", type="b", pch=16, col=col)
  grid(lty=2, col="gray")
  plot(residuals, xlab=xlab, ylab="residuals",
       main="Residuals", type="b", pch=16, col=col)
  grid(lty=2, col="gray")
  acf(residuals, main="Autocorrelations\n of Residuals",
      na.action=na.pass)
  pacf(residuals, main="Parial Autocorrelations\n of Residuals",
       na.action=na.pass)
  title(main=main, outer=TRUE)
}

#' @param object an object of class \code{"ARIMA"}.
#' @describeIn Arima \code{update()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
#'
update.ARIMA <- function(object, ...){
  cl0 <- object$call
  cl <- match.call()
  args <- object$dots
  args[c("formula", "data", "subset", "na.action", "order", "seasonal")] <-
    object[c("formula", "data", "subset", "na.action", "order", "seasonal")]
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

#' @describeIn Arima \code{coef()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
coef.ARIMA <- function(object, ...) coef(object$arima)

#' @describeIn Arima \code{model.matrix()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
model.matrix.ARIMA <- function(object, ...) object$model.matrix

#' @param n.ahead number of future cases to predict.
#' @param newdata data frame with rows containing the
#' predictors (if any) for predicted future cases.
#' @param se.fit if \code{TRUE} (the default is \code{FALSE}), compute
#' the standard errors of the predictions.
#'
#' @returns the \code{\link{cv.ARIMA}()} method returns an object of class
#'   \code{c("cvOrdered", "cv")} (see  \code{\link{cvOrdered}()} and \code{\link{cv}()}).
#'
#'   \code{Arima()} returns an object of class \code{"ARIMA"} with the
#'   following components: \code{formula}, the model formula; \code{data},
#'   the data set to which the model was fit; \code{subset}, the subset
#'   expression (if specified); \code{na.action}, see \code{\link{na.pass}};
#'   \code{order}, the order of the ARIMA model;
#'   \code{seasonal}, the seasonal specification;
#'   \code{call}, the function call;
#'   \code{dots}, any other arguments specified; \code{arima},
#'   the object returned by the \code{\link{arima}()} function,
#'   which \code{Arima()} calls;
#'   \code{response}, the response variable; \code{model.matrix},
#'   the model matrix, if there are predictors in the model.
#'
#' @describeIn Arima \code{fitted()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
fitted.ARIMA <- function(object, ...){
  residuals <- object$arima$residuals
  y <- object$response
  y - residuals
}
#' @describeIn Arima \code{residuals()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
residuals.ARIMA <- function(object, ...) object$arima$residuals
#'
#' @describeIn Arima \code{predict()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
predict.ARIMA <- function(object, n.ahead, newdata = NULL,
                          se.fit = FALSE, ...){
  if (missing(n.ahead) && is.null(newdata)){
    return(fitted(object))
  }
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

#' @param begin.with the number of cases in
#' the first fold. The remaining cases are divided among the subsequent
#' \code{k} - 1 folds.
#' @describeIn Arima \code{cv()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
cv.ARIMA <- function(model,
                   data = model$data,
                   criterion = mse,
                   k = "n",
                   reps,
                   seed,
                   fold.type = c("cumulative", "preceding", "window"),
                   begin.with=max(25, ceiling(n/10)),
                   lead = 1L,
                   criterion.name = deparse(substitute(criterion)),
                   details = n <= 1e4,
                   ncores = 1L,
                   ...) {

  f <- function(i, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i) #, fold.type=fold.type)
    indices.j <- fold(folds, i, predict=TRUE, lead=lead)
    model.i <- update(model, data = data[indices.i, , drop=FALSE])
    fit.i <- predict(model.i, n.ahead=max(lead),
                     newdata=data[indices.j, , drop=FALSE])
    list(
      fit.i = fit.i[lead],
      coef.i = coef(model.i)
    )
  }

  fPara <- function(i,
                    model.function,
                    model.function.name,
                    ...) {
    # helper function to compute cv criterion for each fold
    #  with parallel computations
    indices.i <- fold(folds, i)
    indices.j <- fold(folds, i, predict=TRUE, lead=lead)
    # bring Arima() into scope
    assign(model.function.name, model.function)
    # the following deals with a scoping issue that can
    #   occur with args passed via ... (which is saved in dots)
    predict.args <- c(list(
      object = update(model, data = data[indices.i, , drop=FALSE]),
      newdata = data[indices.j, , drop=FALSE],
      n.ahead = max(lead)), dots)

    fit.i <- do.call(predict, predict.args)
    list(
      fit.i = fit.i[lead],
      coef.i = coef(predict.args$object)
    )
  }

  n <- nrow(data)

  fold.type <- match.arg(fold.type)

  x.names <- colnames(model$model.matrix)

  model.function.name <- "Arima"
  model.function <- Arima

  result <- cvOrdered(
    model = model,
    data = data,
    criterion = criterion,
    criterion.name = criterion.name,
    k = k,
    fold.type=fold.type,
    begin.with = begin.with,
    lead = lead,
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
  result[["full crit"]] <- criterion(na.omit(GetResponse(model)),
                                     na.omit(fitted(model)))
  result[["lead"]] <- lead
  result
}

