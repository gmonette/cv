#' Fit and Cross-Validate ARIMA Models
#'
#' \code{\link{cv}()} methods for time-series models
#' fit by \code{Arima()}, which provides a formula interface to the
#' \code{\link{arima}()}) function.
#'
#' @param mod an object of class \code{"ARIMA"}
#' produced by the \code{Arima()} function.

#' @examples
#' if (require("stats", quietly=TRUE) &&
#'     require("datasets", quietly=TRUE) &&
#'     require(splines) &&
#'     require("car") &&
#'     require("effects")){
#' withAutoprint({
#' # model adapted from help("arima")
#' LH <- data.frame(lh = lh)
#' lh.arima <- Arima(~lh, data=LH)
#' summary(lh.arima)
#' plot(lh.arima)
#' testArima(lh.arima)
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
#' summary(lake.arima)
#' Anova(lake.arima, type=3)
#' linearHypothesis(lake.arima, hypothesis=c("ar1", "ar2"))
#' plot(lake.arima)
#' testArima(lake.arima)
#' summary(cv.lake <- cv(lake.arima, lead=1:5))
#' plot(cv.lake)
#' plot(Effect("year", lake.arima, residuals=TRUE))
#' lake.arima.bs <- update(lake.arima, . ~ ns(year, 3))
#' summary(cv(lake.arima.bs, lead=1:5))
#' plot(Effect("year", lake.arima.bs, residuals=TRUE))
#' lake.arima.quad <- update(lake.arima, . ~ poly(year, 2))
#' summary(cv(lake.arima.quad, lead=1:5, min.ahead=3))
#' plot(Effect("year", lake.arima.quad, residuals=TRUE))
#' plot(cv(models(linear = lake.arima,
#'                bspline = lake.arima.bs,
#'                quadratic = lake.arima.quad),
#'         lead=1:5, min.ahead=3, data=Lake))
#' })
#' }

#' @param formula either a one-sided formula giving the response variable
#' in an ARIMA model with no predictors, or a two-sided formula with the
#' response on the left and terms for the predictors on the right (as in
#' a typical R regression model); if the timeseries is differenced, the intercept
#' is removed from the model even if the formula implies an intercept.
#' For the \code{\link[stats]{model.frame}()}
#' method (to match the generic), an \code{"ARIMA"} object.
#' @param order the \eqn{p, d, q} specification of the ARIMA model;
#' see \code{\link{arima}()} for details; the default is \eqn{p = 1, d = 0, q = 0},
#' an AR(1) model.
#' @param seasonal specification of the seasonal part of the ARIMA model;
#' see \code{\link{arima}()} for details; the default is not to include
#' a seasonal part of the model.
#' @param ... further arguments to be passed to \code{\link{arima}()},
#' \code{Arima()}, or other functions and methods.
#'
#' @describeIn Arima model-formula wrapper for the \code{\link{arima}()} function.
#' @export
Arima <- function(formula, data,
                  order = c(1L, 0L, 0L),
                  seasonal = list(order = c(0L, 0L, 0L), period = NA),
                  ...){
  if (!inherits(data, "ts_data_frame")) {
    data <- as.ts(data)
    message("Note: 'data' coerced to 'ts_data_frame'")
  }
  tsp <- tsp(data[, 1])
  dots <- list(...)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)
  which.int <- which("(Intercept)" == colnames(x))
  if (length(which.int > 0)) {
    has.intercept <- TRUE
    x <- x[, -which.int, drop=FALSE]
  } else {
    has.intercept <- FALSE
  }
  result <- list(formula=formula, data=data,
                 order=order, seasonal=seasonal,
                 call=cl, model=mf, dots=dots)
  if (length(formula) == 2){
    if (!(ncol(x) == 1) && is.numeric(x[, 1]))
      stop("formula must specify a single response")
    response <- x[, 1]
    if (!is.ts(response)) {
      tsp(response) <- tsp
      class(response) <- c("tsp", class(response))
    }
    result$arima <- stats::arima(response, order=order, seasonal=seasonal,
                                 ...)
    result$response <- response
  } else {
    y <- model.response(mf, "numeric")
    if (!(is.vector(y) || is.ts(y)) && is.numeric(y))
      stop("formula must specify a single response")
    if (!is.ts(y)) {
      tsp(y) <- tsp
      class(y) <- c("tsp", class(y))
    }
    result$arima <- stats::arima(y, order=order, seasonal=seasonal, xreg=x,
                                 include.mean=has.intercept, ...)
    result$response <- y
    result$model.matrix <- x
  }
  class(result) <- "ARIMA"
  result
}

#' @param x an object of class \code{"ARIMA"}.
#' @param digits for printing.
#' @describeIn Arima \code{print()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
print.ARIMA <- function (x, digits = max(3L, getOption("digits") - 3L),
                         ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Estimates:\n")
  print(format(coef(x), digits = digits), quote=FALSE)
  cat("\n")
  invisible(x)
}

#' @param correlation compute and report correlations of coefficients.
#' @describeIn Arima \code{summary()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
summary.ARIMA <- function (object, correlation=FALSE, ...) {
  arima <- object$arima
  call <- object$call
  coef <- coef(object)
  vcov <- vcov(object)
  corr <- if (correlation) cov2cor(vcov)
  sigma2 <- arima$sigma2
  logLik <- logLik(object)
  aic <- AIC(object)
  residuals <- residuals(object)
  result <- list(call=call, coef=coef, vcov=vcov, corr=corr,
                 sigma2=sigma2, logLik=logLik, aic=aic,
                 residuals = if (length(residuals) >= 5) residuals)
  class(result) <- "summary.ARIMA"
  result
}

#' @param signif.stars show "significance stars" in coefficient table?
#' @describeIn Arima \code{print()} method for \code{"summary.ARIMA"} objects.
#' @export
print.summary.ARIMA <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  residuals <- x$residuals
  if (!is.null(residuals)){
    quants <- quantile(residuals, c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE)
    names(quants) <- c("Min", "1st Q", "Median", "3rd Q", "Max")
    cat("Residuals:\n")
    print(quants, digits=digits)
    cat("\n")
  }
  coefMat <- cbind(b <- x$coef, se.b <- sqrt(diag(x$vcov)),
                   z <- b/se.b, p = 2*pnorm(abs(z), lower.tail=FALSE))
  rownames(coefMat) <- names(b)
  colnames(coefMat) <- c("Estimate", "Std. Error", "z value",
                         "Pr(>|z|")
  cat("Estimates:\n")
  printCoefmat(coefMat, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  corr <- x$corr
  if (!is.null(corr)){
    cat("\nCorrelations of coefficients:\n")
    print(corr, digits=digits)
  }
  cat("\nResidual standard deviation:", signif(sqrt(x$sigma2), digits))
  cat("\nLog-likelhood =", signif(x$logLik, digits),
      "\nAIC =", signif(x$aic, digits), "\n")
  invisible(x)
}

#' @describeIn Arima \code{plot()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function; calls \code{\link[stats]{tsdiag}()}.
#' @export
plot.ARIMA <- function(x, ...){
  tsdiag(x$arima, ...)
}

#' @describeIn Arima \code{\link[car]{Anova}()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @method Anova ARIMA
#' @exportS3Method car::Anova ARIMA
Anova.ARIMA <- function(mod, type = c("II", "III", 2, 3), ...){
  vc <- vcov(mod)
  coefs <- coef(mod)
  all.names <- names(coefs)
  b.names <- colnames(model.matrix(mod))
  has.intercept <- "(Intercept)" %in% all.names
  if (has.intercept) b.names <- c("(Intercept)", b.names)
  if (has.intercept) {
    model.matrix <- cbind(1, mod$model.matrix)
    colnames(model.matrix)[1] <- "(Intercept)"
    mod$model.matrix <- model.matrix
  }
  mod$terms <- terms(mod$formula)
  mod$arima$coef <- coefs[b.names]
  mod$arima$var.coef <- vc[b.names, b.names]
  NextMethod()
}

#' @param hypothesis.matrix specification of the linear hypothesis;
#' for details, see \code{\link[car]{linearHypothesis}()}.
#' @param rhs optional right-hand-side vector for the linear hypothesis;
#' for details, see \code{\link[car]{linearHypothesis}()}.
#' @describeIn Arima \code{\link[car]{linearHypothesis}()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @method linearHypothesis ARIMA
#' @exportS3Method car::linearHypothesis ARIMA
linearHypothesis.ARIMA <- function(model, hypothesis.matrix, rhs=NULL, ...){
  NextMethod(vcov. = vcov(model), coef. = coef(model),
             suppress.vcov.msg = TRUE)
}

#' @exportS3Method effects::effSources ARIMA
effSources.ARIMA <- function(mod){
  coefs <- coef(mod)
  vc <- vcov(mod)
  all.names <- names(coefs)
  b.names <- colnames(model.matrix(mod))
  has.intercept <- "(Intercept)" %in% all.names
  if (has.intercept) b.names <- c("(Intercept)", b.names)
  args <- list(
    coefficients = coefs[b.names],
    vcov = vc[b.names, b.names])
  args
}

#' @rdname Arima
#' @export
testArima <- function(model, ...){
  UseMethod("testArima")
}

#' @param lag maximum lag to compute residual autocorrelation; if not
#'   specified, the same default maximum lag as \code{\link[stats]{acf}}
#'   is used.
#' @param type test of autocorrelations, either \code{"Box-Pierce"} (the default)
#'   or \code{"Ljung-Box"} (see \code{\link[stats]{Box.test}()});
#'   for the \code{Anova()}, method, "type" of test (see \code{\link[car]{Anova}()}
#'   for details.)
#' @describeIn Arima test autocorrelations of ARIMA model residuals;
#'   the test is performed by \code{\link[stats]{Box.test}()}.
#' @importFrom stats Box.test logLik vcov AIC cov2cor pnorm quantile terms is.ts
#' KalmanForecast deltat ts tsp as.ts time tsdiag tsp<- window
#' @importFrom grDevices n2mfrow
#' @export
testArima.ARIMA <- function(model, lag = floor(10*log10(n)),
                            type = c("Box-Pierce", "Ljung-Box"),
                            ...){
  type <- match.arg(type)
  residuals <- residuals(model)
  n <- length(residuals)
  fitdf <- length(coef(model))
  stats::Box.test(residuals, lag=lag, type=type,
                  fitdf=if (fitdf < lag) fitdf else 0)
}

#' @param object an object of class \code{"ARIMA"}.
#' @describeIn Arima \code{update()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
#'
update.ARIMA <- function(object, formula, ...){
  cl0 <- object$call
  cl <- match.call()
  args <- object$dots
  args[c("data", "order", "seasonal")] <-
    object[c("data", "order", "seasonal")]
  dots <- list(...)
  names <- names(dots)
  for (name in names){
    args[[name]]  <- dots[[name]]
  }
  args$formula <- if (missing(formula)) {
    object$formula
    } else {
      update(object$formula, formula)
    }
  result <- do.call(Arima, args)
  names <- names(cl)
  for (name in names){
    if (name == "" || name == "object") next
    cl0[[name]] <- cl[[name]]
  }
  cl0[["formula"]] <- args$formula
  result$call <- cl0
  result
}

#' @describeIn Arima \code{coef()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
coef.ARIMA <- function(object, ...) {
  coefs <- coef(object$arima)
  names(coefs)[names(coefs) == "intercept"] <- "(Intercept)"
  coefs
}

#' @describeIn Arima \code{vcov()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
vcov.ARIMA <- function(object, ...) {
  vc <- vcov(object$arima)
  where.intercept <- rownames(vc) == "intercept"
  rownames(vc)[where.intercept] <- "(Intercept)"
  colnames(vc)[where.intercept] <- "(Intercept)"
  vc
}

#' @describeIn Arima \code{logLik()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
logLik.ARIMA <- function(object, ...) logLik(object$arima)

#' @describeIn Arima \code{AIC()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
AIC.ARIMA <- function(object, ..., k){
  object$arima$aic
}

#' @describeIn Arima \code{model.matrix()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
model.matrix.ARIMA <- function(object, ...){
  object$model.matrix
}

#' @describeIn Arima \code{model.frame()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
model.frame.ARIMA <- function(formula, ...){
  formula$model
}

#' @returns the \code{\link{cv.ARIMA}()} method returns an object of class
#'   \code{c("cvOrdered", "cv")} (see  \code{\link{cvOrdered}()} and \code{\link{cv}()}).
#'
#'   \code{Arima()} returns an object of class \code{"ARIMA"} with the
#'   following components: \code{formula}, the model formula; \code{data},
#'   the data set to which the model was fit;
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

#' @param n.ahead number of future cases to predict.
#' @param newdata data frame with rows containing the
#' predictors (if any) for predicted future cases.
#' @param se.fit if \code{TRUE} (the default is \code{FALSE}), compute
#' the standard errors of the predictions.
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

#' @param model an object of class \code{"ARIMA"}
#' produced by the \code{Arima()} function.
#' @param mod an object of class \code{"ARIMA"}
#' produced by the \code{Arima()} function.
#' @param data time-series data frame with the data to which the model was fit;
#' can usually be inferred from the model; for \code{Arima()}, a time-series data
#' frame with data to which the model is to be fit. \code{data} must be
#' a time-series data frame, i.e., a \code{"ts_data_frame"} object, whose
#' columns are each \code{\link[stats]{ts}} objects with common times; if it
#' is not, an attempt will be made to coerce it to a  \code{"ts_data_frame"} object
#' by \code{as.ts.data.frame()}.
#' @param criterion function to compute the CV cost criterion
#' (default \code{\link{mse}}).
#' @param k number of folds, an integer \eqn{\gt 2}; the default value of \code{k} depends on
#' \code{fold.type};  if \code{k = "n"},
#' the first fold is determined by \code{begin.with} and subsequent folds
#' each contain a single case; that makes sense if \code{fold.type = "cumulative"}.
#' Ignored for the \code{ATC()} method (to match the generic).
#' @param reps ignored (to match \code{\link{cv}()} generic function).
#' @param seed ignored (to match \code{\link{cv}()} generic function).
#' @param fold.type  if \code{"window"} (the default), folds comprise a moving window
#' of \code{begin.with} cases; if \code{"cumulative"}, predict
#' the response for 1 or more cases after the i-th fold from the model fit to data
#' in the \code{i}th and all preceding folds; if \code{"preceding"}, predict using cases in
#' the \code{i}th fold only.
#' @param begin.with the number of cases in
#' the first fold; disregarded for \code{fold.type="preceding"}.
#' The remaining cases are divided among the subsequent
#' \code{k} - 1 folds.
#' @param lead how far ahead to predict (can be a vector of positive integers);
#' the default is \code{1}.
#' @param min.ahead the minimum number of "future" cases to include in
#' the last fold (default \code{1}); more than one case may be required for
#' some kinds of predictions for terms with data-dependent bases, such as those using
#' \code{\link[stats]{poly}()} (orthogonal polynomials).
#' @param criterion.name name of the CV criterion; can usually be
#' inferred from \code{criterion}.
#' @param details return fold-wise parameter estimates for cases in each fold after the first;
#' the default is \code{TRUE} for \eqn{k \le 10,000}.
#' @param ncores if \code{ncores} \eqn{> 1}, the computation is parallelized.
#' @describeIn Arima \code{cv()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
cv.ARIMA <- function(model,
                   data = model$data,
                   criterion = mse,
                   k = if (fold.type == "preceding") 10 else "n",
                   reps,
                   seed,
                   fold.type = c("window", "cumulative", "preceding"),
                   begin.with=max(25, ceiling(n/10)),
                   lead = 1L,
                   min.ahead = 1L,
                   criterion.name = deparse(substitute(criterion)),
                   details = n <= 1e4,
                   ncores = 1L,
                   ...) {

  f <- function(i, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i)
    indices.j <- fold(folds, i, predict=TRUE, lead=lead,
                      up.to.max.lead=TRUE)
    model.i <- suppressWarnings(update(model, data = window(data, start=indices.i["start"],
                                           end=indices.i["stop"])))
    fit.i <- as.vector(predict(model.i, n.ahead=max(lead),
                     newdata=window(data, start=indices.j["start"],
                                    end=indices.j["stop"])))
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
    indices.j <- fold(folds, i, predict=TRUE, lead=lead,
                      up.to.max.lead=TRUE)
    # bring Arima() into scope
    assign(model.function.name, model.function)
    # the following deals with a scoping issue that can
    #   occur with args passed via ... (which is saved in dots)
    predict.args <- c(list(
      object = suppressWarnings(update(model, data = window(data, start=indices.i["start"],
                                           end=indices.i["stop"]))),
      newdata = window(data, start=indices.j["start"],
                       end=indices.j["stop"]),
      n.ahead = max(lead)), dots)

    fit.i <- as.vector(do.call(predict, predict.args))

    list(
      fit.i = fit.i[lead],
      coef.i = coef(predict.args$object)
    )
  }

  if (!inherits(data, "ts_data_frame")) {
    data <- as.ts(data)
    message("Note: 'data' coerced to 'ts_data_frame'")
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
    min.ahead = min.ahead,
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
  yy <- as.vector(GetResponse(model))
  yyhat <- as.vector(fitted(model))
  result[["full crit"]] <- criterion(na.omit(yy), na.omit(yyhat))
  result[["lead"]] <- lead
  result
}

#' @describeIn Arima \code{as.ts()} method for \code{"data.frame"} objects;
#' produces an object of class \code{"ts_data_frame"} inheriting from \code{"data.frame"}.
#' @param start the beginning time (see \code{link[stats]{ts}()}).
#' @param end the ending time (can be inferred from the data,
#' see \code{link[stats]{ts}()}).
#' @param frequency the period of the time series (e.g., \code{4}
#' for quarterly data, see \code{link[stats]{ts}()}).
#' @exportS3Method stats::as.ts
as.ts.data.frame <- function(x, start=1, end, frequency=1, ...){
  which.ts <- which(sapply(x, is.ts))
  if (missing(start)){
    if (length(which.ts) > 0L) {
      tsp <- tsp(x[[min(which.ts)]])
    } else {
      z <- numeric(nrow(x))
      z <- ts(z)
      tsp <- tsp(z)
    }
  } else {
    z <- numeric(nrow(x))
    z <- if (missing(end)){
      ts(z, start=start, frequency=frequency)
    } else {
      ts(z, start=start, end=end, frequency=frequency)
    }
    tsp <- tsp(z)
  }
  nms <- names(x)
  for (j in 1L:ncol(x)){
    if (j %in% which.ts){
      if (!all(tsp == tsp(x[, j]))){
        stop("incompatible 'ts' definition for ", nms[j])
      }
      next
    }
    xx <- x[, j]
    tsp(xx) <- tsp
    class(xx) <- c("ts", class(xx))
    x[, j] <- xx
  }
  class(x) <- c("ts_data_frame", class(x))
  x
}

#' @describeIn Arima \code{as.ts()} method for \code{"ts_data_frame"} objects;
#' returns the object unchanged.
#' @exportS3Method stats::as.ts
as.ts.ts_data_frame <- function(x, ...) x

#' @describeIn Arima \code{window()} method for \code{"ts_data_frame"} objects;
#' returns a subset of the object, also a \code{"ts_data_frame"}, with correctly subsetted \code{"ts"}
#' objects as columns.
#' @exportS3Method stats::window
window.ts_data_frame <- function(x, start=NULL, end=NULL, ...){
  result <- lapply(x, function(z) window(z, start=start, end=end, ...))
  result <- as.data.frame(result)
  class(result) <- c("ts_data_frame", class(x))
  result
}

#' @describeIn Arima \code{time()} method for \code{"ts_data_frame"} objects;
#' returns the (common) times associated with the columns of the \code{"ts_data_frame"}.
#' @exportS3Method stats::time
time.ts_data_frame <- function(x, ...){
  time(x[, 1])
}


# the following method isn't exported and shadows
# stats::predict.Arima(), from which it is derived
predict.Arima <- function (object, n.ahead = 1L, newxreg = NULL, se.fit = TRUE,
                           ...){
  myNCOL <- function(x) if (is.null(x))
    0
  else NCOL(x)
  # the following modification to stats::predict.Arima()
  # handles newxreg data frames produced e.g. by
  # poly() and ns()
  if (!is.null(newxreg)) newxreg <- as.matrix(newxreg)
  rsd <- object$residuals
  xr <- object$call$xreg
  xreg <- if (!is.null(xr))
    eval.parent(xr)
  else NULL
  ncxreg <- myNCOL(xreg)
  if (myNCOL(newxreg) != ncxreg)
    stop("'xreg' and 'newxreg' have different numbers of columns")
  xtsp <- tsp(rsd)
  n <- length(rsd)
  arma <- object$arma
  coefs <- object$coef
  narma <- sum(arma[1L:4L])
  if (length(coefs) > narma) {
    if (names(coefs)[narma + 1L] == "intercept") {
      newxreg <- cbind(intercept = rep(1, n.ahead), newxreg)
      ncxreg <- ncxreg + 1L
    }
    xm <- if (narma == 0)
      drop(as.matrix(newxreg) %*% coefs)
    else drop(as.matrix(newxreg) %*% coefs[-(1L:narma)])
  }
  else xm <- 0
  if (arma[2L] > 0L) {
    ma <- coefs[arma[1L] + 1L:arma[2L]]
    if (any(Mod(polyroot(c(1, ma))) < 1))
      warning("MA part of model is not invertible")
  }
  if (arma[4L] > 0L) {
    ma <- coefs[sum(arma[1L:3L]) + 1L:arma[4L]]
    if (any(Mod(polyroot(c(1, ma))) < 1))
      warning("seasonal MA part of model is not invertible")
  }
  z <- KalmanForecast(n.ahead, object$model)
  pred <- ts(z[[1L]] + xm, start = xtsp[2L] + deltat(rsd),
             frequency = xtsp[3L])
  if (se.fit) {
    se <- ts(sqrt(z[[2L]] * object$sigma2), start = xtsp[2L] +
               deltat(rsd), frequency = xtsp[3L])
    list(pred = pred, se = se)
  }
  else pred
}
