#' Fit and Cross-Validate ARIMA Models
#'
#' \code{\link{cv}()} methods for time-series models
#' fit by \code{Arima()}, which provides a formula interface to the
#' \code{\link{arima}()}) function.
#'
#' @param model an object of class \code{"ARIMA"}
#' produced by the \code{Arima()} function.
#' @param mod an object of class \code{"ARIMA"}
#' produced by the \code{Arima()} function.
#' @param data data frame with the data to which the model was fit;
#' can usually be inferred from the model; for \code{Arima()}, a data
#' frame with data to which the model is to be fit.
#' @param criterion function to compute the CV cost criterion
#' (default \code{\link{mse}}).
#' @param k number of folds, an integer \eqn{\gt 2}; the default value of \code{k} depends on
#' \code{fold.type};  if \code{k = "n"},
#' the first fold is determined by \code{begin.with} and subsequent folds
#' each contain a single case; that makes sense if \code{fold.type = "cumulative"}.
#' Ignored for the \code{ATC()} method (to match the generic).
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
#' })
#' }

#' @param formula either a one-sided formula giving the response variable
#' in an ARIMA model with no predictors, or a two-sided formula with the
#' response on the left and terms for the predictors on the right (as in
#' a typical R regression model); if the timeseries is differenced, the intercept
#' is removed from the model even if the formula implies an intercept.
#' For the \code{\link[stats]{model.frame}()}
#' method (to match the generic), an \code{"ARIMA"} object.
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
#' @param ... further arguments to be passed to \code{\link{arima}()},
#' \code{Arima()}, or other functions and methods.
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
  if (length(which.int > 0)) {
    has.intercept <- TRUE
    x <- x[, -which.int, drop=FALSE]
  } else {
    has.intercept <- FALSE
  }
  result <- list(formula=formula, data=data, subset=subset,
                 na.action=na.action, order=order, seasonal=seasonal,
                 call=cl, model=mf, dots=dots)
  if (length(formula) == 2){
    if (!(ncol(x) == 1) && is.numeric(x[, 1]))
      stop("formula must specify a single response")
    result$arima <- stats::arima(x[, 1], order=order, seasonal=seasonal,
                                 ...)
    result$response <- x[, 1]
  } else {
    y <- model.response(mf, "numeric")
    if (!is.vector(y) && is.numeric(y))
      stop("formula must specify a single response")
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
    quants <- quantile(residuals, c(0, 0.25, 0.5, 0.75, 1))
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

#' @param y which diagnostic plots to display; the default is
#' \code{c("residuals", "portmanteau", "acf", "pacf")}, where
#' \code{"portmanteau"} is the Box-Pierce or Ljung-Box test of
#' the autocorrelations of the residuals.
#' @param test one of \code{"Box-Pierce"} or \code{"Ljung-Box"},
#' with the former the default.
#' @param max.lag the maximum lag for the portmanteau tests; the
#' default is the same as the maximum lag for the residuals autocorrelations
#' and partial autocorrelations.
#' @param xlab label for horizontal ("time") axis; defaults to
#' \code{"Time"}.
#' @param main title for diagnostic plots.
#' @param col color for points and lines.
#' @describeIn Arima \code{plot()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
plot.ARIMA <- function(x,
                       y=c("fitted", "residuals", "portmantequ", "acf", "pacf"),
                       test=c("Box-Pierce", "Ljung-Box"),
                       max.lag=min(floor(10*log10(n)),  n - 1L),
                       xlab="time",
                       main="Diagnostic Plots",
                       col="blue", ...){

  test <- match.arg(test)

  which.plots <- if (missing(y)){
    c("residuals", "portmanteau", "acf", "pacf")
  } else {
    match.arg(y, several.ok=TRUE)
  }

  residuals <- residuals(x)
  n <- length(residuals)

  mfrow <- n2mfrow(length(which.plots))
  save <- par(mfrow=mfrow, oma=c(0, 0, 1, 0))
  on.exit(par(save))

  if ("fitted" %in% which.plots){
    plot(fitted(x), xlab=xlab, ylab=expression(hat(y)),
       main = "Fitted Values", type="b", pch=16, col=col)
    grid(lty=2, col="gray")
  }

  if ("residuals" %in% which.plots){
    plot(residuals, xlab=xlab, ylab="residuals",
       main="Residuals", type="b", pch=16, col=col)
    grid(lty=2, col="gray")
  }

  if ("portmanteau" %in% which.plots){
    p <- numeric(max.lag)
    for (lag in 1:max.lag)
      p[lag] <- testArima(x, lag=lag, type=test)[["p.value"]]
    plot(p, xlab="Lag", ylab="p-value",
         main=paste(test, "Tests\nfor Residuals"))
    abline(h=0.05, lty=2, col=col)
  }

  if ("acf" %in% which.plots)
    acf(residuals, main="Autocorrelations\n of Residuals",
      na.action=na.pass, ci.col=col)

  if ("pacf" %in% which.plots)
    pacf(residuals, main="Parial Autocorrelations\n of Residuals",
       na.action=na.pass, ci.col=col)

  title(main=main, outer=TRUE)
}

#' @importFrom car Anova linearHypothesis
#' @describeIn Arima \code{\link[car]{Anova}()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
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
#' @export
linearHypothesis.ARIMA <- function(model, hypothesis.matrix, rhs=NULL, ...){
  NextMethod(vcov. = vcov(model), coef. = coef(model),
             suppress.vcov.msg = TRUE)
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
#' @importFrom stats Box.test logLik vcov AIC cov2cor pnorm quantile terms
#' @importFrom grDevices n2mfrow
#' @export
testArima.ARIMA <- function(model, lag=1,
                            type = c("Box-Pierce", "Ljung-Box"),
                            ...){
  type <- match.arg(type)
  residuals <- residuals(model)
  n <- length(residuals)
  # use same lag as acf() and pacf():
  fitdf <- length(coef(model))
  stats::Box.test(residuals, lag=lag, type=type,
                  fitdf=if (fitdf < lag) fitdf else 0)
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
#' the first fold; disregarded for \code{fold.type="preceding"}.
#' The remaining cases are divided among the subsequent
#' \code{k} - 1 folds.
#' @describeIn Arima \code{cv()} method for \code{"ARIMA"} objects
#' created by the \code{\link{Arima}()} function.
#' @export
cv.ARIMA <- function(model,
                   data = model$data,
                   criterion = mse,
                   k = if (fold.type == "preceding") 10 else "n",
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

