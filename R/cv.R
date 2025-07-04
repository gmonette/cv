#' Cross-Validate Regression Models
#'
#' \code{cv()} is a parallelized generic k-fold (including n-fold, i.e., leave-one-out)
#' cross-validation function, with a default method,
#' specific methods for linear and generalized-linear models that can be much
#' more computationally efficient, and a method for robust linear models.
#' There are also \code{cv()} methods for \link[=cv.merMod]{mixed-effects models},
#' for \link[=cv.function]{model-selection procedures},
#' and for \link[=cv.modList]{several models fit to the same data},
#'  which are documented separately.
#'
#' @param model a regression model object (see Details).
#' @param data data frame to which the model was fit (not usually necessary).
#' @param criterion cross-validation criterion ("cost" or lack-of-fit) function of form \code{f(y, yhat)}
#'        where \code{y} is the observed values of the response and
#'        \code{yhat} the predicted values; the default is \code{\link{mse}}
#'        (the mean-squared error).
#' @param criterion.name a character string giving the name of the CV criterion function
#' in the returned \code{"cv"} object (not usually needed).
#' @param k perform k-fold cross-validation (default is \code{10}); \code{k}
#' may be a number or \code{"loo"} or \code{"n"} for n-fold (leave-one-out)
#' cross-validation.
#' @param reps number of times to replicate k-fold CV (default is \code{1}).
#' @param confint if \code{TRUE} (the default if the number of cases is 400
#' or greater), compute a confidence interval for the bias-corrected CV
#' criterion, if the criterion is the average of casewise components;
#' for \code{plot.cvList()}, whether to plot confidence intervals around the
#' biased-adjusted CV criterion, defaulting to \code{TRUE} and applicable only
#' if confidence intervals are included in the \code{"cv"} object.
#' @param level confidence level (default \code{0.95}).
#' @param seed for R's random number generator; optional, if not
#' supplied a random seed will be selected and saved; not needed
#' for n-fold cross-validation.
#' @param details if \code{TRUE} (the default if the number of
#' folds \code{k <= 10}), save detailed information about the value of the
#' CV criterion for the cases in each fold and the regression coefficients
#' with that fold deleted.
#' @param ncores number of cores to use for parallel computations
#'        (default is \code{1}, i.e., computations aren't done in parallel).
#' @param method computational method to apply to a linear (i.e., \code{"lm"}) model
#' or to a generalized linear (i.e., \code{"glm"}) model. See Details for an explanation
#' of the available options.
#' @param type for the default method, value to be passed to the
#' \code{type} argument of \code{predict()};
#' the default is `type="response"`, which is appropriate, e.g., for a \code{"glm"} model
#' and may be recognized or ignored by \code{predict()} methods for other model classes.
#' @param start if \code{TRUE} (the default is \code{FALSE}), the \code{start} argument
#' to \code{\link[stats]{update}()} is set to the vector of regression coefficients for the model fit
#' to the full data, possibly making the CV updates faster, e.g., for a GLM.
#' @param ... to match generic; passed to \code{predict()} for the default \code{cv()} method;
#' passed to the \code{\link[car]{Tapply}()} function in the \pkg{car} package for
#' \code{summary.cvDataFrame()}; passed to default \code{\link{plot}()} method for
#' \code{plot.cvList()} or \code{plot.cv()}.
#' @param model.function a regression function, typically for a new \code{cv()} method that
#' that calls \code{cv.default()} via \code{NextMethod()},
#' residing in a package that's not a declared dependency of the \pkg{cv} package,
#' e.g., \code{nnet::multinom}. It's usually not necessary to specify
#' \code{model.function} to make \code{cv.default()} work.
#'
#' @returns The \code{cv()} methods return an object of class \code{"cv"}, with the CV criterion
#' (\code{"CV crit"}), the bias-adjusted CV criterion (\code{"adj CV crit"}),
#' the criterion for the model applied to the full data (\code{"full crit"}),
#' the confidence interval and level for the bias-adjusted CV criterion (\code{"confint"}),
#' the number of folds (\code{"k"}), and the seed for R's random-number
#' generator (\code{"seed"}). If \code{details=TRUE}, then the returned object
#' will also include a \code{"details"} component, which is a list of two
#' elements: \code{"criterion"}, containing the CV criterion computed for the
#' cases in each fold; and \code{"coefficients"}, regression coefficients computed
#' for the model with each fold deleted.  Some methods may return a
#' subset of these components and may add additional information.
#' If \code{reps} > \code{1}, then an object of class \code{"cvList"} is returned,
#' which is literally a list of \code{"cv"} objects.
#'
#' @seealso \code{\link{cv.merMod}}, \code{\link{cv.function}},
#' \code{\link{cv.modList}}.
#'
#' @details
#' The default \code{cv()} method uses \code{\link[stats]{update}()} to refit the model
#' to each fold, and should work if there are appropriate \code{update()}
#' and \code{\link{predict}()} methods, and if the default method for \code{\link{GetResponse}()}
#' works or if a \code{\link{GetResponse}()} method is supplied. The model must, however,
#' work correctly with \code{update()}, and in particular not have variables
#' in the model formula that aren't in the data to which the model was fit: see the last example.
#'
#' The \code{"lm"} and \code{"glm"} methods can use much faster computational
#' algorithms, as selected by the \code{method} argument. The linear-model
#' method accommodates weighted linear models.
#'
#' For both classes of models, for the leave-one-out (n-fold) case, fitted values
#' for the folds can be computed from the hat-values via
#' \code{method="hatvalues"} without refitting the model;
#' for GLMs, this method is approximate, for LMs it is exact.
#'
#' Again for both classes of models, when more than one case is omitted
#' in each fold, fitted values may be obtained without refitting the
#' model by exploiting the Woodbury matrix identity via \code{method="Woodbury"}.
#' As for hatvalues, this method is exact for LMs and approximate for
#' GLMs.
#'
#' The default for linear models is \code{method="auto"},
#' which is equivalent to \code{method="hatvalues"} for n-fold cross-validation
#' and \code{method="Woodbury"} otherwise; \code{method="naive"} refits
#' the model via \code{update()} and is generally much slower. The
#' default for generalized linear models is \code{method="exact"},
#' which employs \code{update()}. This default is conservative, and
#' it is usually safe to use \code{method="hatvalues"} for n-fold CV
#' or \code{method="Woodbury"} for k-fold CV.
#'
#' There is also a method for robust linear models fit by
#' \code{\link[MASS]{rlm}()} in the \pkg{MASS} package (to avoid
#' inheriting the \code{"lm"} method for which the default \code{"auto"}
#' computational method would be inappropriate).
#'
#' For additional details, see the "Cross-validating regression models"
#' vignette (\code{vignette("cv", package="cv")}).
#'
#' \code{cv()} is designed to be extensible to other classes of regression
#' models; see the "Extending the cv package" vignette
#' (\code{vignette("cv-extend", package="cv")}).
#'
#' @examples
#' if (requireNamespace("ISLR2", quietly=TRUE)){
#' withAutoprint({
#' data("Auto", package="ISLR2")
#' m.auto <- lm(mpg ~ horsepower, data=Auto)
#' cv(m.auto,  k="loo")
#' summary(cv(m.auto,  k="loo"))
#' summary(cv.auto <- cv(m.auto, seed=1234))
#' compareFolds(cv.auto)
#' plot(cv.auto)
#' plot(cv.auto, what="coefficients")
#' summary(cv.auto.reps <- cv(m.auto, seed=1234, reps=3))
#' cvInfo(cv.auto.reps, what="adjusted CV criterion")
#' plot(cv.auto.reps)
#' plot(cv(m.auto, seed=1234, reps=10, confint=TRUE))
#' D.auto.reps <- as.data.frame(cv.auto.reps)
#' head(D.auto.reps)
#' summary(D.auto.reps, mse ~ rep + fold, include="folds")
#' summary(D.auto.reps, mse ~ rep + fold, include = "folds",
#'         subset = fold <= 5) # first 5 folds
#' summary(D.auto.reps, mse ~ rep, include="folds")
#' summary(D.auto.reps, mse ~ rep, fun=sd, include="folds")
#' })
#' } else {
#' cat("\n install 'ISLR2' package to run these examples\n")
#' }
#'
#' if (requireNamespace("carData", quietly=TRUE)){
#' withAutoprint({
#' data("Mroz", package="carData")
#' m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
#' summary(cv.mroz <- cv(m.mroz, criterion=BayesRule, seed=123))
#' cvInfo(cv.mroz)
#' cvInfo(cv.mroz, "adjusted")
#' cvInfo(cv.mroz, "confint")
#'
#' data("Duncan", package="carData")
#' m.lm <- lm(prestige ~ income + education, data=Duncan)
#' m.rlm <- MASS::rlm(prestige ~ income + education,
#'                    data=Duncan)
#' summary(cv(m.lm, k="loo", method="Woodbury"))
#' summary(cv(m.rlm, k="loo"))
#' })
#' } else {
#' cat("\n install 'carData' package to run these examples\n")
#' }
#'
#' # the following (due to Joshua Philipp Entrop)
#' # produces an error:
#' \dontrun{
#' data("Auto", package="ISLR2")
#' Auto$mpg_20 <- as.numeric(Auto$mpg < 20)
#' mlist <- lapply(
#'                 1:3,
#'                 \(p) glm(mpg_20  ~ poly(horsepower, p), data = Auto)
#' )
#' cv(
#'    models(mlist),
#'    data = Auto,
#'    seed = 2120)
#' }
#' @export
cv <- function(model, data, criterion, k, reps = 1L, seed, ...) {
  UseMethod("cv")
}

#' @describeIn cv \code{"default"} method.
#' @importFrom stats coef family fitted getCall
#' lm.wfit lsfit model.frame
#' model.matrix model.response predict qnorm
#' update weighted.mean weights
#' residuals hatvalues printCoefmat sd
#' na.pass
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lme4 lmer
#' @importFrom nlme lme
#' @importFrom MASS rlm
#' @importFrom methods functionBody
#' @export
cv.default <- function(model,
                       data = insight::get_data(model),
                       criterion = mse,
                       k = 10L,
                       reps = 1L,
                       seed = NULL,
                       criterion.name = deparse(substitute(criterion)),
                       details = k <= 10L,
                       confint = n >= 400L,
                       level = 0.95,
                       ncores = 1L,
                       type = "response",
                       start = FALSE,
                       model.function,
                       ...) {
  f <- function(i_, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i_)
    model.i <- if (start) {
      update(model, data = data[-indices.i,], start = b)
    } else {
      update(model, data = data[-indices.i,])
    }
    fit.all.i <- predict(model.i, newdata = data, type = type, ...)
    fit.i <- fit.all.i[indices.i]
    list(
      fit.i = fit.i,
      crit.all.i = criterion(y, fit.all.i),
      coef.i = coef(model.i)
    )
  }

  fPara <-
    function(i_,
             model.function = model.function,
             model.function.name,
             ...) {
      # helper function to compute cv criterion for each fold
      #  with parallel computations
      indices.i <- fold(folds, i_)
      if (!is.null(model.function))
        assign(model.function.name, model.function)
      # the following deals with a scoping issue that can
      #   occur with args passed via ... (which is saved in dots)
      predict.args <- c(list(
        object = if (start) {
          update(model, data = data[-indices.i,], start = b)
        } else {
          update(model, data = data[-indices.i,])
        },
        newdata = data,
        type = type
      ), dots)
      fit.all.i <- do.call(predict, predict.args)
      fit.i <- fit.all.i[indices.i]
      list(
        fit.i = fit.i,
        crit.all.i = criterion(y, fit.all.i),
        coef.i = coef(predict.args$object)
      )
    }

  if (missing(model.function)) {
    model.function <- try(getCall(model)[[1L]])
    if (inherits(model.function, "try-error")) {
      stop(
        "cv.default() cannot extract the call from the model\n",
        "try specifying the model.function argument"
      )
    }
    model.function.name <- as.character(model.function)
    if (make.names(model.function.name) != model.function.name) {
      stop(
        model.function.name,
        " is not a valid object name.\n",
        "try specifying the model.function argument\n",
        " to cv.default()"
      )
    }
    model.function <- eval(model.function)
  } else {
    model.function.name <- sub("^.*\\:\\:", "",
                               deparse(substitute(model.function)))
  }

  n <- nrow(data)
  cvCompute(
    model = model,
    data = data,
    criterion = criterion,
    criterion.name = criterion.name,
    k = k,
    reps = reps,
    seed = seed,
    details = details,
    confint = confint,
    level = level,
    ncores = ncores,
    type = type,
    start = start,
    f = f,
    fPara = fPara,
    model.function = model.function,
    model.function.name = model.function.name,
    ...
  )
}

#' @describeIn cv \code{"lm"} method.
#' @export
cv.lm <- function(model,
                  data = insight::get_data(model),
                  criterion = mse,
                  k = 10L,
                  reps = 1L,
                  seed = NULL,
                  details = k <= 10L,
                  confint = n >= 400L,
                  level = 0.95,
                  method = c("auto", "hatvalues", "Woodbury", "naive"),
                  ncores = 1L,
                  ...) {
  UpdateLM <- function(omit) {
    # compute coefficients with omit cases deleted
    #  uses the Woodbury matrix identity
    # <https://en.wikipedia.org/wiki/Woodbury_matrix_identity>
    x <- X[omit, , drop = FALSE]
    dg <- if (length(omit) > 1L)
      diag(1 / w[omit])
    else
      1 / w[omit]
    XXi.u <-
      XXi + (XXi %*% t(x) %*% solve(dg - x %*% XXi %*% t(x)) %*% x %*% XXi)
    b.u <-
      XXi.u %*% (Xy - t(X[omit, , drop = FALSE]) %*% (w[omit] * y[omit]))
    as.vector(b.u)
  }
  f <- function(i_, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i_)
    b.i <- UpdateLM(indices.i)
    names(b.i) <- coef.names
    fit.all.i <- X %*% b.i
    fit.i <- fit.all.i[indices.i]
    list(
      fit.i = fit.i,
      crit.all.i = criterion(y, fit.all.i),
      coef.i = b.i
    )
  }
  X <- model.matrix(model)
  y <- GetResponse(model)
  w <- weights(model)
  if (is.null(w))
    w <- rep(1, length(y))
  n <- nrow(data)
  if (is.character(k)) {
    if (k == "n" || k == "loo") {
      k <- n
    }
  }
  if (!is.numeric(k) ||
      length(k) > 1L || k > n || k < 2L || k != round(k)) {
    stop('k must be an integer between 2 and n or "n" or "loo"')
  }
  if (k == n) {
    if (reps > 1L)
      stop("reps should not be > 1 for n-fold CV")
    if (!missing(seed) &&
        !is.null(seed))
      message("Note: seed ignored for n-fold CV")
    seed <- NULL
  }

  b <- coef(model)
  p <- length(b)
  coef.names <- names(b)
  if (p > model$rank) {
    message(paste0("The model has ", if (sum(is.na(b)) == 1L)
      "an ",
      "aliased coefficient", if (sum(is.na(b)) > 1L)
        "s", ":"))
    print(b[is.na(b)])
    message("Aliased coefficients removed from the model")
    X <- X[,!is.na(b)]
    p <- ncol(X)
    model <- lm.wfit(X, y, w)
  }
  method <- match.arg(method)
  if (method == "hatvalues" &&
      k != n)
    stop('method="hatvalues" available only when k=n')
  if (method == "auto") {
    method <- if (k == n)
      "hatvalues"
    else
      "Woodbury"
  }
  if (method == "naive")
    return(NextMethod(model.function = NULL))
  if (method == "hatvalues") {
    h <- hatvalues(model)
    if (any(abs(h - 1) < sqrt(.Machine$double.eps)))
      stop("there are hatvalues numerically equal to 1")
    yhat <- y - residuals(model) / (1 - h)
    cv <- criterion(y, yhat)
    result <- list(
      k = "n",
      "CV crit" = cv,
      "method" = method,
      "criterion" = deparse(substitute(criterion))
    )
    class(result) <- "cv"
    return(result)
  }

  XXi <- chol2inv(model$qr$qr[1L:p, 1L:p, drop = FALSE])
  Xy <- t(X) %*% (w * y)

  cvCompute(
    model = model,
    data = data,
    criterion = criterion,
    criterion.name = deparse(substitute(criterion)),
    k = k,
    reps = reps,
    seed = seed,
    details = details,
    confint = confint,
    level = level,
    ncores = ncores,
    start = FALSE,
    method = method,
    f = f,
    locals = list(
      UpdateLM = UpdateLM,
      X = X,
      w = w,
      XXi = XXi,
      Xy = Xy,
      b = b,
      p = p,
      coef.names = coef.names
    ),
    ...
  )

}


#' @describeIn cv \code{"glm"} method.
#' @export
cv.glm <- function(model,
                   data = insight::get_data(model),
                   criterion = mse,
                   k = 10L,
                   reps = 1L,
                   seed = NULL,
                   details = k <= 10L,
                   confint = n >= 400L,
                   level = 0.95,
                   method = c("exact", "hatvalues", "Woodbury"),
                   ncores = 1L,
                   start = FALSE,
                   ...) {
  UpdateIWLS <- function(omit) {
    # compute coefficients with omit cases deleted
    #  uses the Woodbury matrix identity
    # <https://en.wikipedia.org/wiki/Woodbury_matrix_identity>
    x <- X[omit, , drop = FALSE]
    dg <- if (length(omit) > 1L)
      diag(1 / w[omit])
    else
      1 / w[omit]
    XXi.u <-
      XXi + (XXi %*% t(x) %*% solve(dg - x %*% XXi %*% t(x)) %*% x %*% XXi)
    b.u <-
      XXi.u %*% (Xz - t(X[omit, , drop = FALSE]) %*% (w[omit] * z[omit]))
    b.u <- as.vector(b.u)
    names(b.u) <- coef.names
    b.u
  }
  f <- function(i_, ...) {
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i_)
    b.i <- UpdateIWLS(indices.i)
    fit.all.i <- linkinv(X %*% b.i)
    fit.i <- fit.all.i[indices.i]
    list(
      fit.i = fit.i,
      crit.all.i = criterion(y, fit.all.i),
      coef.i = b.i
    )
  }
  n <- nrow(data)
  if (is.character(k)) {
    if (k == "n" || k == "loo") {
      k <- n
    }
  }
  if (!is.numeric(k) ||
      length(k) > 1L || k > n || k < 2L || k != round(k)) {
    stop('k must be an integer between 2 and n or "n" or "loo"')
  }
  if (k == n) {
    if (reps > 1L)
      stop("reps should not be > 1 for n-fold CV")
    if (!missing(seed) &&
        !is.null(seed))
      message("Note: seed ignored for n-fold CV")
    seed <- NULL
  }
  method <- match.arg(method)
  if (method == "hatvalues" &&
      k != n)
    stop('method="hatvalues" available only when k=n')
  if (method == "exact") {
    result <-
      cv.default(
        model = model,
        data = data,
        criterion = criterion,
        k = k,
        reps = reps,
        seed = seed,
        ncores = ncores,
        method = method,
        details = details,
        start = start,
        model.function = NULL,
        ...
      )
    if (inherits(result, "cv"))
      result$"criterion" <- deparse(substitute(criterion))
    return(result)
  } else if (method == "hatvalues") {
    y <- GetResponse(model)
    h <- hatvalues(model)
    if (any(abs(h - 1) < sqrt(.Machine$double.eps)))
      stop("there are hatvalues numerically equal to 1")
    yhat <- y - residuals(model, type = "response") / (1 - h)
    cv <- criterion(y, yhat)
    result <- list(
      k = "n",
      "CV crit" = cv,
      method = method,
      "criterion" = deparse(substitute(criterion))
    )
    class(result) <- "cv"
    return(result)
  } else {
    b <- coef(model)
    coef.names <- names(b)
    p <- length(b)
    w <- weights(model, type = "working")
    X <- model.matrix(model)
    y <- GetResponse(model)
    if (p > model$rank) {
      message(paste0(
        "The model has ",
        if (sum(is.na(b)) == 1L)
          "an ",
        "aliased coefficient",
        if (sum(is.na(b)) > 1L)
          "s",
        ":"
      ))
      print(b[is.na(b)])
      message("Aliased coefficients removed from the model")
      X <- X[,!is.na(b)]
      p <- ncol(X)
    }
    eta <- predict(model)
    mu <-  fitted(model)
    z <- eta + (y - mu) / family(model)$mu.eta(eta)
    mod.lm <- lm.wfit(X, z, w)
    linkinv <- family(model)$linkinv
    XXi <- chol2inv(mod.lm$qr$qr[1L:p, 1L:p, drop = FALSE])
    Xz <- t(X) %*% (w * z)
    if (!is.numeric(k) ||
        length(k) > 1L || k > n || k < 2L || k != round(k)) {
      stop("k must be an integer between 2 and n")
    }

    cvCompute(
      model = model,
      data = data,
      criterion = criterion,
      criterion.name = deparse(substitute(criterion)),
      k = k,
      reps = reps,
      seed = seed,
      details = details,
      confint = confint,
      level = level,
      ncores = ncores,
      start = start,
      method = method,
      f = f,
      locals = list(
        UpdateIWLS = UpdateIWLS,
        linkinv = linkinv,
        X = X,
        w = w,
        z = z,
        XXi = XXi,
        b = b,
        p = p,
        coef.names = coef.names
      ),
      ...
    )
  }
}

#' @describeIn cv \code{"rlm"} method (to avoid inheriting the \code{"lm"} method).
#' @export
cv.rlm <- function(model, data, criterion, k, reps = 1L, seed, ...) {
  result <- NextMethod(method = "naive", model.function = NULL)
  result$method <- NULL
  result
}


#' @describeIn cv \code{print()} method for \code{"cv"} objects.
#' @param x a \code{"cv"}, \code{"cvList"}, or \code{"cvDataFrame"}
#' object to be printed or coerced to a data frame.
#' @param digits significant digits for printing,
#' default taken from the \code{"digits"} option.
#' @exportS3Method base::print
print.cv <- function(x, digits = getOption("digits"), ...) {
  rnd <- function(x) {
    if (round(log10(x)) >= digits)
      round(x)
    else
      signif(x, digits)
  }
  if (!is.null(x[["criterion"]]) &&
      x[["criterion"]] != "criterion"){
    cat(paste0("cross-validation criterion (",
        x[["criterion"]],
        ") = ",
        rnd(x[["CV crit"]])),
        "\n")
  } else {
    cat("cross-validation criterion =", rnd(x[["CV crit"]]), "\n")
  }
  invisible(x)
}

#' @exportS3Method base::summary
#' @describeIn cv \code{summary()} method for \code{"cv"} objects.
#' @param x a \code{"cv"}, \code{"cvList"}, or \code{"cvDataFrame"}
#' object to be plotted or summarized.
summary.cv <- function(object, digits = getOption("digits"), ...) {
  rnd <- function(object) {
    if (round(log10(object)) >= digits)
      round(object)
    else
      signif(object, digits)
  }
  cat(object[["k"]], "-Fold Cross Validation", sep = "")
  if (!is.null(object[["clusters"]])) {
    cat(" based on", object[["n clusters"]],
        paste0("{", paste(object[["clusters"]], collapse = ", "), "}"),
        "clusters")
  }
  if (!is.null(object[["method"]]))
    cat("\nmethod:", object[["method"]])
  if (!is.null(object[["criterion"]]) && object[["criterion"]] != "criterion")
    cat("\ncriterion:", object[["criterion"]])
  if (is.null(object[["SD CV crit"]])) {
    cat("\ncross-validation criterion =", rnd(object[["CV crit"]]))
  } else {
    cat("\ncross-validation criterion = ",
        rnd(object[["CV crit"]]),
        " (",
        rnd(object[["SD CV crit"]]),
        ")",
        sep = "")
  }
  if (!is.null(object[["adj CV crit"]])) {
    if (is.null(object[["SD adj CV crit"]])) {
      cat("\nbias-adjusted cross-validation criterion =", rnd(object[["adj CV crit"]]))
    } else {
      cat("\nbias-adjusted cross-validation criterion = ",
          rnd(object[["adj CV crit"]]),
          " (",
          rnd(object[["SD adj CV crit"]]),
          ")",
          sep = "")
    }
  }
  if (!is.null(object$confint)) {
    cat(
      paste0(
        "\n",
        object$confint["level"],
        "% CI for bias-adjusted CV criterion = (",
        rnd(object$confint["lower"]),
        ", ",
        rnd(object$confint["upper"]),
        ")"
      )
    )
  }
  if (!is.null(object[["full crit"]]))
    cat("\nfull-sample criterion =", rnd(object[["full crit"]]))
  cat("\n")
  invisible(object)
}

#' @describeIn cv \code{print()} method for \code{"cvList"} objects.
#' @exportS3Method base::print
print.cvList <- function(x, ...) {
  xx <- x
  reps <- length(xx)
  names(xx) <- paste("Replicate", 1L:reps)
  xx$Average <- xx[[1L]]
  sumry <-   summarizeReps(xx)
  xx$Average[["CV crit"]] <- sumry[["CV crit"]]

  if (!is.null(xx[[1]][["criterion"]]) &&
      xx[[1]][["criterion"]] != "criterion"){
    cat(paste0("cross-validation criterion (",
               x[[1]][["criterion"]],
               ")\n"))
  } else {
    cat("cross-validation criterion\n")
  }

  for (rep in seq_along(xx)) {
    cat(names(xx)[rep], ": ", sep = "")
    cat(xx[[rep]][["CV crit"]], "\n")
  }
  return(invisible(x))
}

#' @describeIn cv \code{summary()} method for \code{"cvList"} objects.
#' @exportS3Method base::summary
summary.cvList <- function(object, ...) {
  xx <- object
  reps <- length(xx)
  names(xx) <- paste("Replicate", 1L:reps)
  xx$Average <- xx[[1L]]
  sumry <-   summarizeReps(xx)
  xx$Average[["CV crit"]] <- sumry[["CV crit"]]
  xx$Average[["adj CV crit"]] <- sumry[["adj CV crit"]]
  xx$Average[["SD CV crit"]] <- sumry[["SD CV crit"]]
  xx$Average[["SD adj CV crit"]] <- sumry[["SD adj CV crit"]]
  xx$Average$confint <- NULL
  for (rep in seq_along(xx)) {
    cat("\n", names(xx)[rep], ":\n", sep = "")
    summary(xx[[rep]])
  }
  return(invisible(object))
}

#' @param y to match the \code{\link{plot}()} generic function, ignored.
#' @param what for \code{plot()} methods, what to plot: for the \code{"cv"} method, either \code{"CV criterion"}
#' (the default), or \code{"coefficients"};
#' for the \code{"cvList"} method, either \code{"adjusted CV criterion"}
#' (the default if present in the \code{"cv"} object) or \code{"CV object"}.
#'
#' For \code{cvInfo()}, the information to extract from a \code{"cv"},
#' \code{"cvModList"}, or \code{"cvList"} object,
#' one of: \code{"CV criterion"}, \code{"adjusted CV criterion"},
#' \code{"full CV criterion"} (the CV criterion applied to the model fit to the
#' full data set), \code{"SE"} (the standard error of the adjusted CV criterion),
#' \code{"confint"} (confidence interval for the adjusted CV criterion),
#' \code{"k"}, (the number of folds), \code{"seed"} (the seed employed for
#' R's random-number generator), \code{"method"} (the computational method
#' employed, e.g., for a \code{"lm"} model object), or \code{"criterion name"}
#' (the CV criterion employed); not all of these elements may be present, in
#' which case \code{cvInfo()} would return \code{NULL}.
#'
#' Partial matching is supported, so, e.g., \code{cvInfo(cv-object, "adjusted")}
#' is equivalent to \code{cvInfo(cv-object, "adjusted CV criterion")}
#' @describeIn cv \code{plot()} method for \code{"cv"} objects.
#' @exportS3Method base::plot
plot.cv <- function(x, y, what=c("CV criterion", "coefficients"), ...){
  if (is.null(x$details))
    stop("no 'details' element in 'x', nothing to plot")
  what <- match.arg(what)
  if (what == "CV criterion"){
    cv <- x[["CV crit"]]
    cv.folds <- x$details$criterion
    k <- x$k
    ylim <- range(c(cv, cv.folds))
    criterion <- x$criterion
    if (criterion == "criterion") criterion <- NULL
    plot(1:k, cv.folds, axes=FALSE, frame.plot=TRUE,
         pch=16, col=3, xlab="Fold",
         ylab=paste0("CV criterion: ", criterion),
         ...)
    axis(1, at = 1:k)
    axis(2)
    abline(h=cv, lty=2, lwd=2, col=2)
    usr <- par("usr")
    x <- mean(usr[1:2])
    y <- usr[4] + 0.15*(usr[4] - usr[3])
    graphics::legend(x, y, legend="Overall CV Criterion",
           lwd=2, lty=2, col=2, xpd=TRUE, bty="n", xjust=0.5)
    return(invisible(NULL))
  } else {
    coefs <- as.data.frame(x, columns="coefficients")
    coefs.stacked <- utils::stack(coefs[, -1])
    coefs.stacked <- cbind(coefs$fold, coefs.stacked)
    names(coefs.stacked) <- c("Fold", "Coefficient", "coef.name")
    coefs.stacked$coef.name <- sub("coef.", "", coefs.stacked$coef.name)
    lattice::xyplot(Coefficient ~ Fold | coef.name, data=coefs.stacked,
                    xlim=c(1, x$k),
                    scales=list(relation= "free"),
                    panel = function(x, y, ...){
                      lattice::panel.points(x[x != 0], y[x != 0], ...)
                      lattice::llines(x=range(x), y=rep(y[x == 0], 2),
                             lty=2, lwd=2, col=palette()[2])
                    },
                    key=list(text=list(lab="Overall Coefficient:",
                                       col=palette()[1]),
                             lines=2,
                             col=palette()[2], lty=2, lwd=2)
    )
  }
}

#' @describeIn cv \code{plot()} method for \code{"cvList"} objects.
#' @exportS3Method base::plot
plot.cvList <- function(x, y,
                        what=c("adjusted CV criterion", "CV criterion"),
                        confint=TRUE, ...){
  reps <- length(x)
  what <- match.arg(what)
  what <- if (what == "adjusted CV criterion") "adj CV crit" else "CV crit"
  if (is.null(x[[1]]$"adj CV crit")) what <- "CV crit"
  cv <- sapply(x, function(x) x[[what]])
  ci <- if (confint && !is.null(x[[1]]$confint) && what == "adj CV crit") {
    sapply(x, function(x) x[["confint"]])
  } else {
    NULL
  }
  ylim <- if (!is.null(ci)){
    c(min(ci["lower", ]), max(ci["upper", ]))
  } else {
    range(cv)
  }
  criterion <- x[[1]]$criterion
  if (criterion == "criterion") criterion <- NULL
  ylab <- paste0(if(what == "CV crit") "CV criterion: "
                 else "Adjusted CV criterion: ", criterion)
  plot(1:reps, cv, xlab="Replication", ylab=ylab,
       pch=16, col=2, xlim = c(0.75, length(x) + 0.25), ylim=ylim,
       axes=FALSE, frame.plot=TRUE, ...)
  axis(1, at=1:reps)
  axis(2)
  if (!is.null(ci)){
    for (j in 1:reps){
      arrows(x0=j, y0=ci["lower", j], y1=ci["upper", j],
             angle=90, code=3, lwd=2, col=3, length=1/(2*reps))
    }
  }
  invisible(NULL)
}

#' @describeIn cv extract information from a \code{"cv"} object.
#' @export
cvInfo <- function(object, what, ...){
  UseMethod("cvInfo")
}

#' @rdname cv
#' @export
cvInfo.cv <- function(object,  what=c("CV criterion",
                                      "adjusted CV criterion",
                                      "full CV criterion",
                                       "confint", "SE", "k", "seed",
                                       "method", "criterion name"),
                         ...){
                          what <- match.arg(what)
                          elements <- c("CV crit", "adj CV crit", "full crit",
                                        "confint", "SE adj CV crit", "k",
                                        "seed", "method", "criterion")
                          names(elements) <- c("CV criterion",
                                               "adjusted CV criterion",
                                               "full CV criterion",
                                               "confint", "SE", "k", "seed",
                                               "method", "criterion name")
                          info <- object[[elements[[what]]]]
                          nms <- names(info)
                          attributes(info) <- NULL
                          names(info) <- nms
                          info
                         }

#' @rdname cv
#' @export
cvInfo.cvModList <- function(object,  what=c("CV criterion",
                                      "adjusted CV criterion",
                                      "full CV criterion",
                                      "confint", "SE", "k", "seed",
                                      "method", "criterion name"),
                      ...){
  sapply(object, cvInfo, what=what, ...=...)
}

#' @rdname cv
#' @export
cvInfo.cvList <- function(object,  what=c("CV criterion",
                                             "adjusted CV criterion",
                                             "full CV criterion",
                                             "confint", "SE", "k", "seed",
                                             "method", "criterion name"),
                             ...){
  result <- sapply(object, cvInfo, what=what, ...=...)
  if (is.matrix(result)) {
    colnames(result) <- paste0("rep.", seq_along(object))
  } else {
    names(result) <- paste0("rep.", seq_along(object))
  }
  result
}

#' @describeIn cv \code{as.data.frame()} method for \code{"cv"} objects.
#' @param row.names optional row names for the result,
#' defaults to \code{NULL}.
#' @param optional to match the \code{\link{as.data.frame}()} generic function;
#' if \code{FALSE} (the default is \code{TRUE}), then the names of the columns
#' of the returned data frame, including the names of coefficients,
#' are coerced to syntactically correct names.
#' @param rows the rows of the resulting data frame to retain: setting
#' \code{rows="cv"} retains rows pertaining to the overall CV result
#' (marked as "\code{fold 0}" ); setting \code{rows="folds"} retains
#' rows pertaining to individual folds 1 through k; the default is
#' \code{rows = c("cv", "folds")}, which retains all rows.
#' @param columns the columns of the resulting data frame to retain:
#' setting \code{columns="critera"} retains columns pertaining to CV
#' criteria; setting \code{columns="coefficients"} retains columns pertaining
#' to model coefficients (broadly construed); the default is
#' \code{columns = c("criteria", "coefficients")}, which retains both;
#' and the columns \code{"model"}, \code{"rep"}, and \code{"fold"}, if present,
#' are always retained.
#' @exportS3Method base::as.data.frame
as.data.frame.cv <- function(x,
                             row.names = NULL,
                             optional = TRUE,
                             rows = c("cv", "folds"),
                             columns = c("criteria", "coefficients"),
                             ...) {
  rows <- match.arg(rows, several.ok = TRUE)
  columns <- match.arg(columns, several.ok = TRUE)
  D <- data.frame(fold = 0,
                  criterion = x$"CV crit")
  if (!is.null(x$"adj CV crit")) {
    D <- cbind(D, adjusted.criterion = x$"adj CV crit")
  }
  if (!is.null(x$"full crit")) {
    D <- cbind(D, full.criterion = x$"full crit")
  }
  if (!is.null(x$confint)) {
    D <-
      cbind(
        D,
        confint.lower = x$confint[1L],
        confint.upper = x$confint[2L],
        se.cv = x$"SE adj CV crit"
      )
  }
  if (!is.null(x$coefficients)) {
    coefs <- x$coefficients
    if (!is.matrix(coefs)) {
      coef.names <- names(coefs)
      coef.names[coef.names == "(Intercept)"] <- "Intercept"
      coef.names <- paste0("coef.", coef.names)
      names(coefs) <- coef.names
      D <- cbind(D, t(coefs))
    }
  }
  if (!is.null(x$details)) {
    D2 <- data.frame(
      fold = seq_along(x$details$criterion),
      criterion = x$details$criterion
    )
    if (!is.null(x$details$coefficients)) {
      Ds <- lapply(x$details$coefficients, t)
      D3 <- do.call(Merge, Ds)
      colnames <- colnames(D3)
      colnames[colnames == "(Intercept)"] <- "Intercept"
      colnames <- paste0("coef.", colnames)
      colnames(D3) <- colnames
      D2 <- cbind(D2, D3)
    }
    if (nrow(D2) > 0) D <- Merge(D, D2)
  }
  criterion <- x$criterion
  if (!is.null(criterion)) {
    colnames(D)[which(colnames(D) == "criterion")] <- criterion
    colnames(D)[which(colnames(D) == "adjusted.criterion")] <-
      paste0("adjusted.", criterion)
    colnames(D)[which(colnames(D) == "full.criterion")] <-
      paste0("full.", criterion)
    colnames(D)[which(colnames(D) == "confint.lower")] <-
      paste0("adj.", criterion, ".lower")
    colnames(D)[which(colnames(D) == "confint.upper")] <-
      paste0("adj.", criterion, ".upper")
    colnames(D)[which(colnames(D) == "se.cv")] <-
      paste0("SE.adj.", criterion)
  }
  rownames(D) <- row.names

  if (!"cv" %in% rows)
    D <- D[D$fold != 0,]
  if (!"folds" %in% rows)
    D <- D[D$fold == 0,]

  coefs.cols <- grepl("^coef\\.", colnames(D))
  always.cols <- colnames(D) %in% c("fold", "model", "rep")
  criteria.cols <- !(coefs.cols | always.cols)
  if (!"coefficients" %in% columns)
    coefs.cols <- FALSE
  if (!"criteria" %in% columns)
    criteria.cols <- FALSE
  D <- D[, coefs.cols | always.cols | criteria.cols]

  if (!optional) names(D) <- make.names(names(D), unique = TRUE)

  class(D) <- c("cvDataFrame", class(D))
  D
}

#' @describeIn cv \code{as.data.frame()} method for \code{"cvList"} objects.
#' @exportS3Method base::as.data.frame
as.data.frame.cvList <- function(x, row.names = NULL, optional = TRUE,
                                 ...) {
  Ds <- lapply(x, as.data.frame, optional=TRUE, ...)
  for (i in seq_along(Ds)) {
    Ds[[i]] <- cbind(rep = i, Ds[[i]])
  }
  D <- do.call(Merge, Ds)
  rownames(D) <- row.names
  class(D) <- c("cvListDataFrame", "cvDataFrame", class(D))
  if (!optional) names(D) <- make.names(names(D), unique = TRUE)
  D
}

#' @describeIn cv \code{print()} method for \code{"cvDataFrame"} objects.
#' @exportS3Method base::print
print.cvDataFrame <- function(x,
                              digits = getOption("digits") - 2L,
                              ...) {
  NextMethod(digits = digits)
}

#' @describeIn cv \code{summary()} method for \code{"cvDataFrame"} objects.
#' @param object an object to summarize or a \code{"cv"}, \code{"cvModlist"},
#' or \code{"cvList"} object from which to extract information via \code{cvInfo()}.
#' @param formula of the form \code{some.criterion ~ classifying.variable(s)}
#' (see examples).
#' @param subset a subsetting expression; the default (\code{NULL})
#' is not to subset the \code{"cvDataFrame"} object.
#' @param fun summary function to apply, defaulting to \code{mean}.
#' @param include which rows of the \code{"cvDataFrame"} to include in the
#' summary. One of \code{"cv"} (the default), rows representing the overall CV
#' results; \code{"folds"}, rows for individual folds; \code{"all"}, all rows
#' (generally not sensible).
#' @exportS3Method base::summary
summary.cvDataFrame <- function(object,
                                formula,
                                subset = NULL,
                                fun = mean,
                                include = c("cv", "folds", "all"),
                                ...) {
  include <- match.arg(include)
  if (include == "cv") {
    object <- object[object$fold == 0,]
  } else if (include == "folds") {
    object <- object[object$fold != 0,]
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "subset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$data <- object
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)
  object <- eval(mf, parent.frame())
  car::Tapply(formula, fun = fun, data = object, ...)
}
