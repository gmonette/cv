#' Cross-Validate Regression Models
#'
#' A parallelized generic k-fold (including n-fold, i.e., leave-one-out)
#' cross-validation function, with a default method, and
#' specific methods for linear and generalized-linear models that can be much
#' more computationally efficient.
#'
#' @param model a regression model object (see Details).
#' @param data data frame to which the model was fit (not usually necessary).
#' @param criterion cross-validation criterion ("cost" or lack-of-fit) function of form \code{f(y, yhat)}
#'        where \code{y} is the observed values of the response and
#'        \code{yhat} the predicted values; the default is \code{\link{mse}}
#'        (the mean-squared error).
#' @param k perform k-fold cross-validation (default is \code{10}); \code{k}
#' may be a number or \code{"loo"} or \code{"n"} for n-fold (leave-one-out)
#' cross-validation; for \code{folds()}, \code{k} must be a number.
#' @param reps number of times to replicate k-fold CV (default is \code{1}).
#' @param confint if \code{TRUE} (the default if the number of cases is 400
#' or greater), compute a confidence interval for the bias-corrected CV
#' criterion, if the criterion is the average of casewise components.
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
#' @param method computational method to apply to a linear (i.e. \code{"lm"}) model
#' or to a generalized linear (i.e., \code{"glm"}) model. See Details for an explanation
#' of the available options.
#' @param type for the default method, value to be passed to the
#' \code{type} argument of \code{predict()};
#' the default is `type="response"`, which is appropriate, e.g., for a `"glm"` model
#' and may be recognized or ignored by \code{predict()} methods for other model classes.
#' @param start if \code{TRUE} (the default is \code{FALSE}), the \code{start} argument
#' set to the vector of regression coefficients for the model fit to the full data is passed
#' to \code{update()} possibly making the CV updates faster, e.g. for a GLM.
#' @param n number of cases, for constructed folds.
#' @param folds an object of class \code{"folds"}.
#' @param i a fold number for an object of class \code{"folds"}.
#' @param ... to match generic; passed to \code{predict()} for the default method.
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
#' \code{folds()} returns an object of class folds, for which
#' there are \code{fold()} and \code{print()} methods.
#'
#' @seealso \code{\link{cvMixed}}, \code{\link{cvSelect}}.
#'
#' @details
#' The default \code{cv()} method uses \code{\link{update}()} to refit the model
#' to each fold, and should work if there are appropriate \code{update()}
#' and \code{\link{predict}()} methods, and if the default method for \code{\link{GetResponse}()}
#' works or if a \code{GetResponse()} method is supplied.
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
#' which employs \code{update()}.
#'
#' There is also a method for robust linear models fit by
#' \code{\link[MASS]{rlm}()} in the \pkg{MASS} package (to avoid
#' inheriting the \code{"lm"} method for which the default \code{"auto"}
#' computational method would be inappropriate).
#'
#' For additional details, see the "Cross-validation of regression models"
#' vignette (\code{vignette("cv", package="cv")}).
#'
#' \code{cv()} is designed to be extensible to other classes of regression
#' models; see the "Extending the cv package" vignette
#' (\code{vignette("cv-extend", package="cv")}).
#'
#' @examples
#' data("Auto", package="ISLR2")
#' m.auto <- lm(mpg ~ horsepower, data=Auto)
#' cv(m.auto,  k="loo")
#' (cv.auto <- cv(m.auto, seed=1234))
#' compareFolds(cv.auto)
#' cv(m.auto, seed=1234, reps=3)
#'
#' data("Mroz", package="carData")
#' m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
#' cv(m.mroz, criterion=BayesRule, seed=123)
#'
#' data("Duncan", package="carData")
#' m.lm <- lm(prestige ~ income + education, data=Duncan)
#' m.rlm <- MASS::rlm(prestige ~ income + education,
#'                    data=Duncan)
#' cv(m.lm, k="loo", method="Woodbury")
#' cv(m.rlm, k="loo")
#'
#' (ffs <- folds(102, 5))
#' fold(ffs, 2)
#'
#' @export
cv <- function(model, data, criterion, k, reps=1, seed, ...){
  UseMethod("cv")
}

#' @describeIn cv \code{default} method
#' @importFrom stats coef family fitted lm.wfit lsfit model.frame
#' model.matrix model.response predict qnorm
#' update weighted.mean weights
#' residuals hatvalues printCoefmat sd
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lme4 lmer
#' @importFrom nlme lme
#' @importFrom MASS rlm
#' @importFrom methods functionBody
#' @export
cv.default <- function(model, data=insight::get_data(model),
                       criterion=mse, k=10, reps=1, seed,
                       details = k <= 10,
                       confint = n >= 400, level=0.95,
                       ncores=1,
                       type="response",
                       start=FALSE, ...){

  f <- function(i){
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i)
    model.i <- if (start) {
      update(model, data=data[ - indices.i, ], start=b)
    } else {
      update(model, data=data[ - indices.i, ])
      }
    fit.all.i <- predict(model.i, newdata=data, type=type, ...)
    fit.i <- fit.all.i[indices.i]
    list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
         coef.i=coef(model.i))
  }

  fPara <- function(i){
    # helper function to compute cv criterion for each fold
    #  with parallel computations
    indices.i <- fold(folds, i)
    # the following deals with a scoping issue that can
    #   occur with args passed via ... (which is saved in dots)
    predict.args <- c(list(
      object= if (start) {
          update(model, data=data[ - indices.i, ], start=b)
        } else {
          update(model, data=data[ - indices.i, ])
        },
      newdata=data, type=type), dots)
    fit.all.i <- do.call(predict, predict.args)
    fit.i <- fit.all.i[indices.i]
    list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
         coef.i=coef(predict.args$object))
  }

  y <- GetResponse(model)
  b <- coef(model)
  n <- nrow(data)
  if (is.character(k)){
    if (k == "n" || k == "loo") {
      k <- n
    }
  }
  if (!is.numeric(k) || length(k) > 1L || k > n || k < 2 || k != round(k)){
    stop('k must be an integer between 2 and n or "n" or "loo"')
  }
  if (k != n){
    if (missing(seed)) seed <- sample(1e6, 1L)
    set.seed(seed)
    message("R RNG seed set to ", seed)
  } else {
    if (reps > 1) stop("reps should not be > 1 for n-fold CV")
    if (!missing(seed) && !is.null(seed)) message("Note: seed ignored for n-fold CV")
    seed <- NULL
  }
  folds <- folds(n, k)
  yhat <- if (is.factor(y)){
    factor(rep(NA, n), levels=levels(y))
  } else if (is.character(y)) {
    character(n)
  } else {
    numeric(n)
  }

  if (details){
    crit.i <- numeric(k)
    coef.i <- vector(k, mode="list")
    names(crit.i) <- names(coef.i) <- paste("fold", 1:k, sep=".")
  } else {
    crit.i <- NULL
    coef.i <- NULL
  }

  if (ncores > 1){
    dots <- list(...)
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    result <- foreach(i = 1L:k) %dopar% {
      fPara(i)
    }
    stopCluster(cl)
    for (i in 1L:k){
      yhat[fold(folds, i)] <- result[[i]]$fit.i
      if (details){
        crit.i[i] <- criterion(y[fold(folds, i)],
                               yhat[fold(folds, i)])
        coef.i[[i]] <- result[[i]]$coef.i
      }
    }
  } else {
    result <- vector(k, mode="list")
    for (i in 1L:k){
      result[[i]] <- f(i)
      yhat[fold(folds, i)] <- result[[i]]$fit.i
      if (details){
        crit.i[i] <- criterion(y[fold(folds, i)],
                               yhat[fold(folds, i)])
        coef.i[[i]] <- result[[i]]$coef.i
      }
    }
  }
  cv <- criterion(y, yhat)
  cv.full <- criterion(y, predict(model, type=type, ...))
  loss <- getLossFn(cv) # casewise loss function
  if (!is.null(loss)) {
    adj.cv <- cv + cv.full -
      weighted.mean(sapply(result, function(x) x$crit.all.i), folds$folds)
    se.cv <- sd(loss(y, yhat))/sqrt(n)
    halfwidth <- qnorm(1 - (1 - level)/2)*se.cv
    ci <- if (confint) c(lower = adj.cv - halfwidth, upper = adj.cv + halfwidth,
                         level=round(level*100)) else NULL
  } else {
    adj.cv <- NULL
    ci <- NULL
  }
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv,
                 "full crit" = cv.full, "confint"=ci,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 "criterion" = deparse(substitute(criterion)),
                 "details"=list(criterion=crit.i,
                              coefficients=coef.i))
  class(result) <- "cv"
  if (reps == 1) {
    return(result)
  } else {
    res <- cv(model=model, data=data, criterion=criterion,
              k=k, ncores=ncores, reps=reps - 1,
              details=details, ...)
    if (reps  > 2){
      res[[length(res) + 1]] <- result
    } else {
      res <- list(res, result)
    }
    for (i in 1:(length(res) - 1)){
      res[[i]]["criterion"] <- res[[length(res)]]["criterion"]
    }
    class(res) <- "cvList"
    return(res)
  }
}

#' @describeIn cv \code{print()} method
#' @param x a \code{"cv"} or \code{"cvList"} object to be printed
#' @param digits significant digits for printing,
#' default taken from the \code{"digits"} option
#' @export
print.cv <- function(x, digits=getOption("digits"), ...){
  rnd <- function(x){
    if (round(log10(x)) >= digits) round(x)
    else signif(x, digits)
  }
  cat(x[["k"]], "-Fold Cross Validation", sep="")
  if (!is.null(x[["clusters"]])){
    cat(" based on", x[["n clusters"]],
        paste0("{", paste(x[["clusters"]], collapse=", "), "}"),
        "clusters")
  }
  if (!is.null(x[["method"]])) cat("\nmethod:", x[["method"]])
  if (!is.null(x[["criterion"]]) && x[["criterion"]] != "criterion")
    cat("\ncriterion:", x[["criterion"]])
  if (is.null(x[["SD CV crit"]])){
    cat("\ncross-validation criterion =", rnd(x[["CV crit"]]))
  } else {
    cat("\ncross-validation criterion = ",
        rnd(x[["CV crit"]]), " (", rnd(x[["SD CV crit"]]), ")", sep="")
  }
  if (!is.null(x[["adj CV crit"]])){
    if (is.null(x[["SD adj CV crit"]])){
    cat("\nbias-adjusted cross-validation criterion =", rnd(x[["adj CV crit"]]))
    } else {
      cat("\nbias-adjusted cross-validation criterion = ",
          rnd(x[["adj CV crit"]]), " (", rnd(x[["SD adj CV crit"]]), ")", sep="")
    }
  }
  if (!is.null(x$confint)){
    cat(paste0("\n", x$confint["level"],
        "% CI for bias-adjusted CV criterion = (",
        rnd(x$confint["lower"]), ", ", rnd(x$confint["upper"]), ")")
    )
  }
  if (!is.null(x[["full crit"]]))
    cat("\nfull-sample criterion =", rnd(x[["full crit"]]), "\n")
  invisible(x)
}

#' @describeIn cv \code{print()} method
#' @export
print.cvList <- function(x, ...){
  xx <- x
  reps <- length(xx)
  names(xx) <- paste("Replicate", 1L:reps)
  xx$Average <- xx[[1L]]
  sumry <-   summarizeReps(xx)
  xx$Average[["CV crit"]] <- sumry[["CV crit"]]
  xx$Average[["adj CV crit"]] <- sumry[["adj CV crit"]]
  xx$Average[["SD CV crit"]] <- sumry[["SD CV crit"]]
  xx$Average[["SD adj CV crit"]] <- sumry[["SD adj CV crit"]]
  xx$Average$confint <- NULL
  for (rep in seq_along(xx)){
    cat("\n", names(xx)[rep], ":\n", sep="")
    print(xx[[rep]])
  }
  return(invisible(x))
}

#' @describeIn cv \code{"lm"} method
#' @export
cv.lm <- function(model, data=insight::get_data(model),
                  criterion=mse, k=10, reps=1, seed,
                  details = k <= 10,
                  confint = n >= 400, level=0.95,
                  method=c("auto", "hatvalues", "Woodbury", "naive"),
                  ncores=1, ...){
  UpdateLM <- function(omit){
    # compute coefficients with omit cases deleted
    #  uses the Woodbury matrix identity
    # <https://en.wikipedia.org/wiki/Woodbury_matrix_identity>
    x <- X[omit, , drop=FALSE]
    dg <- if (length(omit) > 1L) diag(1/w[omit]) else 1/w[omit]
    XXi.u <- XXi + (XXi %*% t(x) %*% solve(dg - x %*% XXi %*% t(x)) %*% x %*% XXi)
    b.u <- XXi.u %*% (Xy - t(X[omit, , drop=FALSE]) %*% (w[omit] * y[omit]))
    as.vector(b.u)
  }
  f <- function(i){
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i)
    b.i <- UpdateLM(indices.i)
    fit.all.i <- X %*% b.i
    fit.i <- fit.all.i[indices.i]
    list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
         coef.i=b.i)
  }
  X <- model.matrix(model)
  y <- GetResponse(model)
  w <- weights(model)
  if (is.null(w)) w <- rep(1, length(y))
  n <- nrow(data)
  if (is.character(k)){
    if (k == "n" || k == "loo") {
      k <- n
    }
  }
  if (!is.numeric(k) || length(k) > 1L || k > n || k < 2 || k != round(k)){
    stop('k must be an integer between 2 and n or "n" or "loo"')
  }
  if (k == n){
    if (reps > 1) stop("reps should not be > 1 for n-fold CV")
    if (!missing(seed) && !is.null(seed)) message("Note: seed ignored for n-fold CV")
    seed <- NULL
  }

  b <- coef(model)
  p <- length(b)
  if (p > model$rank) {
    message(paste0("The model has ", if (sum(is.na(b)) == 1L) "an ",
                   "aliased coefficient", if (sum(is.na(b)) > 1L) "s", ":"))
    print(b[is.na(b)])
    message("Aliased coefficients removed from the model")
    X <- X[, !is.na(b)]
    p <- ncol(X)
    model <- lm.wfit(X, y, w)
  }
  method <- match.arg(method)
  if (method == "hatvalues" && k !=n ) stop('method="hatvalues" available only when k=n')
  if (method == "auto"){
    method <- if (k == n) "hatvalues" else "Woodbury"
  }
  if (method == "naive") return(NextMethod())
  if (method == "hatvalues"){
    h <- hatvalues(model)
    if (any(abs(h - 1) < sqrt(.Machine$double.eps)))
      stop("there are hatvalues numerically equal to 1")
    yhat <- y - residuals(model)/(1 -h)
    cv <- criterion(y, yhat)
    result <- list(k="n", "CV crit" = cv, "method"=method,
                   "criterion" = deparse(substitute(criterion)))
    class(result) <- "cv"
    return(result)
  }

  if (details){
    names.b <- names(b)
    crit.i <- numeric(k)
    coef.i <- vector(k, mode="list")
    names(crit.i) <- names(coef.i) <- paste("fold", 1:k, sep=".")
  } else {
    crit.i <- NULL
    coef.i <- NULL
  }

  XXi <- chol2inv(model$qr$qr[1L:p, 1L:p, drop = FALSE])
  Xy <- t(X) %*% (w * y)
  if (!is.numeric(k) || length(k) > 1L || k > n || k < 2 || k != round(k)){
    stop("k must be an integer between 2 and n")
  }
  if (k != n){
    if (missing(seed)) seed <- sample(1e6, 1L)
    set.seed(seed)
    message("R RNG seed set to ", seed)
  } else {
    seed <- NULL
  }
  folds <- folds(n, k)
  yhat <- numeric(n)
  if (ncores > 1L){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    result <- foreach(i = 1L:k) %dopar% {
      f(i)
    }
    stopCluster(cl)
    for (i in 1L:k){
      yhat[fold(folds, i)] <- result[[i]]$fit.i
      if (details){
        crit.i[i] <- criterion(y[fold(folds, i)],
                               yhat[fold(folds, i)])
        b.i <- result[[i]]$coef.i
        names(b.i) <- names.b
        coef.i[[i]] <- b.i
      }
    }
  } else {
    result <- vector(k, mode="list")
    for (i in 1L:k){
      result[[i]] <- f(i)
      yhat[fold(folds, i)] <- result[[i]]$fit.i
      if (details){
        crit.i[i] <- criterion(y[fold(folds, i)],
                               yhat[fold(folds, i)])
        b.i <- result[[i]]$coef.i
        names(b.i) <- names.b
        coef.i[[i]] <- b.i
      }
    }
  }
  cv <- criterion(y, yhat)
  cv.full <- criterion(y, fitted(model))
  loss <- getLossFn(cv) # casewise loss function
  if (!is.null(loss)) {
    adj.cv <- cv + cv.full -
      weighted.mean(sapply(result, function(x) x$crit.all.i), folds$folds)
    se.cv <- sd(loss(y, yhat))/sqrt(n)
    halfwidth <- qnorm(1 - (1 - level)/2)*se.cv
    ci <- if (confint) c(lower = adj.cv - halfwidth, upper = adj.cv + halfwidth,
                         level=round(level*100)) else NULL
  } else {
    adj.cv <- NULL
    ci <- NULL
  }
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv,
                 "full crit" = cv.full, "confint"=ci,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 "method"=method,
                 "criterion" = deparse(substitute(criterion)),
                 "details"=list(criterion=crit.i,
                                coefficients=coef.i))
  class(result) <- "cv"
  if (reps == 1) {
    return(result)
  } else {
    res <- cv(model=model, data=data, criterion=criterion,
              k=k, ncores=ncores, method=method, reps=reps - 1,
              details=details, ...)
    if (reps  > 2){
      res[[length(res) + 1]] <- result
    } else {
      res <- list(res, result)
    }
    for (i in 1:(length(res) - 1)){
      res[[i]]["criterion"] <- res[[length(res)]]["criterion"]
    }
    class(res) <- "cvList"
    return(res)
  }
}

#' @describeIn cv \code{"glm"} method
#' @export
cv.glm <- function(model, data=insight::get_data(model),
                   criterion=mse, k=10, reps=1, seed,
                   details = k <= 10,
                   confint = n >= 400, level=0.95,
                   method=c("exact", "hatvalues", "Woodbury"),
                   ncores=1,
                   start=FALSE,
                   ...){
  UpdateIWLS <- function(omit){
    # compute coefficients with omit cases deleted
    #  uses the Woodbury matrix identity
    # <https://en.wikipedia.org/wiki/Woodbury_matrix_identity>
    x <- X[omit, , drop=FALSE]
    dg <- if (length(omit) > 1L) diag(1/w[omit]) else 1/w[omit]
    XXi.u <- XXi + (XXi %*% t(x) %*% solve(dg - x %*% XXi %*% t(x)) %*% x %*% XXi)
    b.u <- XXi.u %*% (Xz - t(X[omit, , drop=FALSE]) %*% (w[omit] * z[omit]))
    as.vector(b.u)
  }
  f <- function(i){
    # helper function to compute cv criterion for each fold
    indices.i <- fold(folds, i)
    b.i <- UpdateIWLS(indices.i)
    fit.all.i <- linkinv(X %*% b.i)
    fit.i <- fit.all.i[indices.i]
    list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
         coef.i=b.i)
  }
  n <- nrow(data)
  if (is.character(k)){
    if (k == "n" || k == "loo") {
      k <- n
    }
  }
  if (!is.numeric(k) || length(k) > 1L || k > n || k < 2 || k != round(k)){
    stop('k must be an integer between 2 and n or "n" or "loo"')
  }
  if (k == n){
    if (reps > 1) stop("reps should not be > 1 for n-fold CV")
    if (!missing(seed) && !is.null(seed)) message("Note: seed ignored for n-fold CV")
    seed <- NULL
  }
  method <- match.arg(method)
  if (k != n){
    if (missing(seed)) seed <- sample(1e6, 1L)
    set.seed(seed)
    if (method != "exact") message("R RNG seed set to ", seed)
  } else {
    seed <- NULL
  }
  if (method == "hatvalues" && k !=n ) stop('method="hatvalues" available only when k=n')
  if (method == "exact"){
    result <- cv.default(model=model, data=data, criterion=criterion, k=k, reps=reps, seed=seed,
               ncores=ncores, method=method, details=details, start=start, ...)
    if (inherits(result, "cv")) result$"criterion" <- deparse(substitute(criterion))
    return(result)
  } else if (method == "hatvalues") {
    y <- GetResponse(model)
    h <- hatvalues(model)
    if (any(abs(h - 1) < sqrt(.Machine$double.eps)))
      stop("there are hatvalues numerically equal to 1")
    yhat <- y - residuals(model, type="response")/(1 - h)
    cv <- criterion(y, yhat)
    result <- list(k="n", "CV crit" = cv, method=method,
                   "criterion" = deparse(substitute(criterion)))
    class(result) <- "cv"
    return(result)
  } else {
    b <- coef(model)
    p <- length(b)
    w <- weights(model, type="working")
    X <- model.matrix(model)
    y <- GetResponse(model)
    if (p > model$rank) {
      message(paste0("The model has ", if (sum(is.na(b)) == 1L) "an ",
                     "aliased coefficient", if (sum(is.na(b)) > 1L) "s", ":"))
      print(b[is.na(b)])
      message("Aliased coefficients removed from the model")
      X <- X[, !is.na(b)]
      p <- ncol(X)
    }
    eta <- predict(model)
    mu <-  fitted(model)
    z <- eta + (y - mu)/family(model)$mu.eta(eta)
    mod.lm <- lm.wfit(X, z, w)
    linkinv <- family(model)$linkinv
    XXi <- chol2inv(mod.lm$qr$qr[1L:p, 1L:p, drop = FALSE])
    Xz <- t(X) %*% (w * z)
    if (!is.numeric(k) || length(k) > 1L || k > n || k < 2 || k != round(k)){
      stop("k must be an integer between 2 and n")
    }
    folds <- folds(n, k)
    yhat <- numeric(n)

    if (details){
      names.b <- names(b)
      crit.i <- numeric(k)
      coef.i <- vector(k, mode="list")
      names(crit.i) <- names(coef.i) <- paste("fold", 1:k, sep=".")
    } else {
      crit.i <- NULL
      coef.i <- NULL
    }

    if (ncores > 1L){
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      result <- foreach(i = 1L:k) %dopar% {
        f(i)
      }
      stopCluster(cl)
      for (i in 1L:k){
        yhat[fold(folds, i)] <- result[[i]]$fit.i
        if (details){
          crit.i[i] <- criterion(y[fold(folds, i)],
                                 yhat[fold(folds, i)])
          b.i <- result[[i]]$coef.i
          names(b.i) <- names.b
          coef.i[[i]] <- b.i
        }
      }
    } else {
      result <- vector(k, mode="list")
      for (i in 1L:k){
        result[[i]] <- f(i)
        yhat[fold(folds, i)] <- result[[i]]$fit.i
        if (details){
          crit.i[i] <- criterion(y[fold(folds, i)],
                                 yhat[fold(folds, i)])
          b.i <- result[[i]]$coef.i
          names(b.i) <- names.b
          coef.i[[i]] <- b.i
        }
      }
    }
    cv <- criterion(y, yhat)
    cv.full <- criterion(y, fitted(model))
    loss <- getLossFn(cv) # casewise loss function
    if (!is.null(loss)) {
      adj.cv <- cv + cv.full -
        weighted.mean(sapply(result, function(x) x$crit.all.i), folds$folds)
      se.cv <- sd(loss(y, yhat))/sqrt(n)
      halfwidth <- qnorm(1 - (1 - level)/2)*se.cv
      ci <- if (confint) c(lower = adj.cv - halfwidth, upper = adj.cv + halfwidth,
                           level=round(level*100)) else NULL
    } else {
      adj.cv <- NULL
      ci <- NULL
    }
    result <- list("CV crit" = cv, "adj CV crit" = adj.cv,
                   "full crit" = cv.full, confint=ci,
                   "k" = if (k == n) "n" else k, "seed" = seed,
                   "method"=method,
                   "criterion" = deparse(substitute(criterion)),
                   "details"=list(criterion=crit.i,
                                  coefficients=coef.i))
    class(result) <- "cv"
    if (reps == 1) {
      return(result)
    } else {
      res <- cv(model=model, data=data, criterion=criterion,
                k=k, ncores=ncores, method=method, reps=reps - 1, ...)
      if (reps  > 2){
        res[[length(res) + 1]] <- result
      } else {
        res <- list(res, result)
      }
      for (i in 1:(length(res) - 1)){
        res[[i]]["criterion"] <- res[[length(res)]]["criterion"]
      }
      class(res) <- "cvList"
      return(res)
    }
  }
}

#' @describeIn cv \code{"rlm"} method (to avoid inheriting the \code{"lm"} method)
#' @export
cv.rlm <- function(model, data, criterion, k, reps = 1, seed, ...){
  result <- NextMethod(method="naive")
  result$method <- NULL
  result
}

#' @describeIn cv used internally by \code{cv()} methods (not for direct use);
#' exported to support new \code{cv()} methods.
#' @export
folds <- function(n, k){
  nk <-  n %/% k # number of cases in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  indices <- if (n > k) sample(n, n)  else 1:n # permute cases
  result <- list(n=n, k=k, folds=folds, starts=starts,
                 ends=ends, indices=indices)
  class(result) <- "folds"
  result
}

#' @describeIn cv to extract a fold from a \code{"folds"} object.
#' @export
fold <- function(folds, i, ...) UseMethod("fold")

#' @describeIn cv \code{fold()} method for \code{"folds"} objects.
#' @export
fold.folds <- function(folds, i, ...) folds$indices[folds$starts[i]:folds$ends[i]]

#' @describeIn cv \code{print()} method for \code{"folds"} objects.
#' @export
print.folds <- function(x, ...){
  if (x$k == x$n){
    cat("LOO:", x$k, "folds for", x$n, "cases")
    return(invisible(x))
  }
  cat(x$k, "folds of approximately", floor(x$n/x$k),
      "cases each")
  for (i in 1:min(x$k, 10)){
    cat("\n fold", paste0(i, ": "))
    fld <- fold(x, i)
    if (length(fld) <= 10)  cat(fld)
    else cat(fld[1:10], "...")
  }
  if (x$k > 10) cat("\n ...")
  cat("\n")
  invisible(x)
}


# not exported

summarizeReps <- function(x){
  CVcrit <- mean(sapply(x, function(x) x[["CV crit"]]))
  CVcritSD <- sd(sapply(x, function(x) x[["CV crit"]]))
  CVcritRange <- range(sapply(x, function(x) x[["CV crit"]]))
  if (!is.null(x[[1]][["adj CV crit"]])) {
    adjCVcrit <- mean(sapply(x, function(x) x[["adj CV crit"]]))
    adjCVcritSD <- sd(sapply(x, function(x) x[["adj CV crit"]]))
    adjCVcritRange <- range(sapply(x, function(x) x[["adj CV crit"]]))
  } else {
    adjCVcrit <- adjCVcritSD <- adjCVcritRange <- NULL
  }
  list("CV crit" = CVcrit, "adj CV crit" = adjCVcrit,
       "CV crit range" = CVcritRange,
       "SD CV crit" = CVcritSD, "SD adj CV crit" = adjCVcritSD,
       "adj CV crit range" = adjCVcritRange)
}

getLossFn <- function(cv){
  fn.body <- attr(cv, "casewise loss")
  if (is.null(fn.body)) return(NULL)
  eval(parse(text=paste0("function(y, yhat) {\n",
                         paste(fn.body, collapse="\n"),
                         "\n}")))
}

