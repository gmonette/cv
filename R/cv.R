#' Cross-Validate Regression Models
#'
#' A parallelized generic k-fold (including n-fold, i.e., leave-one-out)
#' cross-validation function, with a default method, and
#' specific methods for linear and generalized-linear models that can be much
#' more computationally efficient.
#'
#' @param model a regression model object (see Details).
#' @param data data frame to which the model was fit (not usually necessary)
#' @param criterion cross-validation criterion function of form \code{f(y, yhat)}
#'        where \code{y} is the observed values of the response and
#'        \code{yhat} the predicted values; the default is \code{\link{mse}
#'        (the mean-squared error)
#' @param k perform k-fold cross-validation (default is \code{10}); \code{k}
#' may be a number or \code{"loo"} or \code{"n"} for n-fold (leave-one-out)
#' cross-validation.
#' @param seed for R's random number generator; optional, if not
#' supplied a random seed will be selected and saved; not needed
#' for n-fold cross-validation
#' @param parallel do computations in parallel? (default is \code{FALSE})
#' @param ncores number of cores to use for parallel computations
#'           (default is number of physical cores detected)
#' @param method computational method to apply to a linear (i.e. \code{"lm"}) model
#' or to a generalized linear (i.e., \code{"glm"}) model. See Details for an explanation
#' of the available options.
#' @param ... to match generic
#' @returns a \code{"cv"} object with the cross-validation criterion averaged across the folds,
#' the bias-adjusted averaged CV criterion,
#' the criterion applied to the model fit to the full data set,
#' and the initial value of R's RNG seed. Some methods return a
#' subset of these components and may add additional information.
#'
#' @details
#' The default method uses \code{\link{update}()} to refit the model
#' to each fold, and should work if there are appropriate \code{update()}
#' and \code{\link{predict}()} methods and if the default method for \code{\link{getResponse}()}
#' works or if a \code{getResponse()} method is supplied.
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
#' @examples
#' data("Auto", package="ISLR")
#' m.auto <- lm(mpg ~ horsepower, data=Auto)
#' cv(m.auto,  k="loo")
#' cv(m.auto, seed=1234)
#'
#' data("Caravan", package="ISLR")
#' m.caravan <- glm(Purchase ~ ., data=Caravan[1:2500, ], family=binomial)
#' cv(m.caravan, k=5, criterion=BayesRule, seed=123)
#' @export
cv <- function(model, data, criterion, k, seed, ...){
  UseMethod("cv")
}

#' @describeIn cv \code{default} method
#' @importFrom stats coef family fitted lm.wfit model.frame
#' model.matrix model.response predict update weighted.mean weights
#' residuals hatvalues
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lme4 lmer
#' @export
cv.default <- function(model, data=insight::get_data(model),
                       criterion=mse, k=10,
                       seed, parallel=FALSE,
                       ncores=parallelly::availableCores(logical=FALSE),
                       method=NULL, ...){
  f <- function(i){
    # helper function to compute cv criterion for each fold
    indices.i <- indices[starts[i]:ends[i]]
    model.i <- update(model, data=data[ - indices.i, ])
    fit.o.i <- predict(model.i, newdata=data, type="response")
    fit.i <- fit.o.i[indices.i]
    c(criterion(y[indices.i], fit.i), criterion(y, fit.o.i))
  }
  y <- getResponse(model)
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
    seed <- NULL
  }
  nk <-  n %/% k # number of cases in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  indices <- if (n > k) sample(n, n)  else 1:n # permute cases
  if (parallel && ncores > 1){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    result <- foreach(i = 1L:k, .combine=rbind) %dopar% {
      f(i)
    }
    stopCluster(cl)
  } else {
    result <- matrix(0, k, 2L)
    for (i in 1L:k){
      result[i, ] <- f(i)
    }
  }
  cv <- weighted.mean(result[, 1L], folds)
  cv.full <- criterion(y, fitted(model))
  adj.cv <- cv + cv.full - weighted.mean(result[, 2L], folds)
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 "method"=method)
  class(result) <- "cv"
  result
}

#' @describeIn cv \code{print()} method
#' @param x a \code{cv} object to be printed
#' @export
print.cv <- function(x, ...){
  cat(x[["k"]], "-Fold Cross Validation", sep="")
  if (!is.null(x[["clusters"]])){
    cat(" based on", x[["n clusters"]],
        paste0("{", paste(x[["clusters"]], collapse=" ,"), "}"),
        "clusters")
  }
  if (!is.null(x[["method"]])) cat("\nmethod:", x[["method"]])
  cat("\ncross-validation criterion =", x[["CV crit"]])
  if (!is.null(x[["adj CV crit"]]))
    cat("\nbias-adjusted cross-validation criterion =", x[["adj CV crit"]])
  if (!is.null(x[["full crit"]]))
    cat("\nfull-sample criterion =", x[["full crit"]], "\n")
  invisible(x)
}

#' @describeIn cv \code{"lm"} method
#' @export
cv.lm <- function(model, data=insight::get_data(model), criterion=mse, k=10,
                  seed, method=c("auto", "hatvalues", "Woodbury", "naive"),
                  parallel=FALSE,
                  ncores=parallelly::availableCores(logical=FALSE), ...){
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
    indices.i <- indices[starts[i]:ends[i]]
    b.i <- UpdateLM(indices.i)
    fit.o.i <- X %*% b.i
    fit.i <- fit.o.i[indices.i]
    c(criterion(y[indices.i], fit.i), criterion(y, fit.o.i))
  }
  X <- model.matrix(model)
  y <- getResponse(model)
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
    yhat <- y - residuals(model)/(1 - hatvalues(model))
    cv <- criterion(y, yhat)
    result <- list(k="n", "CV crit" = cv, "method"=method)
    class(result) <- "cv"
    return(result)
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
  nk <-  n %/% k # number of cases in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  indices <- if (n > k) sample(n, n)  else 1:n # permute cases
  if (parallel && ncores > 1L){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    result <- foreach(i = 1L:k, .combine=rbind) %dopar% {
      f(i)
    }
    stopCluster(cl)
  } else {
    result <- matrix(0, k, 2L)
    for (i in 1L:k){
      result[i, ] <- f(i)
    }
  }
  cv <- weighted.mean(result[, 1L], folds)
  cv.full <- criterion(y, fitted(model))
  adj.cv <- cv + cv.full - weighted.mean(result[, 2L], folds)
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 "method"=method)
  class(result) <- "cv"
  result
}

#' @describeIn cv \code{"glm"} method
#' @export
cv.glm <- function(model, data=insight::get_data(model), criterion=mse, k=10,
                   seed,
                   method=c("exact", "hatvalues", "Woodbury"),
                   parallel=FALSE,
                   ncores=parallelly::availableCores(logical=FALSE),
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
    indices.i <- indices[starts[i]:ends[i]]
    b.i <- UpdateIWLS(indices.i)
    fit.o.i <- linkinv(X %*% b.i)
    fit.i <- fit.o.i[indices.i]
    c(criterion(y[indices.i], fit.i), criterion(y, fit.o.i))
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
    cv.default(model=model, data=data, criterion=criterion, k=k, seed=seed,
               parallel=parallel, ncores=ncores, method=method, ...)
  } else if (method == "hatvalues") {
    y <- getResponse(model)
    yhat <- y - residuals(model, type="response")/(1 - hatvalues(model))
    cv <- criterion(y, yhat)
    result <- list(k="n", "CV crit" = cv, method="method")
    class(result) <- "cv"
    return(result)
  } else {
    b <- coef(model)
    p <- length(b)
    w <- weights(model, type="working")
    X <- model.matrix(model)
    y <- getResponse(model)
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
    nk <-  n %/% k # number of cases in each fold
    rem <- n %% k  # remainder
    folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
    ends <- cumsum(folds) # end of each fold
    starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
    indices <- if (n > k) sample(n, n)  else 1:n # permute cases
    if (parallel && ncores > 1L){
      cl <- makeCluster(ncores)
      registerDoParallel(cl)
      result <- foreach(i = 1L:k, .combine=rbind) %dopar% {
        f(i)
      }
      stopCluster(cl)
    } else {
      result <- matrix(0, k, 2L)
      for (i in 1L:k){
        result[i, ] <- f(i)
      }
    }
    cv <- weighted.mean(result[, 1L], folds)
    cv.full <- criterion(y, fitted(model))
    adj.cv <- cv + cv.full - weighted.mean(result[, 2L], folds)
    result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                   "k" = if (k == n) "n" else k, "seed" = seed,
                   "method"=method)
    class(result) <- "cv"
    result
  }
}
