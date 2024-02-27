getLossFn <- cv:::getLossFn

cvCompute <- function(model, data=insight::get_data(model),
                       criterion=mse, k=10, reps=1, seed,
                       details = k <= 10,
                       confint, level=0.95,
                       ncores=1,
                       type="response",
                       start=FALSE,
                       f,
                       fPara=f,
                       localFunctions=list(),
                       localVariables=list(),
                       ...){

  # put function and variable args in the local environment
  env <- environment()
  environment(f) <- env
  environment(fPara) <- env
  fnNames <- names(localFunctions)
  for (i in seq_along(localFunctions)){
    assign(fnNames[i], localFunctions[[i]])
  }
  varNames <- names(localVariables)
  for (i in seq_along(localVariables)){
    assign(varNames[i], localVariables[[i]])
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
    if (missing(seed) || is.null(seed)) seed <- sample(1e6, 1L)
    set.seed(seed)
    message("R RNG seed set to ", seed)
  } else {
    if (reps > 1) stop("reps should not be > 1 for n-fold CV")
    if (!missing(seed) && !is.null(seed)) message("Note: seed ignored for n-fold CV")
    seed <- NULL
  }
  folds <- folds(n, k)
  indices <- if (n > k) sample(n, n)  else 1:n # permute cases
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

cv.default <- function(model, data=insight::get_data(model),
                       criterion=mse, k=10, reps=1, seed=NULL,
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
    indices.i <- cv::fold(folds, i)
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

  n <- nrow(data)

  cvCompute(model=model, data=data, criterion=criterion, k=k,
            reps=reps, seed=seed, details=details, confint=confint,
            level=level, ncores=ncores, type=type, start=start,
            f=f, fPara=fPara, ...)
}


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
    indices.i <- cv::fold(folds, i)
    b.i <- UpdateLM(indices.i)
    names(b.i) <- coef.names
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
  coef.names <- names(b)
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

  XXi <- chol2inv(model$qr$qr[1L:p, 1L:p, drop = FALSE])
  Xy <- t(X) %*% (w * y)

  cvCompute(model=model, data=data, criterion=criterion, k=k,
            reps=reps, seed=seed, details=details, confint=confint,
            level=level, ncores=ncores, method=method,
            f=f,
            localFunctions=list(UpdateLM=UpdateLM),
            localVariables=list(X=X, w=w, XXi=XXi, Xy=Xy, b=b, p=p,
                                coef.names=coef.names),
            ...)

}



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
    b.u <- as.vector(b.u)
    names(b.u) <- coef.names
    b.u
  }
  f <- function(i){
    # helper function to compute cv criterion for each fold
    indices.i <- cv::fold(folds, i)
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
  if (method == "hatvalues" && k !=n ) stop('method="hatvalues" available only when k=n')
  if (method == "exact"){
    result <- cv:::cv.default(model=model, data=data, criterion=criterion, k=k, reps=reps, seed=seed,
                         ncores=ncores, method=method, details=details, start=start, ...)
    if (inherits(result, "cv")) result$"criterion" <- deparse(substitute(criterion))
    return(result)
  } else if (method == "hatvalues") {
    y <- GetResponse(model)
    h <- hatvalues(model)
    if (any(abs(h - 1) < sqrt(.Machine$double.eps)))
      stop("there are hatvalues numerically equal to 1")
    yhat <- y - residuals(model, type="response")/(1 - h)
    cv <- criterion(y, yhat) # mean(mapply(criterion, y=y, yhat=yhat))
    result <- list(k="n", "CV crit" = cv, method=method,
                   "criterion" = deparse(substitute(criterion)))
    class(result) <- "cv"
    return(result)
  } else {
    b <- coef(model)
    coef.names <- names(b)
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
    indices <- if (n > k) sample(n, n)  else 1:n # permute cases
    yhat <- numeric(n)

    cvCompute(model=model, data=data, criterion=criterion, k=k,
              reps=reps, seed=seed, details=details, confint=confint,
              level=level, ncores=ncores, start=start, method=method,
              f=f,
              localFunctions=list(UpdateIWLS=UpdateIWLS, linkinv=linkinv),
              localVariables=list(X=X, w=w, z=z, XXi=XXi, b=b, p=p,
                                  coef.names=coef.names),
              ...)
  }
}

