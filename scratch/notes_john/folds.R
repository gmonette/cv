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

cv.lm1 <- function(model, data=insight::get_data(model),
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
    indices.i <- indices[starts[i]:ends[i]]
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
  nk <-  n %/% k # number of cases in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  indices <- if (n > k) sample(n, n)  else 1:n # permute cases
  yhat <- numeric(n)
  if (ncores > 1L){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    result <- foreach(i = 1L:k) %dopar% {
      f(i)
    }
    stopCluster(cl)
    for (i in 1L:k){
      yhat[indices[starts[i]:ends[i]]] <- result[[i]]$fit.i
      if (details){
        crit.i[i] <- criterion(y[indices[starts[i]:ends[i]]],
                               yhat[indices[starts[i]:ends[i]]])
        b.i <- result[[i]]$coef.i
        names(b.i) <- names.b
        coef.i[[i]] <- b.i
      }
    }
  } else {
    result <- vector(k, mode="list")
    for (i in 1L:k){
      result[[i]] <- f(i)
      yhat[indices[starts[i]:ends[i]]] <- result[[i]]$fit.i
      if (details){
        crit.i[i] <- criterion(y[indices[starts[i]:ends[i]]],
                               yhat[indices[starts[i]:ends[i]]])
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
      weighted.mean(sapply(result, function(x) x$crit.all.i), folds)
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

cv.lm2 <- function(model, data=insight::get_data(model),
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
    indices.i <- fold(folds, i) # indices[starts[i]:ends[i]]
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
  # nk <-  n %/% k # number of cases in each fold
  # rem <- n %% k  # remainder
  # folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  # ends <- cumsum(folds) # end of each fold
  # starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  # indices <- if (n > k) sample(n, n)  else 1:n # permute cases
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
     # yhat[indices[starts[i]:ends[i]]] <- result[[i]]$fit.i
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
      weighted.mean(sapply(result, function(x) x$crit.all.i), folds)
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




ffs <- folds(102, 10)

fold <- function(folds, ...) UseMethod("fold")

fold.folds <- function(f, i) f$indices[f$starts[i]:f$ends[i]]
fold(ffs, 2)

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


