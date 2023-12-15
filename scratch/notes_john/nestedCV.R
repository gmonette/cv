nestedCV <- function(model, data, loss, k=10, reps=200,
                     seed, level=0.90, ...){
  UseMethod("nestedCV")
}

squaredErrorLoss <- function(y, yhat) (y - yhat)^2

nestedCV.default <- function(model, data=insight::get_data(model),
                             loss=squaredErrorLoss, k=10, reps=200,
                             seed, level=0.90, ...){

  innerCV <- function(j.omit){
    y <- getResponse(model)
    folds <- (1:k)[-j.omit]
    indices <- as.vector(unlist(mapply(function(s, e) s:e,
                                       starts[folds], ends[folds])))
    e.in <- numeric(length(indices))
    e.in.start <- 1
    for (j in folds){
      folds.j <- setdiff(folds, j)
      indices.j <- as.vector(unlist(mapply(function(s, e) s:e,
                                           starts[folds.j], ends[folds.j])))
      model.j <- update(model, data=data[indices.j, ])
      indices.j <- starts[j]:ends[j]
      yhat <- predict(model.j, newdata=data[indices.j, ])
      e.in.end <- e.in.start + length(yhat) - 1
      e.in[e.in.start:e.in.end] <- loss(y[names(yhat)], yhat)
      e.in.start <- e.in.end + 1
    }
    e.in
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
  if (missing(seed)) seed <- sample(1e6, 1L)
  set.seed(seed)
  message("R RNG seed set to ", seed)
  nk <-  n %/% k # number of cases in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  y <- getResponse(model)

  # ordinary cv:
  data <- data[sample(n, n), ] # permute cases
  e <- numeric(n)
  for(j in 1:k){
    indices.j <- as.vector(unlist(mapply(function(s, e) s:e,
                                         starts[-j], ends[-j])))
    mod <- update(model, data=data[indices.j, ])
    yhat <- predict(mod, data[-indices.j, ])
    # browser()
    e[-indices.j] <- loss(y[names(yhat)], yhat)
  }
  err.cv <- mean(e)
  se.cv <- sd(e)/sqrt(n)

  es <- numeric(sum(rep(n, k) - folds)*reps) # numeric()
  a <- numeric(reps*k)
  b <- numeric(reps*k)
  i.es.start <- 1
  i.ab <- 0
  for (r in 1:reps){
    data <- data[sample(n, n), ]
    for (j in 1:k){
      e.in <- innerCV(j)
      indices.j <- as.vector(unlist(mapply(function(s, e) s:e,
                                           starts[-j], ends[-j])))
      mod <- update(model, data[indices.j, ])
      indices.j <- starts[j]:ends[j]
      yhat <- predict(mod, newdata=data[indices.j, ])
      e.out <- loss(y[names(yhat)], yhat)
      i.ab <- i.ab + 1
      a[i.ab] <- (mean(e.in) - mean(e.out))^2
      b[i.ab] <-  var(e.out)/length(indices.j)
      i.es.end <- (i.es.start + length(e.in)) - 1
      es[i.es.start:i.es.end] <- e.in
      i.es.start <- i.es.end + 1
    }
  }
  mse <- mean(a) - mean(b)
  adj.mse <- ((k-1)/k)*mse
  if (adj.mse < se.cv^2) adj.mse <- se.cv^2
  if (adj.mse > k*se.cv^2) adj.mse <- k*se.cv^2
  err.ncv <- mean(es)
  bias <- (1 + (k - 2)/k)*(err.ncv - err.cv)
  halfwidth <- qnorm(1 - (1 - level)/2)*sqrt(mse)
  ci.lower.ncv <- err.ncv - bias - halfwidth
  ci.upper.ncv <- err.ncv - bias + halfwidth
  halfwidth <- qnorm(1 - (1 - level)/2)*se.cv
  ci.lower.cv <- err.cv - halfwidth
  ci.upper.cv <- err.cv + halfwidth
  result <- c(mse=mse, adj.mse=adj.mse, err.ncv=err.ncv, err.cv=err.cv, se.cv=se.cv,
    bias=bias, ci.lower.ncv=ci.lower.ncv, ci.upper.ncv=ci.upper.ncv,
    ci.lower.cv=ci.lower.cv, ci.upper.cv=ci.upper.cv, level=level,
    k=k, reps=reps)
  class(result) <- "nestedCV"
  result
}

summary.nestedCV <- function(object, digits=getOption("digits"), ...){
  cat("", paste0(object["k"], "-fold Nested Cross Validation with"),
      object["reps"], "Replications")
  cat("\n Nested CV estimate of error:",
      signif(object["err.ncv"], digits=digits))
  cat("\n Estimated MSE of NCV estimate:",
      signif(object["mse"], digits=digits),
      paste0("(RMSE = ", signif(sqrt(object["mse"]), digits=digits), ")"))
  cat("\n Adjusted estimated MSE of NCV estimate:",
      signif(object["adj.mse"], digits=digits),
      paste0("(RMSE = ", signif(sqrt(object["adj.mse"]), digits=digits), ")"))
  cat("\n CV estimate of error:",
      signif(object["err.cv"], digits=digits))
  cat("\n Std. error of CV estimate of error:",
      signif(object["se.cv"], digits=digits))
  cat("\n Estimated bias of NCV estimate of error:",
      signif(object["bias"], digits=digits))
  cat(" \n", paste0(round(100*object["level"]), "%"),
      "conf. int. for NCV estimate:",
      paste0("(", signif(object["ci.lower.ncv"], digits=digits),
             ", ", signif(object["ci.upper.ncv"], digits=digits), ")"))
  cat(" \n", paste0(round(100*object["level"]), "%"),
      "conf. int. for CV estimate:",
      paste0("(", signif(object["ci.lower.cv"], digits=digits),
             ", ", signif(object["ci.upper.cv"], digits=digits), ")\n"))
}

print.nestedCV <- function(x, digits=getOption("digits"), ...){
  cat("", paste0(x["k"], "-fold Nested Cross Validation with"),
      x["reps"], "Replications")
  cat("\n Nested CV estimate of error:",
      signif(x["err.ncv"], digits=digits))
  invisible(x)
}

nestedCV.lm <- function(model, data=insight::get_data(model),
                             loss=squaredErrorLoss, k=10, reps=200,
                             seed, level=0.90, ...){

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

  innerCV <- function(j.omit){
    y <- getResponse(model)
    folds <- (1:k)[-j.omit]
    indices <- as.vector(unlist(mapply(function(s, e) s:e,
                                       starts[folds], ends[folds])))
    e.in <- numeric(length(indices))
    e.in.start <- 1
    for (j in folds){
      indices.j <- starts[j]:ends[j]
      omit <- c(starts[j.omit]:ends[j.omit], indices.j)
      b.j <- UpdateLM(omit)
      yhat <- (X[indices.j, ] %*% b.j)[, 1]
      e.in.end <- e.in.start + length(yhat) - 1
      e.in[e.in.start:e.in.end] <- loss(y[names(yhat)], yhat)
      e.in.start <- e.in.end + 1
    }
    e.in
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
  if (missing(seed)) seed <- sample(1e6, 1L)
  set.seed(seed)
  message("R RNG seed set to ", seed)

  nk <-  n %/% k # number of cases in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold

  # ordinary cv:
  data <- data[order <- sample(n, n), ] # permute cases
  model <- update(model, data=data)
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

  X <- model.matrix(model)
  y <- getResponse(model)
  w <- weights(model)
  if (is.null(w)) w <- rep(1, length(y))
  XXi <- chol2inv(model$qr$qr[1L:p, 1L:p, drop = FALSE])
  Xy <- t(X) %*% (w * y)

  e <- numeric(n)
  for(j in 1:k){
    indices.j <- starts[j]:ends[j]
    b.j <- UpdateLM(indices.j)
    yhat <- (X[indices.j, ] %*% b.j)[, 1]
    e[indices.j] <- loss(y[names(yhat)], yhat)
  }
  err.cv <- mean(e)
  se.cv <- sd(e)/sqrt(n)

  es <- numeric(sum(rep(n, k) - folds)*reps) # numeric()
  a <- numeric(reps*k)
  b <- numeric(reps*k)
  i.es.start <- 1
  i.ab <- 0
  for (r in 1:reps){
    data <- data[sample(n, n), ]
    model <- update(model, data=data)

    X <- model.matrix(model)
    y <- getResponse(model)
    w <- weights(model)
    if (is.null(w)) w <- rep(1, length(y))
    XXi <- chol2inv(model$qr$qr[1L:p, 1L:p, drop = FALSE])
    Xy <- t(X) %*% (w * y)

    for (j in 1:k){
      e.in <- innerCV(j)
      indices.j <- as.vector(unlist(mapply(function(s, e) s:e,
                                           starts[-j], ends[-j])))
      mod <- update(model, data[indices.j, ])
      indices.j <- starts[j]:ends[j]
      b.j <- UpdateLM(indices.j)
      yhat <- predict(mod, newdata=data[indices.j, ])
      e.out <- loss(y[names(yhat)], yhat)
      i.ab <- i.ab + 1
      a[i.ab] <- (mean(e.in) - mean(e.out))^2
      b[i.ab] <-  var(e.out)/length(indices.j)
      i.es.end <- (i.es.start + length(e.in)) - 1
      es[i.es.start:i.es.end] <- e.in
      i.es.start <- i.es.end + 1
    }
  }
  mse <- mean(a) - mean(b)
  adj.mse <- ((k-1)/k)*mse
  if (adj.mse < se.cv^2) adj.mse <- se.cv^2
  if (adj.mse > k*se.cv^2) adj.mse <- k*se.cv^2
  err.ncv <- mean(es)
  bias <- (1 + (k - 2)/k)*(err.ncv - err.cv)
  halfwidth <- qnorm(1 - (1 - level)/2)*sqrt(adj.mse)
  ci.lower.ncv <- err.ncv - bias - halfwidth
  ci.upper.ncv <- err.ncv - bias + halfwidth
  halfwidth <- qnorm(1 - (1 - level)/2)*se.cv
  ci.lower.cv <- err.cv - halfwidth
  ci.upper.cv <- err.cv + halfwidth
  result <- c(mse=mse, adj.mse=adj.mse, err.ncv=err.ncv, err.cv=err.cv, se.cv=se.cv,
              bias=bias, ci.lower.ncv=ci.lower.ncv, ci.upper.ncv=ci.upper.ncv,
              ci.lower.cv=ci.lower.cv, ci.upper.cv=ci.upper.cv, level=level,
              k=k, reps=reps)
  class(result) <- "nestedCV"
  result
}
