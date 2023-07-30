#' Generic Function to Cross-Validate a Model-Selection Procedure
#'
#' @param procedure a model-selection procedure function that accepts the following arguments:
#' `data`, set to the data argument to `cvSelect()`;
#' `indices`, the indices of the rows of data defining the current fold; if missing,
#' the model-selection procedure is applied to the full data; other arguments, to be
#' passed via `...`. `procedure()` should return a two-element vector with the result
#' of applying a cross-validation criterion to the cases in
#' the current fold for the model deleting that fold, and to
#' all of the cases again for the model deleting the current fold.
#' When the indices argument is missing, procedure() returns the cross-validation criterion for all of the cases based on
#' the model fit to all of the cases.
#' @param data full data frame for model selection
#' @param k perform k-fold cross-validation (default is n-fold)
#' @param seed for R's random number generator
#' @param parallel do computations in parallel? (default is FALSE)
#' @param ncores number of cores to use for parallel computations
#'           (default is number of physical cores detected)
#' @param ... arguments to be passed to `procedure()`.
#' @returns a "cv" object with the cv criterion averaged across the folds,
#' the bias-adjusted averaged cv criterion,
#' the criterion applied to the model fit to the full data set,
#' and the initial value of R's RNG seed
#' @importFrom MASS stepAIC
#' @export
cvSelect <- function(procedure,
                     data,
                     k=nrow(data),
                     seed, parallel=FALSE,
                     ncores=parallelly::availableCores(logical=FALSE),
                     ...){
  n <- nrow(data)
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
  if (parallel && ncores > 1){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    result <- foreach(i = 1L:k, .combine=rbind) %dopar% {
    #  require(cv) # !!! temporary
      procedure(data, indices[starts[i]:ends[i]], ...)
    }
    stopCluster(cl)
  } else {
    result <- matrix(0, k, 2L)
    for (i in 1L:k){
      result[i, ] <- procedure(data, indices[starts[i]:ends[i]], ...)
    }
  }
  cv <- weighted.mean(result[, 1L], folds)
  cv.full <- procedure(data, ...)
  adj.cv <- cv + cv.full - weighted.mean(result[, 2L], folds)
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                 "k" = if (k == n) "n" else k, "seed" = seed)
  class(result) <- "cv"
  result
}

#' @describeIn cvSelect select a model using the `stepAIC()` function in the
#' *MASS* package.
#' @param data full data frame
#' @param indices indices of cases in data defining the current fold
#' @param model a regression model object fit to data
#' @param criterion a CV criterion function
#' @param k. the `k` argument to `stepAIC()` (note the period in `k.`).
#' @export
selectStepAIC <- function(data, indices,
                          model, criterion=mse, k.=2, ...){
  #
  #
  #
  #
  y <- getResponse(model)
  if (missing(indices)) {
    model.i <- MASS::stepAIC(model, trace=FALSE, ...)
    fit.o.i <- predict(model.i, newdata=data, type="response")
    return(criterion(y, fit.o.i))
  }
  model <- update(model, data=data[-indices, ])
  model.i <- MASS::stepAIC(model, trace=FALSE, k=k., ...)
  fit.o.i <- predict(model.i, newdata=data, type="response")
  fit.i <- fit.o.i[indices]
  c(criterion(y[indices], fit.i), criterion(y, fit.o.i))
}

# selectAllSubsets <- function(data, indices, formula, ...){
#   if (missing(indices)){
#     regsubsets(formula, data=data, nbest=1)
#   }
# }
