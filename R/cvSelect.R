#' Cross-Validate a Model-Selection Procedure
#'
#' A generic function to cross-validate a model-selection procedure,
#' along with a procedure that applies the \code{\link[MASS]{stepAIC}()}
#' model-selection function in the \pkg{MASS} package.
#'
#'
#' @param procedure a model-selection procedure function (see Details).
#' @param data full data frame for model selection.
#' @param k perform k-fold cross-validation (default is 10); \code{k}
#' may be a number or \code{"loo"} or \code{"n"} for n-fold (leave-one-out)
#' cross-validation.
#' @param seed for R's random number generator; not used for n-fold cross-validation.
#' @param parallel do computations in parallel? (default is \code{FALSE}),
#' @param ncores number of cores to use for parallel computations
#'           (default is number of physical cores detected).
#' @param ... arguments to be passed to \code{procedure()}.
#' @returns \code{cvSelect()} return a \code{"cv"} object with the CV criterion averaged across the folds,
#' the bias-adjusted averaged CV criterion,
#' the criterion applied to the model fit to the full data set,
#' and the initial value of R's RNG seed.
#' @importFrom MASS stepAIC
#' @describeIn cvSelect apply cross-validation to a model-selection procedure.
#'
#' @details
#' The model-selection function supplied as the \code{procedure} argument
#' to \code{cvSelect()} should accept the following arguments:
#' \describe{
#'  \item{\code{data}}{set to the \code{data} argument to \code{cvSelect()}.}
#'  \item{\code{indices}}{the indices of the rows of \code{data} defining the current fold; if missing,
#'  the model-selection procedure is applied to the full \code{data}.}
#'   \item{other arguments}{to be passed via the \code{...}
#'   to \code{cvSelect()}.}
#' }
#' \code{procedure()} should return a two-element vector with the result
#' of applying a cross-validation criterion to the cases in
#' the current fold for the model deleting that fold, and to
#' all of the cases again for the model deleting the current fold.
#'
#' When the \code{indices} argument is missing, \code{procedure()} returns the cross-validation criterion for all of the cases based on
#' the model fit to all of the cases.
#'
#' For an example of a model-selection function for the \code{procedure}
#' argument, see the code for \code{selectStepAIC()}.
#'
#' @seealso \code{\link[MASS]{stepAIC}()}
#'
#' @export
cvSelect <- function(procedure,
                     data,
                     k=10,
                     seed, parallel=FALSE,
                     ncores=parallelly::availableCores(logical=FALSE),
                     ...){
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
    result <- foreach(i = 1L:k, .combine=rbind) %dopar%
      procedure(data, indices[starts[i]:ends[i]], ...)
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

#' @describeIn cvSelect select a model using the \code{\link[MASS]{stepAIC}()} function in the
#' \pkg{MASS} package.
#' @param indices indices of cases in data defining the current fold.
#' @param model a regression model object fit to data.
#' @param criterion a CV criterion function.
#' @param k. the \code{k} argument to \code{\link[MASS]{stepAIC}()} (note the period in \code{k.}).
#' @examples
#' data("Auto", package="ISLR2")
#' m.auto <- lm(mpg ~ . - name - origin, data=Auto)
#' cvSelect(selectStepAIC, Auto, seed=123, model=m.auto)
#' cvSelect(selectStepAIC, Auto, seed=123, model=m.auto,
#'          k.=log(nrow(Auto))) # via BIC
#' @export
selectStepAIC <- function(data, indices,
                          model, criterion=mse, k.=2, ...){
  y <- getResponse(model)
  if (missing(indices)) {
    model.i <- MASS::stepAIC(model, trace=FALSE, k=k., ...)
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
