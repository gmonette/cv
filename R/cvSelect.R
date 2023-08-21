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
#' @param save.coef save the coefficients from the selected models? (default is \code{TRUE} if
#' \code{k} is 10 or smaller, \code{FALSE} otherwise)
#' @param reps number of times to replicate k-fold CV (default is \code{1})
#' @param seed for R's random number generator; not used for n-fold cross-validation.
#' @param ncores number of cores to use for parallel computations
#'        (default is \code{1}, i.e., computations aren't done in parallel)
#' @param ... arguments to be passed to \code{procedure()}.
#' @importFrom MASS stepAIC
#' @describeIn cvSelect apply cross-validation to a model-selection procedure.
#' @returns An object of class , inheriting from class \code{"cv"}, with the averaged CV criterion
#' (\code{"CV crit"}), the adjusted average CV criterion (\code{"adj CV crit"}),
#' the criterion for the model applied to the full data (\code{"full crit"}),
#' the number of folds (\code{"k"}), the seed for R's random-number
#' generator (\code{"seed"}), and (optionally) a list of coefficients for the selected models
#' for each fold (\code{"coefs"}).
#' #' If \code{reps} > \code{1}, then an object of class \code{c("cvSelectList", "cvList")} is returned,
#' which is literally a list of \code{c("cvSelect", "cv")} objects.
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
cvSelect <- function(procedure, data, k=10, reps=1,
                     save.coef = k <= 10,
                     seed, ncores=1, ...){
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
    seed <- NULL
  }
  nk <-  n %/% k # number of cases in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  indices <- if (n > k) sample(n, n)  else 1:n # permute cases
  if (ncores > 1){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    arglist <- c(list(data=data, indices=1, save.coef=save.coef),
                 list(...))
    selection <- foreach(i = 1L:k, .combine=c) %dopar% {
      # the following deals with a scoping issue that can
      #   occur with args passed via ...
      arglist$indices <- indices[starts[i]:ends[i]]
      selection <- do.call(procedure, arglist)
    }
    if (save.coef){
      is <- seq(1L:(2*k))
      # CV criteria saved in odd-numbered elements
      #   coefficients in even-numbered elements
      result <- do.call(rbind, selection[is %% 2 == 1])
      coefs <- do.call(list, selection[is %% 2 == 0])
    } else {
      result <- do.call(rbind, selection)
      coefs <- NULL
    }
    stopCluster(cl)
  } else {
    result <- matrix(0, k, 2L)
    coefs <- vector(k, mode="list")
    for (i in 1L:k){
      selection <- procedure(data, indices[starts[i]:ends[i]],
                             save.coef=save.coef, ...)
      result[i, ] <- selection[[1]]
      if (save.coef) coefs[[i]] <- selection[[2]]
    }
  }
  cv <- weighted.mean(result[, 1L], folds)
  cv.full <- procedure(data, ...)
  adj.cv <- cv + cv.full - weighted.mean(result[, 2L], folds)
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 coefs = if (save.coef) coefs else NULL)
  class(result) <- c("cvSelect", "cv")
  if (reps == 1) {
    return(result)
  } else {
    res <- cvSelect(procedure=procedure, data=data, k=k,
                    reps = reps - 1, save.coef = save.coef,
                    ncores=ncores, ...)
    if (reps  > 2){
      res[[length(res) + 1]] <- result
    } else {
      res <- list(res, result)
    }
    class(res) <- c("cvSelectList", "cvList")
    return(res)
  }
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
#'          k.=log(nrow(Auto)), k=5, reps=3) # via BIC
#' @export
selectStepAIC <- function(data, indices,
                          model, criterion=mse, k.=2,
                          save.coef=TRUE, ...){
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
  list(criterion=c(criterion(y[indices], fit.i),
         criterion(y, fit.o.i)),
       if (save.coef) coefs=coef(model.i) else NULL)
}

#' @describeIn cvSelect print the coefficients from the selected models
#' for the several folds.
#' @param object an object of class \code{"cvSelect"}.
#' @param digits significant digits for printing coefficients,
#' (default \code{3}).
#' @export
compareFolds <- function(object, digits=3, ...){
  UseMethod("compareFolds")
}

#' @export
compareFolds.cvSelect <- function(object, digits=3, ...){
  coefficients <- object$coefs
  if (is.null(coefficients))
    stop("coefficients for folds not available")
  names <- unlist(lapply(coefficients, names))
  counts <- table(names)
  counts <- sort(counts, decreasing=TRUE)
  table <- matrix(NA, length(coefficients), length(counts))
  colnames(table) <- names(counts)
  rownames(table) <- paste("Fold", seq(along=coefficients))
  for (i in seq(along=coefficients)){
    table[i, names(coefficients[[i]])] <- coefficients[[i]]
  }
  printCoefmat(table, na.print="", digits=digits)
}


# selectAllSubsets <- function(data, indices, formula, ...){
#   if (missing(indices)){
#     regsubsets(formula, data=data, nbest=1)
#   }
# }
