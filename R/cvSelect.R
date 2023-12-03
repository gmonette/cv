#' Cross-Validate a Model-Selection Procedure
#'
#' \code{cvSelect()} is a general function to cross-validate a model-selection procedure;
#' \code{selectStepAIC()} is a procedure that applies the \code{\link[MASS]{stepAIC}()}
#' model-selection function in the \pkg{MASS} package; \code{selectTrans()} is a procedure
#' for selecting predictor and response transformations in regression, which
#' uses the \code{\link[car]{powerTransform}()} function in the
#' \pkg{car} package; and \code{selectTransAndStepAIC()} combines predictor and response
#' transformation with predictor selection.
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
#' @param ... for \code{cvSelect()}, arguments to be passed to \code{procedure()};
#' for \code{selectStepAIC()}, arguments to be passed to \code{stepAIC()}.
#' @importFrom MASS stepAIC
#' @describeIn cvSelect apply cross-validation to a model-selection procedure.
#' @returns An object of class \code{"cvSelect"}, inheriting from class \code{"cv"}, with the averaged CV criterion
#' (\code{"CV crit"}), the adjusted average CV criterion (\code{"adj CV crit"}),
#' the criterion for the model applied to the full data (\code{"full crit"}),
#' the number of folds (\code{"k"}), the seed for R's random-number
#' generator (\code{"seed"}), and (optionally) a list of coefficients
#' (or, in the case of \code{selectTrans()}, estimated transformation
#' parameters, and in the case of \code{selectTransAndStepAIC()}, both regression coefficients
#' and transformation parameters) for the selected models
#' for each fold (\code{"coefficients"}).
#' If \code{reps} > \code{1}, then an object of class \code{c("cvSelectList", "cvList")} is returned,
#' which is literally a list of \code{c("cvSelect", "cv")} objects.
#' @details
#' The model-selection function supplied as the \code{procedure} argument
#' to \code{cvSelect()} should accept the following arguments:
#' \describe{
#'  \item{\code{data}}{set to the \code{data} argument to \code{cvSelect()}.}
#'  \item{\code{indices}}{the indices of the rows of \code{data} defining the current fold; if missing,
#'  the model-selection procedure is applied to the full \code{data}.}
#'   \item{other arguments}{to be passed via \code{...}
#'   from \code{cvSelect()}.}
#' }
#' \code{procedure()} should return a two-element vector with the result
#' of applying a cross-validation criterion to the cases in
#' the current fold for the model deleting that fold, and to
#' all of the cases, again for the model deleting the current fold.
#'
#' When the \code{indices} argument is missing, \code{procedure()} returns the cross-validation criterion for all of the cases based on
#' the model fit to all of the cases.
#'
#' For examples of model-selection functions for the \code{procedure}
#' argument, see the code for \code{selectStepAIC()},
#' \code{selectTrans()}, and \code{selectTransAndStepAIC()}.
#'
#' For additional information, see the "Cross-validation of regression
#' models" vignette (\code{vignette("cv", package="cv")})
#' and the "Extending the cv package" vignette
#' (\code{vignette("cv-extend", package="cv")}).
#'
#' @seealso \code{\link[MASS]{stepAIC}}, \code{\link[car]{bcPower}},
#' \code{\link[car]{powerTransform}}
#'
#' @export
cvSelect <- function(procedure, data, criterion=mse,
                     model, response,
                     k=10, reps=1,
                     save.coef = k <= 10,
                     seed, ncores=1, ...){
  n <- nrow(data)
  y <- if (!missing(model)) getResponse(model) else data[[response]]
  if (missing(model)) model <- NULL
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
    if (k == n){
      if (reps > 1) stop("reps should not be > 1 for n-fold CV")
      if (!missing(seed) && !is.null(seed)) message("Note: seed ignored for n-fold CV")
      seed <- NULL
    }
    seed <- NULL
  }
  nk <-  n %/% k # number of cases in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  indices <- if (n > k) sample(n, n)  else 1:n # permute cases
  yhat <- numeric(n)
  crit.all.i <- numeric(k)

  if (ncores > 1){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    arglist <- c(list(data=data, indices=1, save.coef=save.coef,
                      criterion=criterion, model=model),
                 list(...))
    selection <- foreach(i = 1L:k) %dopar% {
      # the following deals with a scoping issue that can
      #   occur with args passed via ...
      arglist$indices <- indices[starts[i]:ends[i]]
      do.call(procedure, arglist)
    }
    stopCluster(cl)
    for (i in 1L:k){
      yhat[indices[starts[i]:ends[i]]] <- selection[[i]]$fit.i
      crit.all.i[i] <- selection[[i]]$crit.all.i
    }
    coefs <- if (save.coef){
      lapply(selection, function(x) x$coefficients)
    } else {
      NULL
    }
  } else {
    coefs <- vector(k, mode="list")
    for (i in 1L:k){
      indices.i <- indices[starts[i]:ends[i]]
      selection <- procedure(data, indices.i,
                             save.coef=save.coef,
                             criterion=criterion, model=model, ...)
      crit.all.i[i] <- selection$crit.all.i
      yhat[indices.i] <- selection$fit.i
      if (save.coef) coefs[[i]] <- selection$coefficients
    }
  }
  cv <- criterion(y, yhat)
  cv.full <- procedure(data, model=model, criterion=criterion, ...)
  adj.cv <- cv + cv.full - weighted.mean(crit.all.i, folds)
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 coefficients = if (save.coef) coefs else NULL)
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
#' @param criterion a CV criterion ("cost" or lack-of-fit) function.
#' @param AIC if \code{TRUE} (the default) use the AIC as the
#' model-selection criterion; if \code{FALSE}, use the BIC.
#' The \code{k} argument to \code{\link[MASS]{stepAIC}()}
#' is set accordingly (note that this is distinct from the number of
#' folds \code{k}).
#' @examples
#' data("Auto", package="ISLR2")
#' m.auto <- lm(mpg ~ . - name - origin, data=Auto)
#' cvSelect(selectStepAIC, Auto, seed=123, model=m.auto)
#' cvSelect(selectStepAIC, Auto, seed=123, model=m.auto,
#'          AIC=FALSE, k=5, reps=3) # via BIC
#' @export
selectStepAIC <- function(data, indices,
                          model, criterion=mse, AIC=TRUE,
                          save.coef=TRUE, ...){
  y <- getResponse(model)
  if (missing(indices)) {
    k. = if (AIC) 2 else log(nrow(data))
    model.i <- MASS::stepAIC(model, trace=FALSE, k=k., ...)
    fit.all.i <- predict(model.i, newdata=data, type="response")
    return(criterion(y, fit.all.i))
  }
  k. <- if (AIC) 2 else log(nrow(data) - length(indices))
  model <- update(model, data=data[-indices, ])
  model.i <- MASS::stepAIC(model, trace=FALSE, k=k., ...)
  fit.all.i <- predict(model.i, newdata=data, type="response")
  fit.i <- fit.all.i[indices]
  list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
       coefficients=if (save.coef) coef(model.i) else NULL)
}


transX <- function(model, predictors, family="bcPower",
                   rounded=TRUE, data=insight::get_data(model)){
  # Find transformations of predictors towards normality
  # Returns: a list of the selected transformations, with $lamdas and
  #          $gammas components (the latter may be NULL)
  # Args:
  #   model: e.g., an "lm" model object
  #   predictors: names of predictors to transform
  #   rounded: use "rounded" transformations, see ?car::powerTransform
  #   family: transformation family recognized by car::powerTransform()
  #   data: data to which the model was fit
  trans <- car::powerTransform(data[, predictors], family=family)
  lambdas <- if(rounded) trans$roundlam else trans$lambda
  names(lambdas) <- paste0("lam.", predictors)
  gammas <- trans$gamma
  if (!is.null(gammas)){
    names(gammas) <- paste0("gam.", predictors)
  }
  list(lambdas=lambdas, gammas=gammas) # gammas may be NULL
}

transy <-  function(model, family="bcPower", rounded=TRUE){
  trans <- car::powerTransform(model, family=family)
  lambda <- if(rounded) trans$roundlam else trans$lambda
  gamma <- trans$gamma
  c(lambda=as.vector(lambda), gamma=as.vector(gamma))
}

bcPowerInverse <- function (y, lambda){
  if (abs(lambda) < sqrt(.Machine$double.eps)) {
    exp(y)
  } else {
    (y*lambda + 1)^(1/lambda)
  }
}

basicPowerInverse <- function (y, lambda){
  if (abs(lambda) < sqrt(.Machine$double.eps)) {
    exp(y)
  } else {
    y^(1/lambda)
  }
}

yjPowerInverse <- function(y, lambda) {
  neg <- y < 0
  y[!neg] <- if(abs(lambda) < sqrt(.Machine$double.eps)) {
    exp(y[!neg]) - 1
  } else {
    (y[!neg]*lambda + 1)^(1/lambda) - 1
  }
  y[neg] <- if(abs(lambda - 2) < sqrt(.Machine$double.eps)) {
    -expm1(-y[neg])
  } else {
    1 - (-(2 - lambda)*y[neg] + 1)^(1/(2 - lambda))
  }
  y
}

#' @describeIn cvSelect select transformations of the predictors and response.
#' @param predictors character vector of names of the predictors in the model
#' to transform; if missing, no predictors will be transformed.
#' @param response name of the response variable; if missing, the response
#' won't be transformed.
#' @param family transformation family for the predictors, one of
#' \code{"bcPower", "bcnPower", "yjPower", "basicPower"},
#' with \code{"bcPower"} as the default. These are the names of transformation
#' functions in the \pkg{car} package; see \code{\link[car]{bcPower}()}
#' @param family.y transformation family for the response,
#' with \code{"bcPower"} as the default.
#' @param rounded if \code{TRUE} (the default) use nicely rounded versions
#' of the estimated transformation parameters (see \code{\link[car]{bcPower}()}).
#' @examples
#'
#' data("Prestige", package="carData")
#' m.pres <- lm(prestige ~ income + education + women,
#'              data=Prestige)
#' cvt <- cvSelect(selectTrans, data=Prestige, model=m.pres, seed=123,
#'                 predictors=c("income", "education", "women"),
#'                 response="prestige", family="yjPower")
#' cvt
#' compareFolds(cvt)
#' coef(cvt, average=median, NAs=1) # NAs not really needed here
#' cv(m.pres, seed=123)
#' @export
selectTrans <- function(data, indices, save.coef=TRUE, model,
                        criterion=mse, predictors, response,
                        family=c("bcPower", "bcnPower", "yjPower", "basicPower"),
                        family.y=c("bcPower", "bcnPower", "yjPower", "basicPower"),
                        rounded=TRUE,
                        ...){
  if (missing(predictors) && missing(response))
    stop("'predictors' and 'response' arguments both missing;",
         "\n no transformations specified")
  y <- getResponse(model)
  family <- match.arg(family)
  family.y <- match.arg(family.y)
  powertrans <- switch(family,
                       bcPower = car::bcPower,
                       bcnPower = car::bcnPower,
                       yjPower = car:: yjPower,
                       basicPower = car::basicPower
  )

  powertrans.y <- switch(family.y,
                         bcPower = car::bcPower,
                         bcnPower = car::bcnPower,
                         yjPower = car:: yjPower,
                         basicPower = car::basicPower
  )

  inverse.pt.y <- switch(family.y,
                         bcPower = bcPowerInverse,
                         bcnPower = car::bcnPowerInverse,
                         yjPower = yjPowerInverse,
                         basicPower = basicPowerInverse
  )

  if (missing(indices)) {
    indices <- nrow(data) + 1 # will use full sample
    full.sample <- TRUE
  } else {
    # refit model omitting current fold
    model <- update(model, data = data[-indices, ])
    full.sample <- FALSE
  }

  # find predictor transformations:
  if (!missing(predictors)){
    trans <- transX(model, predictors, family, rounded)
    lambdas <- trans$lambdas
    gammas <- trans$gammas # nb: there may be no gammas
    # transform predictors:
    for (pred in predictors){
      data[, pred] <- if (is.null(gammas)){
        powertrans(data[, pred], lambdas[paste0("lam.", pred)])
      } else {
        powertrans(data[, pred],
                   lambda=lambdas[paste0("lam.", pred)],
                   gamma=gammas[paste0("gam.", pred)])
      }
    }
    model <- update(model, data=data[-indices, ])
  } else {
    lambdas <- gammas <- NULL
  }

  # transform response:
  if (!missing(response)){
    transy <- transy(model, family.y)
    data[, response] <- if (is.na(transy["gamma"])){
      powertrans.y(data[, response], lambda=transy["lambda"])
    } else {
      powertrans.y(data[, response],
                   lambda=transy["lambda"],
                   gamma=transy["gamma"])
    }
    # refit model with transformed predictors and response,
    #   omitting current fold
    model <- update(model, data=data[-indices, ])
  } else {
    transy <- NULL
  }

  # get predicted values for *all* cases:

  fit.o.i <- if (!missing(response)){
    if (is.na(transy["gamma"])){
      inverse.pt.y(predict(model, newdata = data),
                   lambda=transy["lambda"])
    } else {
      inverse.pt.y(predict(model, newdata = data),
                   lambda=transy["lambda"],
                   gamma=transy["gamma"])
    }
  } else {
    predict(model, newdata = data)
  }

  if (full.sample) return(criterion(y, fit.o.i))
  # ... and for current fold only:
  fit.i <- fit.o.i[indices]
  # compute and return CV criteria and transformation parameters:
  list(criterion = c(criterion(y[indices], fit.i),
                     criterion(y, fit.o.i)),
       coefficients = if (save.coef) c(lambdas, gammas, transy)
  )
}


#' @describeIn cvSelect select transformations of the predictors and response,
#' and then select predictors.
#' @examples
#'
#' Auto$year <- as.factor(Auto$year)
#' Auto$origin <- factor(Auto$origin,
#'                       labels=c("America", "Europe", "Japan"))
#' rownames(Auto) <- make.names(Auto$name, unique=TRUE)
#' Auto$name <- NULL
#' m.auto <- lm(mpg ~ . , data=Auto)
#' cvs <- cvSelect(selectTransStepAIC, data=Auto, seed=76692, model=m.auto,
#'                 criterion=medAbsErr,
#'                 predictors=c("cylinders", "displacement", "horsepower",
#'                              "weight", "acceleration"),
#'                 response="mpg", AIC=FALSE)
#' cvs
#' compareFolds(cvs)
#' @export
selectTransStepAIC <- function(data,
                               indices,
                               save.coef = TRUE,
                               model,
                               criterion = mse,
                               predictors,
                               response,
                               family = c("bcPower", "bcnPower", "yjPower", "basicPower"),
                               family.y = c("bcPower", "bcnPower", "yjPower", "basicPower"),
                               rounded = TRUE,
                               AIC = TRUE,
                               ...) {

  family <- match.arg(family)
  family.y <- match.arg(family.y)

  powertrans <- switch(family,
                       bcPower = car::bcPower,
                       bcnPower = car::bcnPower,
                       yjPower = car:: yjPower,
                       basicPower = car::basicPower
  )

  powertrans.y <- switch(family.y,
                         bcPower = car::bcPower,
                         bcnPower = car::bcnPower,
                         yjPower = car:: yjPower,
                         basicPower = car::basicPower
  )

  inverse.pt.y <- switch(family.y,
                         bcPower = bcPowerInverse,
                         bcnPower = car::bcnPowerInverse,
                         yjPower = yjPowerInverse,
                         basicPower = basicPowerInverse
  )

  y <- getResponse(model) # untransformed response

  # find tranformations of predictors and/or response:

  # if indices is missing, use full data set;
  # otherwise remove cases in current fold (indices)
  inds <- if (missing(indices)) length(y) + 1 else indices

  trans <- if (!missing(predictors) && !missing(response)) {
    selectTrans(
      data = data,
      indices = inds,
      save.coef = TRUE,
      model = model,
      criterion = criterion,
      predictors = predictors,
      response = response,
      family = family,
      family.y = family.y,
      rounded = rounded,
      ...
    )

  } else if (!missing(predictors)) {
    selectTrans(
      data = data,
      indices = inds,
      save.coef = TRUE,
      model = model,
      criterion = criterion,
      predictors = predictors,
      family = family,
      family.y = family.y,
      rounded = rounded,
      ...
    )

  } else if (!missing(response)) {
    selectTrans(
      data = data,
      indices = inds,
      save.coef = TRUE,
      model = model,
      criterion = criterion,
      response = response,
      family = family,
      family.y = family.y,
      rounded = rounded,
      ...
    )

  } else {
    stop(
      "'predictors' and 'response' arguments both missing;",
      "\n no transformations specified"
    )
  }

  # apply transformations to the data:

  powers <- coef(trans)
  if (!missing(predictors)) {
    if (family == "bcnPower") {
      for (predictor in predictors) {
        lambda <- powers[paste0("lam.", predictor)]
        gamma <- powers[paste0("gam.", predictor)]
        data[, predictor] <- car::bcnPower(data[, predictor],
                                           lambda = lambda, gamma = gamma)
      }
    } else {
      for (predictor in predictors) {
        lambda <- powers[paste0("lam.", predictor)]
        data[, predictor] <- do.call(powertrans, list(U = data[, predictor],
                                                      lambda = lambda))
      }
    }
  }

  if (!missing(response)) {
    if (family.y == "bcnPower") {
      lambda <- powers["lambda"]
      gamma <- powers["gamma"]
      data[, response] <- car::bcnPower(data[, response],
                                        lambda = lambda, gamma = gamma)
    } else {
      lambda <- powers["lambda"]
      data[, response] <- do.call(powertrans.y, list(U = data[, response],
                                                     lambda = lambda))
    }
  }

  # re-estimate with transformed data:
  model <- update(model, data = data)

  # perform variable selection:

  if (missing(indices)) {
    k. = if (AIC) 2 else log(nrow(data))
    model.i <- MASS::stepAIC(model, trace=FALSE, k=k., ...)
    fit.all.i <- predict(model.i, newdata=data, type="response")
  } else {
    k. <- if (AIC) 2 else log(nrow(data) - length(indices))
    model <- update(model, data=data[-indices, ])
    model.i <- MASS::stepAIC(model, trace=FALSE, k=k., ...)
    fit.all.i <- predict(model.i, newdata=data, type="response")
  }

  # express fitted values on original response scale
  if (!missing(response)){
    if (is.na(powers["gamma"])){
      fit.all.i  <- inverse.pt.y(fit.all.i,
                                 lambda=powers["lambda"])
    } else {
      fit.all.i  <- inverse.pt.y(fit.all.i,
                                 lambda=powers["lambda"],
                                 gamma=powers["gamma"])
    }
  }

  if (missing(indices)) return(criterion(y, fit.all.i))

  list(criterion=c(criterion(y[indices], fit.all.i[indices]),
                   criterion(y, fit.all.i)),
       coefficients=if (save.coef) c(powers, coef(model.i)) else NULL)
}


#' @describeIn cvSelect print the coefficients from the selected models
#' for the several folds.
#' @param object an object of class \code{"cvSelect"}.
#' @param digits significant digits for printing coefficients
#' (default \code{3}).
#' @export
compareFolds <- function(object, digits=3, ...){
  UseMethod("compareFolds")
}

#' @export
compareFolds.cvSelect <- function(object, digits=3, ...){
  coefficients <- object$coefficients
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

#' @describeIn cvSelect extract the coefficients from the selected models
#' for the several folds and possibly average them.
#' @param average if supplied, a function, such as \code{mean} or \code{median},
#' to use us in averaging estimates across folds; if missing, the
#' estimates for each fold are returned.
#' @param NAs values to substitute for \code{NA}s in calculating
#' averaged estimates; the default, \code{0}, is appropriate, e.g.,
#' for regression coefficients; the value \code{1} might be appropriate
#' for power-transformation estimates.
#' @importFrom utils capture.output
#' @export
coef.cvSelect <- function(object, average, NAs=0, ...){
  capture.output(coef <- compareFolds(object))
  if (missing(average)){
    return(coef)
  }
  coef[is.na(coef)] <- NAs
  apply(coef, 2, average)
}
