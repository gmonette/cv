#' Cross-Validate a Model-Selection Procedure
#'
#' The \code{cv()} \code{"function"} method
#' is a general function to cross-validate a model-selection procedure,
#' such as the following:
#' \code{selectStepAIC()} is a procedure that applies the \code{\link[MASS]{stepAIC}()}
#' model-selection function in the \pkg{MASS} package; \code{selectTrans()} is a procedure
#' for selecting predictor and response transformations in regression, which
#' uses the \code{\link[car]{powerTransform}()} function in the
#' \pkg{car} package; \code{selectTransAndStepAIC()} combines predictor and response
#' transformations with predictor selection; and \code{selectModelList()}
#' uses cross-validation to select a model from a list of models created by
#' \code{\link{models}()} and employs (recursive) cross-validation to assess the predictive
#' accuracy of this procedure.
#'
#' @param data full data frame for model selection.
#' @param y.expression normally the response variable is found from the
#' \code{model} or \code{working.model} argument; but if, for a particular selection procedure, the
#' \code{model} or \code{working.model} argument is absent, or if the response can't be inferred from the
#' model, the response can be specified by an expression, such as \code{expression(log(income))},
#' to be evaluated within the data set provided by the \code{data} argument.
#' @param k perform k-fold cross-validation (default is 10); \code{k}
#' may be a number or \code{"loo"} or \code{"n"} for n-fold (leave-one-out)
#' cross-validation.
#' @param confint if \code{TRUE} (the default if the number of cases is 400
#' or greater), compute a confidence interval for the bias-corrected CV
#' criterion, if the criterion is the average of casewise components.
#' @param level confidence level (default \code{0.95}).
#' @param details if \code{TRUE}, save detailed information about the value of the
#' CV criterion for the cases in each fold and the regression coefficients
#' (and possibly other information)
#' with that fold deleted; default is \code{TRUE} if \code{k} is 10 or smaller,
#' \code{FALSE} otherwise.
#' @param save.model save the model that's selected using the \emph{full} data set
#' (default, \code{FALSE}).
#' @param reps number of times to replicate k-fold CV (default is \code{1})
#' @param seed for R's random number generator; not used for n-fold cross-validation.
#' If not explicitly set, a seed is randomly generated and saved to make the results
#' reproducible. In some cases, for internal use only, \code{seed} is set to
#' \code{FALSE} to suppress automatically setting the seed.
#' @param ncores number of cores to use for parallel computations
#'        (default is \code{1}, i.e., computations aren't done in parallel)
#' @param ... for \code{cvSelect()} and the \code{cv()} \code{"function"} method,
#' arguments to be passed to \code{procedure()};
#' for \code{selectStepAIC()} and \code{selectTransStepAIC()},
#' arguments to be passed to \code{stepAIC()}.
#' @importFrom MASS stepAIC
#' @returns An object of class \code{"cvSelect"},
#' inheriting from class \code{"cv"}, with the CV criterion
#' (\code{"CV crit"}), the bias-adjusted CV criterion (\code{"adj CV crit"}),
#' the criterion for the model applied to the full data (\code{"full crit"}),
#' the confidence interval and level for the bias-adjusted CV criterion (\code{"confint"}),
#' the number of folds (\code{"k"}), the seed for R's random-number
#' generator (\code{"seed"}), and (optionally) a list of coefficients
#' (or, in the case of \code{selectTrans()}, estimated transformation
#' parameters, and in the case of \code{selectTransAndStepAIC()}, both regression coefficients
#' and transformation parameters) for the selected models
#' for each fold (\code{"coefficients"}).
#' If \code{reps} > \code{1}, then an object of class \code{c("cvSelectList", "cvList")} is returned,
#' which is literally a list of \code{c("cvSelect", "cv")} objects.
#' @details
#' The model-selection function supplied as the \code{procedure} (for \code{cvSelect()})
#' or \code{model} (for \code{cv()}) argument
#' should accept the following arguments:
#' \describe{
#'  \item{\code{data}}{set to the \code{data} argument to \code{cvSelect()} or \code{cv()}.}
#'  \item{\code{indices}}{the indices of the rows of \code{data} defining the current fold; if missing,
#'  the model-selection procedure is applied to the full \code{data}.}
#'   \item{other arguments}{to be passed via \code{...}
#'   from \code{cvSelect()} or \code{cv()}.}
#' }
#' \code{procedure()} or \code{model()} should return a list with the following
#' named elements: \code{fit.i}, the vector of predicted values for the cases in
#' the current fold computed from the model omitting these cases;
#' \code{crit.all.i}, the CV criterion computed for all of the cases using
#' the model omitting the current fold; and (optionally) \code{coefficients},
#' parameter estimates from the model computed omitting the current fold.
#'
#' When the \code{indices} argument is missing, \code{procedure()} returns the cross-validation criterion for all of the cases based on
#' the model fit to all of the cases.
#'
#' For examples of model-selection functions for the \code{procedure}
#' argument, see the code for \code{selectStepAIC()},
#' \code{selectTrans()}, and \code{selectTransAndStepAIC()}.
#'
#' For additional information, see the "Cross-validating model selection"
#' vignette (\code{vignette("cv-select", package="cv")})
#' and the "Extending the cv package" vignette
#' (\code{vignette("cv-extend", package="cv")}).
#'
#' @seealso \code{\link[MASS]{stepAIC}}, \code{\link[car]{bcPower}},
#' \code{\link[car]{powerTransform}}, \code{\link{cv}}.
#'
#' @param indices indices of cases in data defining the current fold.
#' @param model a regression model object fit to data, or for the
#' \code{cv()} \code{"function"} method, a model-selection procedure function
#' (see Details).
#' @param working.model a regression model object fit to data, typically
#' to begin a model-selection process; for use with \code{selectModelList()},
#' a list of competing models created by \code{\link{models}()}.
#' @param criterion a CV criterion ("cost" or lack-of-fit) function.
#' @param AIC if \code{TRUE} (the default) use the AIC as the
#' model-selection criterion; if \code{FALSE}, use the BIC.
#' The \code{k} argument to \code{\link[MASS]{stepAIC}()}
#' is set accordingly (note that this is distinct from the number of
#' folds \code{k}).
#' @examples
#' if (requireNamespace("ISLR2", quietly=TRUE)){
#' withAutoprint({
#' data("Auto", package="ISLR2")
#' m.auto <- lm(mpg ~ . - name - origin, data=Auto)
#' cv(selectStepAIC, Auto, seed=123, working.model=m.auto)
#' cv(selectStepAIC, Auto, seed=123, working.model=m.auto,
#'          AIC=FALSE, k=5, reps=3) # via BIC
#' })
#' } else {
#' cat("\n install the 'ISLR2' package to run these examples\n")
#' }

#' @describeIn cv.function \code{cv()} method for applying a model
#' model-selection (or specification) procedure.
#' @export
cv.function <- function(model,
                        data,
                        criterion = mse,
                        k = 10L,
                        reps = 1L,
                        seed = NULL,
                        working.model = NULL,
                        y.expression = NULL,
                        confint = n >= 400L,
                        level = 0.95,
                        details = k <= 10L,
                        save.model = FALSE,
                        ncores = 1L,
                        ...) {
  n <- nrow(data)
  cvSelect(
    procedure = model,
    data = data,
    criterion = criterion,
    criterion.name = deparse(substitute(criterion)),
    model = working.model,
    y.expression = y.expression,
    k = k,
    confint = confint,
    level = level,
    reps = reps,
    details = details,
    save.model = save.model,
    seed = seed,
    ncores = ncores,
    ...
  )
}

#' @describeIn cv.function select a regression model using the
#' \code{\link[MASS]{stepAIC}()} function in the \pkg{MASS} package.
#' @export
selectStepAIC <- function(data,
                          indices,
                          model,
                          criterion = mse,
                          AIC = TRUE,
                          details = TRUE,
                          save.model = FALSE,
                          ...) {
  y <- GetResponse(model)
  if (missing(indices)) {
    k. = if (AIC)
      2
    else
      log(nrow(data))
    model <- MASS::stepAIC(model, trace = FALSE, k = k., ...)
    yhat <- predict(model, newdata = data, type = "response")
    return(list(
      criterion = criterion(y, yhat),
      model = if (save.model)
        model
      else
        NULL
    ))
  }
  k. <- if (AIC)
    2
  else
    log(nrow(data) - length(indices))
  model <- update(model, data = data[-indices,])
  model.i <- MASS::stepAIC(model, trace = FALSE, k = k., ...)
  fit.all.i <- predict(model.i, newdata = data, type = "response")
  fit.i <- fit.all.i[indices]
  list(
    fit.i = fit.i,
    crit.all.i = criterion(y, fit.all.i),
    coefficients = if (details)
      coef(model.i)
    else
      NULL
  )
}


transX <- function(model,
                   predictors,
                   family = "bcPower",
                   rounded = TRUE,
                   data = insight::get_data(model)) {
  # Find transformations of predictors towards normality
  # Returns: a list of the selected transformations, with $lamdas and
  #          $gammas components (the latter may be NULL)
  # Args:
  #   model: e.g., an "lm" model object
  #   predictors: names of predictors to transform
  #   rounded: use "rounded" transformations, see ?car::powerTransform
  #   family: transformation family recognized by car::powerTransform()
  #   data: data to which the model was fit
  trans <- car::powerTransform(data[, predictors], family = family)
  lambdas <- if (rounded)
    trans$roundlam
  else
    trans$lambda
  names(lambdas) <- paste0("lam.", predictors)
  gammas <- trans$gamma
  if (!is.null(gammas)) {
    names(gammas) <- paste0("gam.", predictors)
  }
  list(lambdas = lambdas, gammas = gammas) # gammas may be NULL
}

transy <-  function(model,
                    family = "bcPower",
                    rounded = TRUE) {
  trans <- car::powerTransform(model, family = family)
  lambda <- if (rounded)
    trans$roundlam
  else
    trans$lambda
  gamma <- trans$gamma
  c(lambda = as.vector(lambda), gamma = as.vector(gamma))
}

bcPowerInverse <- function (y, lambda) {
  if (abs(lambda) < sqrt(.Machine$double.eps)) {
    exp(y)
  } else {
    (y * lambda + 1) ^ (1 / lambda)
  }
}

basicPowerInverse <- function (y, lambda) {
  if (abs(lambda) < sqrt(.Machine$double.eps)) {
    exp(y)
  } else {
    y ^ (1 / lambda)
  }
}

yjPowerInverse <- function(y, lambda) {
  neg <- y < 0
  y[!neg] <- if (abs(lambda) < sqrt(.Machine$double.eps)) {
    exp(y[!neg]) - 1
  } else {
    (y[!neg] * lambda + 1) ^ (1 / lambda) - 1
  }
  y[neg] <- if (abs(lambda - 2) < sqrt(.Machine$double.eps)) {
    -expm1(-y[neg])
  } else {
    1 - (-(2 - lambda) * y[neg] + 1) ^ (1 / (2 - lambda))
  }
  y
}

#' @param predictors character vector of names of the predictors in the model
#' to transform; if missing, no predictors will be transformed.
#' @param response name of the response variable; if missing, the response
#' won't be transformed.
#' @param family transformation family for the predictors, one of
#' \code{"bcPower", "bcnPower", "yjPower", "basicPower"},
#' with \code{"bcPower"} as the default. These are the names of transformation
#' functions in the \pkg{car} package; see \code{\link[car]{bcPower}()}.
#' @param family.y transformation family for the response,
#' with \code{"bcPower"} as the default.
#' @param rounded if \code{TRUE} (the default) use nicely rounded versions
#' of the estimated transformation parameters (see \code{\link[car]{bcPower}()}).
#' @examples
#' if (requireNamespace("carData", quietly=TRUE)){
#' withAutoprint({
#' data("Prestige", package="carData")
#' m.pres <- lm(prestige ~ income + education + women,
#'              data=Prestige)
#' cvt <- cv(selectTrans, data=Prestige, working.model=m.pres, seed=123,
#'           predictors=c("income", "education", "women"),
#'           response="prestige", family="yjPower")
#' cvt
#' compareFolds(cvt)
#' coef(cvt, average=median, NAs=1) # NAs not really needed here
#' cv(m.pres, seed=123)
#' })
#' } else {
#' cat("install the 'carData' package to run these examples\n")
#' }
#' @describeIn cv.function select transformations of the predictors and response
#' using \code{\link[car]{powerTransform}()} in the \pkg{car} package.
#' @export
selectTrans <- function(data,
                        indices,
                        details = TRUE,
                        save.model = FALSE,
                        model,
                        criterion = mse,
                        predictors,
                        response,
                        family = c("bcPower", "bcnPower", "yjPower", "basicPower"),
                        family.y = c("bcPower", "bcnPower", "yjPower", "basicPower"),
                        rounded = TRUE,
                        ...) {
  if (missing(predictors) && missing(response))
    stop(
      "'predictors' and 'response' arguments both missing;",
      "\n no transformations specified"
    )
  y <- GetResponse(model)
  family <- match.arg(family)
  family.y <- match.arg(family.y)
  powertrans <- switch(
    family,
    bcPower = car::bcPower,
    bcnPower = car::bcnPower,
    yjPower = car::yjPower,
    basicPower = car::basicPower
  )

  powertrans.y <- switch(
    family.y,
    bcPower = car::bcPower,
    bcnPower = car::bcnPower,
    yjPower = car::yjPower,
    basicPower = car::basicPower
  )

  inverse.pt.y <- switch(
    family.y,
    bcPower = bcPowerInverse,
    bcnPower = car::bcnPowerInverse,
    yjPower = yjPowerInverse,
    basicPower = basicPowerInverse
  )

  if (missing(indices)) {
    indices <- nrow(data) + 1L # will use full sample
    full.sample <- TRUE
  } else {
    # refit model omitting current fold
    model <- update(model, data = data[-indices,])
    full.sample <- FALSE
  }

  # find predictor transformations:
  if (!missing(predictors)) {
    trans <- transX(model, predictors, family, rounded)
    lambdas <- trans$lambdas
    gammas <- trans$gammas # nb: there may be no gammas
    # transform predictors:
    for (pred in predictors) {
      data[, pred] <- if (is.null(gammas)) {
        powertrans(data[, pred], lambdas[paste0("lam.", pred)])
      } else {
        powertrans(data[, pred],
                   lambda = lambdas[paste0("lam.", pred)],
                   gamma = gammas[paste0("gam.", pred)])
      }
    }
    model <- update(model, data = data[-indices,])
  } else {
    lambdas <- gammas <- NULL
  }

  # transform response:
  if (!missing(response)) {
    transy <- transy(model, family.y)
    data[, response] <- if (is.na(transy["gamma"])) {
      powertrans.y(data[, response], lambda = transy["lambda"])
    } else {
      powertrans.y(data[, response],
                   lambda = transy["lambda"],
                   gamma = transy["gamma"])
    }
    # refit model with transformed predictors and response,
    #   omitting current fold
    model <- update(model, data = data[-indices,])
  } else {
    transy <- NULL
  }

  # get predicted values for *all* cases:

  fit.o.i <- if (!missing(response)) {
    if (is.na(transy["gamma"])) {
      inverse.pt.y(predict(model, newdata = data),
                   lambda = transy["lambda"])
    } else {
      inverse.pt.y(predict(model, newdata = data),
                   lambda = transy["lambda"],
                   gamma = transy["gamma"])
    }
  } else {
    predict(model, newdata = data)
  }

  if (full.sample) {
    return(list(
      criterion = criterion(y, fit.o.i),
      model = if (save.model){
        if (details){
          model$additional.coefficients <- c(lambdas, gammas, transy)
        }
        model
      }
      else
        NULL
    ))
  }
  # ... and for current fold only:
  # compute and return CV info and transformation parameters:
  list(
    fit.i = fit.o.i[indices],
    crit.all.i = criterion(y, fit.o.i),
    coefficients = if (details)
      c(lambdas, gammas, transy)
  )
}


#' @describeIn cv.function select transformations of the predictors and response,
#' and then select predictors.
#' @examples
#' if (requireNamespace("ISLR2", quietly=TRUE)){
#' withAutoprint({
#' Auto$year <- as.factor(Auto$year)
#' Auto$origin <- factor(Auto$origin,
#'                       labels=c("America", "Europe", "Japan"))
#' rownames(Auto) <- make.names(Auto$name, unique=TRUE)
#' Auto$name <- NULL
#' m.auto <- lm(mpg ~ . , data=Auto)
#' cvs <- cv(selectTransStepAIC, data=Auto, seed=76692, working.model=m.auto,
#'           criterion=medAbsErr,
#'           predictors=c("cylinders", "displacement", "horsepower",
#'                        "weight", "acceleration"),
#'           response="mpg", AIC=FALSE)
#' cvs
#' compareFolds(cvs)
#' })
#' }
#' @export
selectTransStepAIC <- function(data,
                               indices,
                               details = TRUE,
                               save.model = FALSE,
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

  powertrans <- switch(
    family,
    bcPower = car::bcPower,
    bcnPower = car::bcnPower,
    yjPower = car::yjPower,
    basicPower = car::basicPower
  )

  powertrans.y <- switch(
    family.y,
    bcPower = car::bcPower,
    bcnPower = car::bcnPower,
    yjPower = car::yjPower,
    basicPower = car::basicPower
  )

  inverse.pt.y <- switch(
    family.y,
    bcPower = bcPowerInverse,
    bcnPower = car::bcnPowerInverse,
    yjPower = yjPowerInverse,
    basicPower = basicPowerInverse
  )

  y <- GetResponse(model) # untransformed response

  # find tranformations of predictors and/or response:

  # if indices is missing, use full data set;
  # otherwise remove cases in current fold (indices)
  inds <- if (missing(indices))
    length(y) + 1L
  else
    indices

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
        data[, predictor] <-
          do.call(powertrans, list(U = data[, predictor],
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
      data[, response] <-
        do.call(powertrans.y, list(U = data[, response],
                                   lambda = lambda))
    }
  }

  # re-estimate with transformed data:
  model <- update(model, data = data)

  # perform variable selection:

  if (missing(indices)) {
    k. = if (AIC)
      2
    else
      log(nrow(data))
    model.i <- MASS::stepAIC(model, trace = FALSE, k = k., ...)
    fit.all.i <- predict(model.i, newdata = data, type = "response")
  } else {
    k. <- if (AIC)
      2
    else
      log(nrow(data) - length(indices))
    model <- update(model, data = data[-indices,])
    model.i <- MASS::stepAIC(model, trace = FALSE, k = k., ...)
    fit.all.i <- predict(model.i, newdata = data, type = "response")
  }

  # express fitted values on original response scale
  if (!missing(response)) {
    if (is.na(powers["gamma"])) {
      fit.all.i  <- inverse.pt.y(fit.all.i,
                                 lambda = powers["lambda"])
    } else {
      fit.all.i  <- inverse.pt.y(fit.all.i,
                                 lambda = powers["lambda"],
                                 gamma = powers["gamma"])
    }
  }

  # for full sample:
  if (missing(indices)) {
    return(list(
      criterion = criterion(y, fit.all.i),
      model = if (save.model){
        model.i$additional.coefficients <- powers
        model.i
      }
      else
        NULL
    ))
  }

  # ... and for current fold only:
  # compute and return CV info, transformation parameters,
  #   and regression coefficients:
  list(
    fit.i = fit.all.i[indices],
    crit.all.i = criterion(y, fit.all.i),
    coefficients = if (details)
      c(powers, coef(model.i))
    else
      NULL
  )
}

#' @describeIn cv.function select a model using (recursive) CV.
#' @param quietly if \code{TRUE} (the default), simple messages (for example about the
#' value to which the random-number generator seed is set), but not warnings or
#' errors, are suppressed.
#' @param k.recurse the number of folds for recursive CV; defaults
#' to the value of \code{k}; may be specified as \code{"loo"} or
#' \code{"n"} as well as an integer.
#' @export
selectModelList <-
  function(data,
           indices,
           model,
           criterion = mse,
           k = 10L,
           k.recurse = k,
           details =  k <= 10L,
           save.model = FALSE,
           seed = FALSE,
           quietly = TRUE,
           ...) {
    if (missing(indices)) {
      if ((!isFALSE(seed)) && is.null(seed)) {
        seed <- sample(1e6, 1L)
        set.seed(seed)
        message("R RNG seed set to ", seed)
      }
      result <-
        cv(
          model,
          data,
          criterion = criterion,
          k = k.recurse,
          details = FALSE,
          quietly = quietly,
          seed = FALSE,
          ...
        )
      cv.min <- which.min(sapply(result, function(x)
        x$"CV crit"))
      return(list(criterion = result[[cv.min]][["CV crit"]],
                  model = if (save.model)
                    model[[cv.min]]
                  else
                    NULL))
    }
    y <- GetResponse(model[[1L]])
    for (i in seq_along(model)) {
      mod <- model[[i]]
      model[[i]] <- update(mod, data = data[-indices,])
    }
    result <-
      cv(
        model,
        data[-indices,],
        criterion = criterion,
        k = k.recurse,
        details = details,
        quietly = quietly,
        seed = FALSE,
        ...
      )
    cv.min <- which.min(sapply(result, function(x)
      x$"CV crit"))
    model.name <- names(result)[cv.min]
    fit.all.i <-
      predict(model[[cv.min]], newdata = data, type = "response")
    fit.i <- fit.all.i[indices]
    list(
      fit.i = fit.i,
      crit.all.i = criterion(y, fit.all.i),
      coefficients = if (details)
        coef(model[[cv.min]])
      else
        NULL,
      model.name = model.name
    )
  }
#' @examples
#' data("Duncan", package="carData")
#' m1 <- lm(prestige ~ income + education, data=Duncan)
#' m2 <- lm(prestige ~ income + education + type, data=Duncan)
#' m3 <- lm(prestige ~ (income + education)*type, data=Duncan)
#' summary(cv.sel <- cv(selectModelList, data=Duncan, seed=5963,
#'                      working.model=models(m1, m2, m3),
#'                      save.model=TRUE)) # recursive CV
#' selectedModel(cv.sel)
#'
#' @describeIn cv.function print the coefficients from the selected models
#' for the several folds.
#' @param object an object of class \code{"cvSelect"}.
#' @param digits significant digits for printing coefficients
#' (default \code{3}).
#' @export
compareFolds <- function(object, digits = 3, ...) {
  UseMethod("compareFolds")
}

#' @export
compareFolds.default <- function(object, digits = 3, ...) {
  coefficients <- object$details$coefficients
  if (is.null(coefficients))
    stop("details for folds not available")
  cat("CV criterion by folds:\n")
  print(object$details$criterion)
  names <- unlist(lapply(coefficients, names))
  counts <- table(names)
  counts <- sort(counts, decreasing = TRUE)
  table <- matrix(NA, length(coefficients), length(counts))
  colnames(table) <- names(counts)
  rownames(table) <- paste("Fold", seq(along = coefficients))
  for (i in seq(along = coefficients)) {
    table[i, names(coefficients[[i]])] <- coefficients[[i]]
  }
  cat("\nCoefficients by folds:\n")
  printCoefmat(table, na.print = "", digits = digits)
}

#' @describeIn cv.function extract the coefficients from the selected models
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
coef.cvSelect <- function(object, average, NAs = 0, ...) {
  capture.output(coef <- compareFolds(object))
  if (missing(average)) {
    return(coef)
  }
  coef[is.na(coef)] <- NAs
  apply(coef, 2, average)
}

#' @describeIn cv.function extract the selected model from a
#' \code{selectModel} object.
#' @export
selectedModel <- function(object, ...){
  UseMethod("selectedModel")
}

#' @rdname cv.function
#' @export
selectedModel.cvSelect <- function(object, ...){
  object$selected
}

