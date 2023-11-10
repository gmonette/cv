

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
                               response.scale = c("original", "transformed"),
                               AIC = TRUE,
                               ...) {

  family <- match.arg(family)
  family.y <- match.arg(family.y)
  response.scale <- match.arg(response.scale)

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
                         bcPower = cv:::bcPowerInverse,
                         bcnPower = car::bcnPowerInverse,
                         yjPower = cv:::yjPowerInverse,
                         basicPower = cv:::basicPowerInverse
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
      data[, predictor] <- car::bcnPower(data[, response],
                                         lambda = lambda, gamma = gamma)
    } else {
      lambda <- powers["lambda"]
      data[, response] <- do.call(powertrans.y, list(U = data[, response],
                                                     lambda = lambda))
    }
  }

  # re-estimate with transformed data:
  model <- update(model, data = data)
  y.trans <- if (response.scale == "transformed")
    getResponse(model) else NULL

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
  fit.all.i  <- if (!missing(response) && response.scale == "original"){
    if (is.na(powers["gamma"])){
      inverse.pt.y(predict(model.i, newdata = data),
                   lambda=powers["lambda"])
    } else {
      inverse.pt.y(predict(model.i, newdata = data),
                   lambda=powers["lambda"],
                   gamma=powers["gamma"])
    }
  } else {
    predict(model.i, newdata = data)
  }

  if (response.scale == "transformed") y <- y.trans

  if (missing(indices)) return(criterion(y, fit.all.i))

  list(criterion=c(criterion(y[indices], fit.all.i[indices]),
                   criterion(y, fit.all.i)),
       coefficients=if (save.coef) c(powers, coef(model.i)) else NULL)
}
