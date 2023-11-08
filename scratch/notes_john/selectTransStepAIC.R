

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

    # to get transformation estimates for full data:
    all <- missing(indices)
    if (all) indices <- nrow(data) + 1 # will omit nothing

    # find tranformations of predictors and/or response:

    trans <- if (!missing(predictors) && !missing(response)) {
      selectTrans(
        data = data,
        indices = indices,
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
        indices = indices,
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
        indices = indices,
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
          data[, predictor] <- bcnPower(data[, predictor],
                                        lambda = lambda, gamma = gamma)
        }
      } else {
        for (predictor in predictors) {
          lambda <- powers[paste0("lam.", predictor)]
          data[, predictor] <-
            do.call(family, list(U = data[, predictor],
                                 lambda = lambda))
        }
      }
    }

    if (!missing(response)) {
      if (family.y == "bcnPower") {
        lambda <- powers["lambda"]
        gamma <- powers["gamma"]
        data[, predictor] <- bcnPower(data[, response],
                                      lambda = lambda, gamma = gamma)
      } else {
        lambda <- powers["lambda"]
        data[, response] <- do.call(family.y, list(U = data[, response],
                                                   lambda = lambda))
      }
    }

    # re-estimate with transformed data:
    model <- update(model, data = data)

    # perform variable selection:

    sel <- selectStepAIC(
      data = data,
      indices = indices,
      model = model,
      criterion = criterion,
      AIC = AIC,
      save.coef = save.coef,
      ...
    )
    if (save.coef) sel$coefficients <- c(powers, coef(sel))
    if (all) # transformation & selection on full data set
      return(sel$criterion[2])
    else
      return(sel)
  }
