# version in vignette
cv.multinom <- function (model, data, criterion=BayesRuleMulti, k, reps,
                         seed, ...){
  NextMethod(type="class", criterion=criterion,
             criterion.name=deparse(substitute(criterion)))
}

# new version
cv.multinom <- function(model,
                        data=insight::get_data(model),
                        criterion=BayesRuleMulti,
                        k=10,
                        reps=1,
                        seed=NULL,
                        details = k <= 10,
                        confint = n >= 400,
                        level=0.95,
                        ncores=1,
                        start=FALSE,
                        ...){

  f <- function(i){
    # helper function to compute to compute fitted values,
    #  etc., for each fold i

    indices.i <- fold(folds, i)
    model.i <- if (start) {
      update(model, data=data[ - indices.i, ], start=b, trace=FALSE)
    } else {
      update(model, data=data[ - indices.i, ], trace=FALSE)
    }
    fit.all.i <- predict(model.i, newdata=data, type="class")
    fit.i <- fit.all.i[indices.i]
    # returns:
    #  fit.i: fitted values for the i-th fold
    #  crit.all.i: CV criterion for all cases based on model with
    #              i-th fold omitted
    #  coef.i: coefficients for the model with i-th fold omitted
    list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
         coef.i=coef(model.i))
  }

  fPara <- function(i, multinom, ...){
    # helper function for parallel computation
    #   argument multinom makes multinom() locally available
    #   ... is necessary but not used
    indices.i <- fold(folds, i)
    model.i <- if (start) {
      update(model, data=data[ - indices.i, ], start=b, trace=FALSE)
    } else {
      update(model, data=data[ - indices.i, ], trace=FALSE)
    }
    fit.all.i <- predict(model.i, newdata=data, type="class")
    fit.i <- fit.all.i[indices.i]
    list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
         coef.i=coef(model.i))
  }

  n <- nrow(data)

  # see ?cvCompute for definitions of arguments
  cvCompute(model=model,
            data=data,
            criterion=criterion,
            criterion.name=deparse(substitute(criterion)),
            k=k,
            reps=reps,
            seed=seed,
            details=details,
            confint=confint,
            level=level,
            ncores=ncores,
            type="class",
            start=start,
            f=f,
            fPara=fPara,
            multinom=nnet::multinom)
}
