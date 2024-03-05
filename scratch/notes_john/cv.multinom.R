# version in vignette
cv.multinom <- function (model, data, criterion=BayesRuleMulti, k, reps,
                         seed, ...){
  NextMethod(type="class", criterion=criterion,
             criterion.name=deparse(substitute(criterion)),
             reg.fn=nnet::multinom, reg.fn.name="multinom")
}

# new version
# cv.multinom <- function(model, data=insight::get_data(model),
#          criterion=BayesRuleMulti, k=10, reps=1, seed=NULL,
#          details = k <= 10,
#          confint = n >= 400, level=0.95,
#          ncores=1,
#          start=FALSE,  ...){
#
#   f <- function(i){
#     # helper function to compute cv criterion for each fold
#     indices.i <- fold(folds, i)
#     model.i <- if (start) {
#       update(model, data=data[ - indices.i, ], start=b, trace=FALSE)
#     } else {
#       update(model, data=data[ - indices.i, ], trace=FALSE)
#     }
#     fit.all.i <- predict(model.i, newdata=data, type="class")
#     fit.i <- fit.all.i[indices.i]
#     list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
#          coef.i=coef(model.i))
#   }
#
#   fPara <- function(i, multinom, ...){
#     # helper function to compute cv criterion for each fold
#     indices.i <- fold(folds, i)
#     model.i <- if (start) {
#       update(model, data=data[ - indices.i, ], start=b, trace=FALSE)
#     } else {
#       update(model, data=data[ - indices.i, ], trace=FALSE)
#     }
#     fit.all.i <- predict(model.i, newdata=data, type="class")
#     fit.i <- fit.all.i[indices.i]
#     list(fit.i=fit.i, crit.all.i=criterion(y, fit.all.i),
#          coef.i=coef(model.i))
#   }
#
#   n <- nrow(data)
#
#   cvCompute(model=model, data=data, criterion=criterion,
#             criterion.name=deparse(substitute(criterion)),
#             k=k, reps=reps, seed=seed, details=details, confint=confint,
#             level=level, ncores=ncores, type="class", start=start,
#             f=f, fPara=fPara, multinom=nnet::multinom)
# }
