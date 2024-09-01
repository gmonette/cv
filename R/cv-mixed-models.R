#' Cross-Validate Mixed-Effects Model
#'
#' \code{\link{cv}()} methods for mixed-effect models of class \code{"merMod"}, fit
#' by the \code{\link[lme4]{lmer}()} and \code{\link[lme4]{glmer}()} functions
#' in the \pkg{lme4} package; for models of class \code{"lme"}
#' fit by the \code{\link[nlme]{lme}()} function in the \pkg{nlme}
#' package; and for models of class \code{"glmmTMB"} fit by the
#' \code{\link[glmmTMB]{glmmTMB}()} function in the \pkg{glmmTMB} package.
#'
#' @param model a mixed-effects model object for which a \code{cv()} method is available.
#' @param data data frame to which the model was fit (not usually necessary)
#' @param criterion cross-validation ("cost" or lack-of-fit) criterion function of form \code{f(y, yhat)}
#'        where \code{y} is the observed values of the response and
#'        \code{yhat} the predicted values; the default is \code{\link{mse}}
#'        (the mean-squared error).
#' @param k perform k-fold cross-validation; \code{k}
#' may be a number or \code{"loo"} or \code{"n"} for n-fold (leave-one-out)
#' cross-validation; the default is \code{10} if cross-validating individual
#' cases and \code{"loo"} if cross-validating clusters.
#' @param reps number of times to replicate k-fold CV (default is \code{1}),
#' or greater), compute a confidence interval for the bias-corrected CV
#' criterion, if the criterion is the average of casewise components.
#' @param seed for R's random number generator; optional, if not
#' supplied a random seed will be selected and saved; not needed
#' for n-fold cross-validation
#' @param details if \code{TRUE} (the default if the number of
#' folds \code{k <= 10}), save detailed information about the value of the
#' CV criterion for the cases in each fold and the regression coefficients
#' with that fold deleted.
#' @param ncores number of cores to use for parallel computations
#'        (default is \code{1}, i.e., computations aren't done in parallel)
#' @param clusterVariables a character vector of names of the variables
#' defining clusters for a mixed model with nested or crossed random effects;
#' if missing, cross-validation is performed for individual cases rather than
#' for clusters
#' @param blups a function to be used to compute BLUPs for
#' case-based CV when \code{details = TRUE}.
#' @param fixed.effects a function to be used to compute fixed-effect
#' coefficients for cluster-based CV when \code{details = TRUE}.
#' @param ... for \code{cv()} methods, to match generic,
#' and for \code{cvMixed()}, arguments to be passed to \code{update()}.
#'
#' @details
#' For mixed-effects models, cross-validation can be done by "clusters" or by
#' individual observations. If the former, predictions are based only on fixed
#' effects; if the latter, predictions include the random effects (i.e., are the
#' best linear unbiased predictors or "BLUPS").
#'
#' @seealso \code{\link{cv}}, \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
#' \code{\link[nlme]{lme}}, \code{\link[glmmTMB]{glmmTMB}}
#' @examples
#' library("lme4")
#' # from ?lmer:
#' (fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
#' cv(fm1, clusterVariables="Subject") # LOO CV of clusters
#' cv(fm1, seed=447) # 10-fold CV of cases
#' cv(fm1, clusterVariables="Subject", k=5,
#'    seed=834, reps=3) # 5-fold CV of clusters, repeated 3 times
#'
#' library("nlme")
#' # from ?lme
#' (fm2 <- lme(distance ~ age + Sex, data = Orthodont,
#'             random = ~ 1))
#' cv(fm2) # LOO CV of cases
#' cv(fm2, clusterVariables="Subject", k=5, seed=321) # 5-fold CV of clusters
#'
#' library("glmmTMB")
#' # from ?glmmTMB
#' (m1 <- glmmTMB(count ~ mined + (1|site),
#'                zi=~mined,
#'                family=poisson, data=Salamanders))
#' cv(m1, seed=97816, k=5, clusterVariables="site") # 5-fold CV of clusters
#' cv(m1, seed=34506, k=5) # 5-fold CV of cases

#' @returns
#' The methods \code{cv.merMod()}, \code{cv.lme()}, and \code{cv.glmmTMB()},
#' return objects of class \code{"cv"}, or,
#' if \code{reps > 1}, of class \code{"cvList"} (see \code{\link{cv}()}).

#' @describeIn cv.merMod \code{cv()} method for \code{\link[lme4]{lmer}()} and
#' \code{\link[lme4]{glmer}()} models from the \pkg{lme4} package.
#' @export
cv.merMod <-
  function(model,
           data = insight::get_data(model),
           criterion = mse,
           k = NULL,
           reps = 1L,
           seed,
           details = NULL,
           ncores = 1L,
           clusterVariables,
           blups = coef,
           fixed.effects = lme4::fixef,
           ...) {
    cvMixed(
      model,
      package = "lme4",
      data = data,
      criterion = criterion,
      criterion.name = deparse(substitute(criterion)),
      k = k,
      reps = reps,
      seed = seed,
      details = details,
      ncores = ncores,
      clusterVariables = clusterVariables,
      predict.cases.args = list(
        object = model,
        newdata = data,
        type = "response",
        allow.new.levels = TRUE
      ),
      predict.clusters.args = list(
        object = model,
        newdata = data,
        type = "response",
        re.form = NA,
        allow.new.levels = TRUE
      ),
      blups = blups,
      fixed.effects = fixed.effects,
      ...
    )
  }

#' @describeIn cv.merMod \code{cv()} method for \code{\link[nlme]{lme}()}
#' models from the \pkg{nlme} package.
#' @export
cv.lme <-
  function(model,
           data = insight::get_data(model),
           criterion = mse,
           k = NULL,
           reps = 1L,
           seed,
           details = NULL,
           ncores = 1L,
           clusterVariables,
           blups = coef,
           fixed.effects = nlme::fixef,
           ...) {
    cvMixed(
      model,
      package = "nlme",
      data = data,
      criterion = criterion,
      criterion.name = deparse(substitute(criterion)),
      k = k,
      reps = reps,
      seed = seed,
      details = details,
      ncores = ncores,
      clusterVariables = clusterVariables,
      predict.clusters.args = list(
        object = model,
        newdata = data,
        level = 0
      ),
      predict.cases.args = list(
        object = model,
        newdata = data,
        level = 1
      ),
      blups = blups,
      fixed.effects = fixed.effects,
      ...
    )
  }

#' @describeIn cv.merMod \code{cv()} method for \code{\link[glmmTMB]{glmmTMB}()}
#' models from the \pkg{glmmTMB} package.
#' @export
cv.glmmTMB <-
  function(model,
           data = insight::get_data(model),
           criterion = mse,
           k = NULL,
           reps = 1L,
           seed,
           details = NULL,
           ncores = 1L,
           clusterVariables,
           blups = coef,
           fixed.effects = flattenFixefGlmmTMB,
           ...) {
    if(isFALSE(model$call$doFit)) model <- update(model, doFit = TRUE)
    cvMixed(
      model,
      package = "glmmTMB",
      data = data,
      criterion = criterion,
      criterion.name = deparse(substitute(criterion)),
      k = k,
      reps = reps,
      seed = seed,
      details = details,
      ncores = ncores,
      clusterVariables = clusterVariables,
      predict.cases.args = list(
        object = model,
        newdata = data,
        type = "response",
        allow.new.levels = TRUE
      ),
      predict.clusters.args = list(
        object = model,
        newdata = data,
        type = "response",
        re.form = NA,
        allow.new.levels = TRUE
      ),
      blups = blups,
      fixed.effects = fixed.effects,
      ...
    )
  }

# not exported (or registered):

coef.merMod <- function(object, ...) lme4::fixef(object)

coef.lme <- function(object, ...) nlme::fixef(object)

coef.glmmTMB <- function(object, ...) flattenFixefGlmmTMB(object)

flattenFixefGlmmTMB <- function(model, ...){

  coefs <- glmmTMB::fixef(model)

  cond <- if (length(coefs$cond) > 0){
    cond.coefs <- coefs$cond
    coef.names <- names(cond.coefs)
    coef.names[coef.names == "(Intercept)"] <- "Intercept"
    coef.names <- paste0("cond.", coef.names)
    names(cond.coefs) <- coef.names
    cond.coefs
  } else {
    numeric(0)
  }

  zi <- if (length(coefs$zi) > 0){
    zi.coefs <- coefs$zi
    coef.names <- names(zi.coefs)
    coef.names[coef.names == "(Intercept)"] <- "Intercept"
    coef.names <- paste0("zi.", coef.names)
    names(zi.coefs) <- coef.names
    zi.coefs
  } else {
    numeric(0)
  }

  disp <- if (length(coefs$disp) > 0){
    disp.coefs <- coefs$disp
    coef.names <- names(disp.coefs)
    coef.names[coef.names == "(Intercept)"] <- "Intercept"
    coef.names <- paste0("disp.", coef.names)
    names(disp.coefs) <- coef.names
    disp.coefs
  } else {
    numeric(0)
  }

  c(cond, zi, disp)
}


defineClusters <- function(variables, data) {
  all.variables <- names(data)
  if (any(bad <- !variables %in% all.variables)) {
    stop(
      "The following cluster variable",
      if (sum(bad) > 1L)
        "s are"
      else
        " is",
      " not in the data: ",
      paste(variables[bad], collapse = ", ")
    )
  }
  unique(data[, variables, drop = FALSE])
}

selectCluster <- function(cluster, data) {
  result <- apply(data[, names(cluster), drop = FALSE], 1L,
                  function(x)
                    all(x == cluster))
  if (!any(result))
    stop("there is no such cluster: ",
         paste(paste0(names(cluster), "=", cluster), collapse =
                 ", "))
  result
}

selectClusters <- function(clusters, data) {
  result <- apply(clusters, 1L, selectCluster, data = data)
  apply(result, 1L, any)
}
