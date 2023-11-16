#' Cross-Validate Mixed-Effects Model
#'
#' \code{\link{cv}()} methods for models of class \code{"merMod"}, fit
#' by the \code{\link[lme4]{lmer}()} and \code{\link[lme4]{glmer}()} functions
#' in the \pkg{lme4} package; for models of class \code{"lme"}
#' fit by the \code{\link[nlme]{lme}()} function in the \pkg{nlme}
#' package; and for models of class \code{"glmmTMB"} fit by the
#' \code{\link[glmmTMB]{glmmTMB}()} function in the \pkg{glmmTMB} package.
#' The implementations here should be regarded
#' as experimental. The \code{cvMixed()} function is meant to be called by
#' \code{cv()} methods for mixed-effect models and not directly by the user.
#' It can be used to extend \code{cv()} to other classes of mixed-effects models.
#'
#' @param model a mixed-effects model object for which a \code{cv()} method is available.
#' @param fun the mixed-modeling function employed.
#' @param data data frame to which the model was fit (not usually necessary)
#' @param criterion cross-validation ("cost" or lack-of-fit) criterion function of form \code{f(y, yhat)}
#'        where \code{y} is the observed values of the response and
#'        \code{yhat} the predicted values; the default is \code{\link{mse}}
#'        (the mean-squared error)
#' @param k perform k-fold cross-validation; \code{k}
#' may be a number or \code{"loo"} or \code{"n"} for n-fold (leave-one-out)
#' cross-validation; the default is \code{10} if cross-validating individual
#' cases and \code{"loo"} if cross-validating clusters.
#' @param reps number of times to replicate k-fold CV (default is \code{1})
#' @param seed for R's random number generator; optional, if not
#' supplied a random seed will be selected and saved; not needed
#' for n-fold cross-validation
#' @param ncores number of cores to use for parallel computations
#'        (default is \code{1}, i.e., computations aren't done in parallel)
#' @param clusterVariables a character vector of names of the variables
#' defining clusters for a mixed model with nested random effects;
#' if missing, cross-validation is performed for individual cases rather than
#' for clusters
#' @param predict.clusters.args a list of arguments to be used to predict
#' the whole data set from a mixed model when performing CV on clusters;
#' the first two elements should be
#' \code{model} and \code{newdata}; see the "Extending the cv package" vignette
#' (\code{vignette("cv-extend", package="cv")}).
#' @param predict.cases.args a list of arguments to be used to predict
#' the whole data set from a mixed model when performing CV on cases;
#' the first two elements should be
#' \code{model} and \code{newdata}; see the "Extending the cv package" vignette
#' (\code{vignette("cv-extend", package="cv")}).
#'
#' @param ... for \code{cv()} methods, to match generic,
#' and for \code{cvMixed()}, arguments to be passed to \code{update()}.
#'
#' @details
#' For mixed-effects models, cross-validation can be done by "clusters" or by
#' individual observations. If the former, predictions are based only on fixed
#' effects; if the latter, predictions include the random effects (i.e., are the
#' best linear unbiased predictors or "BLUPS").
#' Only mixed models with fully nested random effects are supported.
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
#' library(nlme)
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
#' cv(m1, seed=97816, k=5) # 5-fold CV of clusters
#' cv(m1, seed=34506, k=5) # 5-fold CV of cases
#' @describeIn cvMixed not to be called directly
#' @export
cvMixed <- function(model,
                    fun,
                    data=insight::get_data(model),
                    criterion=mse,
                    k,
                    reps=1,
                    seed,
                    ncores=1,
                    clusterVariables,
                    predict.clusters.args=list(object=model, newdata=data),
                    predict.cases.args=list(object=model, newdata=data),
                    ...){

  f.clusters <- function(i, predict.clusters.args, predict.cases.args, ...){
    indices.i <- indices[starts[i]:ends[i]]
    index <- selectClusters(clusters[- indices.i, , drop=FALSE], data=data)
    update.args <- list(...)
    update.args$object <- model
    update.args$data <- data[index, ]
    predict.clusters.args$object <- do.call(update, update.args, envir=environment(fun))
    fit.all.i <- do.call(predict, predict.clusters.args)
    fit.i <- fit.all.i[!index]
    c(criterion(y[!index], fit.i), criterion(y, fit.all.i))
  }

  f.cases <- function(i, predict.clusters.args, predict.cases.args, ...){
    indices.i <- indices[starts[i]:ends[i]]
    update.args <- list(...)
    update.args$object <- model
    update.args$data <- data[ - indices.i, ]
    predict.cases.args$object <- do.call(update, update.args, envir=environment(fun))
    fit.all.i <- do.call(predict, predict.cases.args)
    fit.i <- fit.all.i[indices.i]
    c(criterion(y[indices.i], fit.i), criterion(y, fit.all.i))
  }

  y <- getResponse(model)

  if (missing(clusterVariables)) clusterVariables <- NULL
  if (is.null(clusterVariables)){
    n <- nrow(data)
    if (missing(k)) k <- 10
    if (is.character(k)){
      if (k == "n" || k == "loo") {
        k <- n
      }
    }
    f <- f.cases
  } else {
    clusters <- defineClusters(clusterVariables, data)
    n <- nrow(clusters)
    if (missing(k)) k <- nrow(clusters)
    f <- f.clusters
  }

  if (!is.numeric(k) || length(k) > 1L || k > n || k < 2 || k != round(k)){
    stop("k must be an integer between 2 and number of",
         if (is.null(clusterVariables)) "cases" else "clusters")
  }

  if (k != n){
    if (missing(seed)) seed <- sample(1e6, 1L)
    set.seed(seed)
    message("R RNG seed set to ", seed)
  } else {
    if (reps > 1) stop("reps should not be > 1 for n-fold CV")
    if (!missing(seed) && !is.null(seed)) message("Note: seed ignored for n-fold CV")
    seed <- NULL
  }

  nk <-  n %/% k # number in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  indices <- if (n > k) sample(n, n)  else 1:n # permute clusters/cases

  if (ncores > 1L){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    result <- foreach(i = 1L:k, .combine=rbind) %dopar% {
      f(i, predict.clusters.args, predict.cases.args, ...)
    }
    stopCluster(cl)
  } else {
    result <- matrix(0, k, 2L)
    for (i in 1L:k){
      result[i, ] <- f(i, predict.clusters.args, predict.cases.args, ...)
    }
  }
  cv <- weighted.mean(result[, 1L], folds)
  args <-
    cv.full <- criterion(y,
                         do.call(predict,
                                 if (is.null(clusterVariables)) {
                                   predict.cases.args} else {
                                     predict.clusters.args
                                   }
                         )
    )
  adj.cv <- cv + cv.full - weighted.mean(result[, 2L], folds)
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 clusters = clusterVariables,
                 "n clusters" = if (!is.null(clusterVariables)) n else NULL
  )
  class(result) <- "cv"
  if (reps == 1) {
    return(result)
  } else {
    res <- cv(model=model, data=data, criterion=criterion,
              k=k, ncores=ncores, reps=reps - 1,
              clusterVariables=clusterVariables,
              ...)
    if (reps  > 2){
      res[[length(res) + 1]] <- result
    } else {
      res <- list(res, result)
    }
    for (i in 1:(length(res) - 1)){
      res[[i]]["criterion"] <- res[[length(res)]]["criterion"]
    }
    class(res) <- "cvList"
    return(res)
  }
}

#' @returns \code{cvMixed()}, and functions based on it, such as the methods
#' \code{cv.merMod()}, \code{cv.lme()}, and \code{cv.glmmTMB()}, return objects of class \code{"cv"}, or,
#' if \code{reps > 1}, of class \code{"cvList"} (see \code{\link{cv}()}).

#' @describeIn cvMixed \code{cv()} method
#' @export
cv.merMod <- function(model, data = insight::get_data(model), criterion = mse,
                      k, reps = 1, seed, ncores = 1, clusterVariables, ...){
  cvMixed(
    model,
    fun=lme4::lmer, # works for glmer() too
    data=data,
    criterion=criterion,
    k=k,
    reps=reps,
    seed=seed,
    ncores=ncores,
    clusterVariables=clusterVariables,
    predict.cases.args=list(object=model,
                            newdata=data,
                            type="response",
                            allow.new.levels=TRUE),
    predict.clusters.args=list(object=model,
                               newdata=data,
                               type="response",
                               re.form=NA,
                               allow.new.levels=TRUE),
    ...)
}

#' @describeIn cvMixed \code{cv()} method
#' @export
cv.lme <- function(model, data = insight::get_data(model), criterion = mse,
                   k, reps = 1, seed, ncores = 1, clusterVariables, ...){
  cvMixed(
    model,
    fun=nlme::lme,
    data=data,
    criterion=criterion,
    k=k,
    reps=reps,
    seed=seed,
    ncores=ncores,
    clusterVariables=clusterVariables,
    predict.clusters.args=list(object=model,
                               newdata=data,
                               level=0),
    predict.cases.args=list(object=model,
                            newdata=data,
                            level=1),
    ...)
}

#' @describeIn cvMixed \code{cv()} method
#' @export
cv.glmmTMB <- function(model, data = insight::get_data(model), criterion = mse,
                       k, reps = 1, seed, ncores = 1, clusterVariables, ...){
  cvMixed(
    model,
    fun=glmmTMB::glmmTMB,
    data=data,
    criterion=criterion,
    k=k,
    reps=reps,
    seed=seed,
    ncores=ncores,
    clusterVariables=clusterVariables,
    predict.cases.args=list(object=model,
                            newdata=data,
                            type="response",
                            allow.new.levels=TRUE),
    predict.clusters.args=list(object=model,
                               newdata=data,
                               type="response",
                               re.form=NA,
                               allow.new.levels=TRUE),
    ...)
}


defineClusters <- function(variables, data){
  all.variables <- names(data)
  if (any(bad <- !variables %in% all.variables)){
    stop("The following cluster variable", if (sum(bad) > 1) "s are" else " is",
         " not in the data: ", paste(variables[bad], collapse=", "))
  }
  unique(data[, variables, drop=FALSE])
}

selectCluster <- function(cluster, data){
  result <- apply(data[, names(cluster), drop=FALSE], 1,
                  function(x) all(x == cluster))
  if (!any(result)) stop("there is no such cluster: ",
                         paste(paste0(names(cluster), "=", cluster), collapse=", "))
  result
}

selectClusters <- function(clusters, data) {
  result <- apply(clusters, 1, selectCluster, data=data)
  apply(result, 1, any)
}

