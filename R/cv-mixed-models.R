#' Cross-Validate Mixed-Effects Model
#'
#' A \code{\link{cv}()} method for models of class \code{"merMod"}, fit
#' by the \code{\link[lme4]{lmer}()} and \code{\link[lme4]{glmer}()} functions
#' in the \pkg{lme4} package. The implementation here should be regarded
#' as experimental.
#'
#' @param model a regression model object of class \code{"merMod"}.
#' @param data data frame to which the model was fit (not usually necessary)
#' @param criterion cross-validation criterion function of form \code{f(y, yhat)}
#'        where \code{y} is the observed values of the response and
#'        \code{yhat} the predicted values; the default is \code{\link{mse}}
#'        (the mean-squared error)
#' @param k perform k-fold cross-validation; \code{k}
#' may be a number or \code{"loo"} or \code{"n"} for n-fold (leave-one-out)
#' cross-validation; the default is \code{10} if cross-validating individual
#' cases and \code{"loo"} if cross-validating clusters.
#' @param reps number of times to replicate k-fold CV (default is \code{1})
#' FIXME! not yet implemented
#' @param seed for R's random number generator; optional, if not
#' supplied a random seed will be selected and saved; not needed
#' for n-fold cross-validation
#' @param ncores number of cores to use for parallel computations
#'        (default is \code{1}, i.e., computations aren't done in parallel)
#' @param clusterVariables a character vector of names of the variables
#' defining clusters for a mixed model with nested random effects;
#' if missing, cross-validation is performed for individual cases rather than
#' for clusters
#' @param includeRandom include the random effects in predicting cases
#' in the clusters used to fit the model with each fold deleted (default
#' is \code{TRUE})? Predictions for cases in the omitted clusters for
#' each fold are always based only on the fixed effects.
#' @param ... to match generic
#'
#' @details
#' For mixed-effects models fit by the \code{lmer()} or \code{glmer{}} functions
#' in the **nlme** package, cross-validation can be done by "clusters" or by
#' individual observations. If the former, predictions are based only on fixed
#' effects; if the latter, predictions include the random effects. Only mixed
#' models with fully nested random effects are supported.
#'
#' @seealso \code{\link{cv}}, \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}}
#' @examples
#' library(lme4)
#' # from ?lmer:
#' (fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
#' cv(fm1, clusterVariables="Subject") # LOO CV of clusters
#' cv(fm1, seed=447) # 10-fold CV of cases
#'
#' @export
cv.merMod <- function(model,
                      data=insight::get_data(model),
                      criterion=mse, k, reps=1,
                      seed,
                      ncores=1,
                      clusterVariables,
                      includeRandom=TRUE,
                      ...){

  defineClusters <- function(variables){
    all.variables <- names(data)
    if (any(bad <- !variables %in% all.variables)){
      stop("The following cluster variable", if (sum(bad) > 1) "s are" else " is",
           " not in the data: ", paste(variables[bad], collapse=", "))
    }
    unique(data[, variables, drop=FALSE])
  }

  selectCluster <- function(cluster){
    result <- apply(data[, names(cluster), drop=FALSE], 1,
                    function(x) all(x == cluster))
    if (!any(result)) stop("there is no such cluster: ",
                           paste(paste0(names(cluster), "=", cluster), collapse=", "))
    result
  }

  selectClusters <- function(clusters) {
    result <- apply(clusters, 1, selectCluster)
    apply(result, 1, any)
  }

  f.clusters <- function(i){
    indices.i <- indices[starts[i]:ends[i]]
    index <- selectClusters(clusters[- indices.i, , drop=FALSE])
    model.i <- update(model, data=data[index, ])
    fit.o.i <- predict(model.i, newdata=data, type="response",
                       re.form=if (includeRandom) NULL else NA,
                       allow.new.levels=TRUE)
    fit.i <- fit.o.i[!index]
    c(criterion(y[!index], fit.i), criterion(y, fit.o.i))
  }

  f.cases <- function(i){
    indices.i <- indices[starts[i]:ends[i]]
    model.i <- update(model, data=data[ - indices.i, ])
    fit.o.i <- predict(model.i, newdata=data, type="response",
                       allow.new.levels=TRUE)
    fit.i <- fit.o.i[indices.i]
    c(criterion(y[indices.i], fit.i), criterion(y, fit.o.i))
  }

  y <- getResponse(model)

  if (missing(clusterVariables)){
    n <- nrow(data)
    if (missing(k)) k <- 10
    if (is.character(k)){
      if (k == "n" || k == "loo") {
        k <- n
      }
    }
    f <- f.cases
  } else {
    clusters <- defineClusters(clusterVariables)
    n <- nrow(clusters)
    if (missing(k)) k <- nrow(clusters)
    f <- f.clusters
  }

  if (!is.numeric(k) || length(k) > 1L || k > n || k < 2 || k != round(k)){
    stop("k must be an integer between 2 and number of",
         if (missing(clusterVariables)) "cases" else "clusters")
  }

  if (k != n){
    if (missing(seed)) seed <- sample(1e6, 1L)
    set.seed(seed)
    message("R RNG seed set to ", seed)
  } else {
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
      f(i)
    }
    stopCluster(cl)
  } else {
    result <- matrix(0, k, 2L)
    for (i in 1L:k){
      result[i, ] <- f(i)
    }
  }
  cv <- weighted.mean(result[, 1L], folds)
  cv.full <- criterion(y, predict(model, type="response",
                                  re.form=if (missing(clusterVariables)) NULL else NA))
  adj.cv <- cv + cv.full - weighted.mean(result[, 2L], folds)
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 clusters = if (!missing(clusterVariables)) clusterVariables else NULL,
                 "n clusters" = if (!missing(clusterVariables)) n else NULL
                 )
  class(result) <- "cv"
  result
}
