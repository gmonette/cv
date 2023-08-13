#' @describeIn cv merMod method
#' @param clusterVariables a character vector of names of the variables
#' defining clusters for a mixed model with nested random effects
#' @export
cv.merMod <- function(model,
                      data=insight::get_data(model),
                      criterion=mse, k=nrow(clusters),
                      seed,
                      ncores=1,
                      clusterVariables,
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

  f <- function(i){
    indices.i <- indices[starts[i]:ends[i]]
    index <- selectClusters(clusters[- indices.i, , drop=FALSE])
    model.i <- update(model, data=data[index, ])
    fit.o.i <- predict(model.i, newdata=data, re.form=NA) # allow.new.levels=TRUE)
    fit.i <- fit.o.i[!index]
    c(criterion(y[!index], fit.i), criterion(y, fit.o.i))
  }

  clusters <- defineClusters(clusterVariables)

  y <- insight::get_response(model)
  n <- nrow(clusters)
  if (!is.numeric(k) || length(k) > 1L || k > n || k < 2 || k != round(k)){
    stop("k must be an integer between 2 and number of clussters")
  }
  if (k != n){
    if (missing(seed)) seed <- sample(1e6, 1L)
    set.seed(seed)
    message("R RNG seed set to ", seed)
  } else {
    seed <- NULL
  }
  nk <-  n %/% k # number of clusters in each fold
  rem <- n %% k  # remainder
  folds <- rep(nk, k) + c(rep(1, rem), rep(0, k - rem)) # allocate remainder
  ends <- cumsum(folds) # end of each fold
  starts <- c(1, ends + 1)[-(k + 1)] # start of each fold
  indices <- if (n > k) sample(n, n)  else 1:n # permute clusters

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
  cv.full <- criterion(y, predict(model, re.form=NA))
  adj.cv <- cv + cv.full - weighted.mean(result[, 2L], folds)
  result <- list("CV crit" = cv, "adj CV crit" = adj.cv, "full crit" = cv.full,
                 "k" = if (k == n) "n" else k, "seed" = seed,
                 clusters = clusterVariables, "n clusters" = n)
  class(result) <- "cv"
  result
}
