#' Cross-Validate Several Models Fit to the Same Data
#'
#' A \code{\link{cv}()} method for an object of class  \code{"modlist"},
#' created by the \code{models()} function. This \code{cv()} method simplifies
#' the process of cross-validating several models on the same set of CV folds.
#' \code{models()} performs some
#' "sanity" checks, warning if the models are of different classes, and
#' reporting an error if they are fit to apparently different data sets or
#' different response variables.
#' @param model a list of regression model objects,
#' created by \code{models()}.
#' @param data (required) the data set to which the models were fit.
#' @param criterion the CV criterion ("cost" or lack-of-fit) function, defaults to
#' \code{\link{mse}}.
#' @param k the number of CV folds; may be omitted, in which case the value
#' will depend on the default for the \code{cv()} method invoked for the
#' individual models.
#' @param reps number of replications of CV for each model.
#' @param seed (optional) seed for R's pseudo-random-number generator,
#' to be used to create the same set of CV folds for all of the models;
#' if omitted, a seed will be randomly generated and saved.
#' @param quietly If \code{TRUE} (the default), simple messages (for example about the
#' value to which the random-number generator seed is set), but not warnings or
#' errors, are suppressed.
#' @param ... for \code{models()}, two or more competing models fit to the
#' the same data; the several models may be named.
#' For \code{cv()}, additional arguments to be passed to the \code{cv()} method
#' applied to each model. For the \code{print()}
#' method, arguments to be passed to the \code{print()} method for
#' the individual model cross-validations. For the \code{plot()},
#' method, arguments to be passed to the base \code{\link[base]{plot}()}
#' function.
#' @param x an object of class \code{"cvModList"} to be printed or plotted.
#' @param y the name of the element in each \code{"cv"} object to be
#' plotted; defaults to \code{"adj CV crit"}, if it exists, or to
#' \code{"CV crit"}.
#' @param xlab label for the x-axis (defaults to blank).
#' @param ylab label for the y-axis (if missing, a label is constructed).
#' @param main main title for the graph (if missing, a label is constructed).
#' @param spread if \code{"range"}, the default, show the range of CV criteria
#' for each model along with their average; if \code{"sd"}, show the average
#' plus or minus 1 standard deviation.
#' @param confint if \code{TRUE} (the default) and if confidence intervals are
#' in any of the \code{"cv"} objects, then plot the confidence intervals around the
#' CV criteria.
#' @param axis.args a list of arguments for the \code{\link{axis}()}
#' function, used to draw the horizontal axis. In addition to
#' the axis arguments given explicitly, \code{side=1} (the horizontal
#' axis) and \code{at=seq(along=x)} (i.e., 1 to the number of models)
#' are used and can't be modified.
#' @param col color for the line and points, defaults to the second
#' element of the color palette; see \code{\link{palette}()}.
#' @param lwd line width for the line (defaults to 2).
#' @return \code{models()} returns a \code{"modList"} object, the
#' \code{cv()} method for which returns a \code{"cvModList"} object.
#' @seealso \code{\link{cv}}, \code{\link{cvMixed}}.
#' @examples
#' data("Duncan", package="carData")
#' m1 <- lm(prestige ~ income + education, data=Duncan)
#' m2 <- lm(prestige ~ income + education + type, data=Duncan)
#' m3 <- lm(prestige ~ (income + education)*type, data=Duncan)
#' (cv.models <- cv(models(m1=m1, m2=m2, m3=m3),
#'                  data=Duncan, seed=7949, reps=5))
#' plot(cv.models)
#' (cv.models.ci <- cv(models(m1=m1, m2=m2, m3=m3),
#'                     data=Duncan, seed=5962, confint=TRUE, level=0.50))
#'                  # nb: n too small for accurate CIs
#' plot(cv.models.ci)
#' @describeIn models create a list of models
#' @export
models <- function(...){
  models <- list(...)
  if (length(models) < 2L) stop("fewer than 2 models to be compared")
  classes <- sapply(models, function(m) class(m)[1L])
  n <- sapply(models, function(m) nrow(insight::get_data(m)))
  if (!all(n[1L] == n[-1L])) {
    stop("models are fit to data sets of differing numbers of cases")
  }
  response <- getResponse(models[[1L]])
  for (i in 2L:length(models)){
    if (!isTRUE(all.equal(response, getResponse(models[[i]]),
                          check.attributes=FALSE))){
      stop("models are not all fit to the same response variable")
    }
  }
  if (length(unique(classes)) > 1L)
    warning("models are not all of the same primary class")
  nms <- names(models)
  if (is.null(nms)) {
    names(models) <- paste0("model.", seq_along(models))
  } else {
    unnamed <- which(nms == "")
    names(models)[unnamed] <- paste0("model.", seq_along(unnamed))
  }
  class(models) <- "modList"
  models
}

#' @describeIn models \code{cv()} method for \code{"modList"} objects
#' @exportS3Method
cv.modList <- function(model, data, criterion=mse, k, reps=1, seed, quietly=TRUE, ...){
  n.models <- length(model)
  if (missing(seed)) seed <- sample(1e6, 1L)
  result <- vector(n.models, mode="list")
  names(result) <- names(model)
  class(result) <- "cvModList"
  for (i in 1L:n.models){
    result[[i]] <- if (missing(k)){
      if (quietly){
        suppressMessages(cv(model[[i]], data=data, criterion=criterion,
                            seed=seed, reps=reps, ...))
      } else {
        cv(model[[i]], data=data, criterion=criterion, seed=seed,
           reps=reps, ...)
      }
    } else {
      if (quietly){
        suppressMessages(cv(model[[i]], data=data, criterion=criterion,
                            k=k, seed=seed, reps=reps, ...))
      } else {
        cv(model[[i]], data=data, criterion=criterion, k=k, seed=seed,
           reps=reps, ...)
      }
    }
  }
  result
}

#' @describeIn models \code{print()} method for \code{"cvModList"} objects
#' @exportS3Method
print.cvModList <- function(x, ...){
  nms <- names(x)
  if (inherits(x[[1L]], "cvList")){
    reps <- length(x[[1L]])
    nms <- paste0(nms, " averaged across ", reps, " replications (with SDs)")
  }
  for (i in seq_along(x)){
    cat(paste0("\nModel ", nms[i], ":\n"))
    if (inherits(x[[i]], "cvList")){
      sumry <- summarizeReps(x[[i]])
      xi <- x[[i]][[1L]]
      nms.sumry <- names(sumry)
      xi[nms.sumry] <- sumry
      print(xi)
    } else {
      print(x[[i]], ...)
    }
  }
  return(invisible(x))
}

#' @describeIn models \code{plot()} method for \code{"cvModList"} objects
#' @importFrom grDevices palette
#' @importFrom graphics abline arrows axis box par points strwidth
#' @importFrom stats na.omit
#' @exportS3Method
plot.cvModList <- function(x, y,
                           spread=c("range", "sd"),
                           confint=TRUE,
                           xlab="",
                           ylab,
                           main,
                           axis.args = list(labels=names(x), las=3L),
                           col=palette()[2L], lwd=2, ...){
  spread <- match.arg(spread)
  if (missing(ylab)){
    ylab <- if (inherits(x[[1L]], "cvList")){
      if (spread == "range"){
        "Cross-Validation Criterion (Average and Range)"
      } else {
        expression("Cross-Validation Criterion (Average"%+-%"SD)")
      }
    } else {
      "Cross-Validation Criterion"
    }
  }
  if (miss.main <- missing(main)){
    main <- "Model Comparison"
    if (inherits(x[[1L]], "cvList")){
      main <- paste(main, "\nAveraged Across",
                    length(x[[1L]]), "Replications")
    }
  }
  if (missing(y)){
    nms <- if (inherits(x[[1L]], "cvList")){
      names(x[[1L]][[1L]])
    } else {
      names(x[[1L]])
    }
    which.crit <- na.omit(match(c("adj CV crit", "CV crit"), nms))
    if (length(which.crit) == 0L){
      stop('can\'t find "adj CV crit" or "CV crit"\n',
           'specify the y argument')
    }
    y <- nms[which.crit[1L]]
  }
  if (isTRUE(axis.args$las == 3L)){
    mai <- par("mai")
    mai[1] <- max(strwidth(axis.args$labels, units="inches")) + 0.5
    save.mai <- par(mai = mai)
    on.exit(par(save.mai))
  }
  if (inherits(x[[1L]], "cvList")){
    ynm <- if (spread == "range") {
      paste(y, "range")
    } else {
      paste("SD", y)
    }
    sumry <- lapply(x, summarizeReps)
    crit <- sapply(sumry, function (x) x[[y]])
    if (spread == "sd"){
      sds <- sapply(sumry, function(x) x[[ynm]])
      min.y <- crit - sds
      max.y <- crit + sds
    } else {
      min.y <- sapply(sumry, function(x) x[[ynm]][1])
      max.y <- sapply(sumry, function(x) x[[ynm]][2])
    }
    plot(c(1L, length(x)), c(min(min.y), max(max.y)),
         xlab=xlab, ylab=ylab,
         main=main, axes=FALSE, type="n")
    xs <- seq(along=x)
    points(xs, crit, type="b", col=col, lwd=lwd)
    arrows(xs, min.y, xs, max.y, length=0.125, angle=90,
           col=col, code=3, lty=1, lwd=1)
  } else {
    crit <- sapply(x, function (x) x[[y]])
    if (y == "adj CV crit"){
      cis <- lapply(x, function(x) x[["confint"]])
      if (confint && !all(sapply(cis, function(ci) is.null(ci)))){
        plot.cis <- TRUE
        if (miss.main){
          main <- "Model Comparison\nwith Confidence Intervals"
        }
        lowers <- sapply(cis, function(ci) {
          low <- ci["lower"]
          if (is.null(low)) NA else low
        })
        uppers <- sapply(cis, function(ci) {
          up <- ci["upper"]
          if (is.null(up)) NA else up
        })
        min.y <- min(lowers, na.rm=TRUE)
        max.y <- max(uppers, na.rm=TRUE)
      }
      else {
        plot.cis <- FALSE
        min.y <- min(crit)
        max.y <- max(crit)
      }
    } else {
      plot.cis <- FALSE
      min.y <- min(crit)
      max.y <- max(crit)
    }
    plot(seq(along=crit), crit, xlab=xlab, ylab=ylab, main=main,
         axes=FALSE, type="b", col=col, lwd=lwd,
         ylim=c(min.y, max.y), ...)
    if (plot.cis){
      xs <- seq_along(cis)
      arrows(xs, lowers, xs, uppers, length=0.125, angle=90,
             col=col, code=3, lty=1, lwd=1)
    }
  }
  abline(h=min(crit), lty=2L, col=col)
  box()
  axis(2)
  axis.args$side <- 1L
  axis.args$at <- seq(along=x)
  do.call(axis, axis.args)
}

#' @export
`[.cvModList` <- function(x, ...){
  result <- NextMethod()
  class(result) <- "cvModList"
  result
}
