as.data.frame.cv <- function(x, row.names, optional, ...) {
  D <- data.frame(fold = 0,
                  criterion = x$"CV crit")
  if (!is.null(x$"adj CV crit")) {
    D <- cbind(D, adjusted.criterion = x$"adj CV crit")
  }
  if (!is.null(x$"full crit")) {
    D <- cbind(D, full.criterion = x$"full crit")
  }
  if (!is.null(x$confint)) {
    D <-
      cbind(D,
            confint.lower = x$confint[1L],
            confint.upper = x$confint[2L])
  }
  if (!is.null(x$coefficients)) {
    coefs <- x$coefficients
    if (!is.matrix(coefs)) {
      coef.names <- names(coefs)
      coef.names[coef.names == "(Intercept)"] <- "Intercept"
      coef.names <- paste0("coef.", coef.names)
      coef.names <- make.names(coef.names, unique = TRUE)
      names(coefs) <- coef.names
      D <- cbind(D, t(coefs))
    }
  }
  if (!is.null(x$details)) {
    D2 <- data.frame(
      fold = seq_along(x$details$criterion),
      criterion = x$details$criterion
    )
    if (!is.null(x$details$coefficients)) {
      D3 <- t(x$details$coefficients[[1L]])
      colnames <- colnames(D3)
      colnames[colnames == "(Intercept)"] <- "Intercept"
      colnames <- paste0("coef.", colnames)
      colnames <- make.names(colnames, unique = TRUE)
      colnames(D3) <- colnames
      rownames(D3) <- 1
      for (i in 2L:length(x$details$coefficients)) {
        D4 <- t(x$details$coefficients[[i]])
        colnames <- colnames(D4)
        colnames[colnames == "(Intercept)"] <- "Intercept"
        colnames <- paste0("coef.", colnames)
        colnames <- make.names(colnames, unique = TRUE)
        colnames(D4) <- colnames
        rownames(D4) <- i
        D3 <- Merge(D3, D4)
      }
      D2 <- cbind(D2, D3)
    }
    D <- Merge(D, D2)
  }
  rownames(D) <- NULL
  criterion <- x$criterion
  if (!is.null(criterion)) {
    colnames(D)[which(colnames(D) == "criterion")] <- criterion
    colnames(D)[which(colnames(D) == "adjusted.criterion")] <-
      paste0("adjusted.", criterion)
    colnames(D)[which(colnames(D) == "full.criterion")] <-
      paste0("full.", criterion)
  }
  class(D) <- c("cvDataFrame", class(D))
  D
}

as.data.frame.cvList <- function(x, row.names, optional, ...) {
  Ds <- lapply(x, as.data.frame)
  D <- cbind(rep = 1, Ds[[1L]])
  for (i in 2L:length((Ds))) {
    D <- Merge(D, cbind(rep = i, Ds[[i]]))
  }
  class(D) <- c("cvListDataFrame", "cvDataFrame", class(D))
  D
}

as.data.frame.cvModList <- function(x, row.names, optional, ...) {
  Ds <- lapply(x, as.data.frame)
  model.names <- names(x)
  D <- cbind(model = model.names[1L], Ds[[1L]])
  for (i in 2L:length((Ds))) {
    D <- Merge(D, cbind(model = model.names[i], Ds[[i]]))
  }
  class(D) <-
    c("cvModListDataFrame",
      "cvListDataFrame",
      "cvDataFrame",
      class(D))
  D
}

Merge <- function(...) {
  Ds <- lapply(list(...), as.data.frame)
  names <- unique(unlist(sapply(Ds, colnames)))
  D <- Ds[[1L]]
  missing <- names[which(!(names %in% colnames(D)))]
  D[, missing] <- NA
  for (i in 2L:length(Ds)) {
    DD <- Ds[[i]]
    missing <- names[which(!(names %in% colnames(DD)))]
    DD[, missing] <- NA
    D <- rbind(D, DD)
  }
  D
}

## Georges: you might start by modifying this:
as.data.frame.coef.mer <- function(x, row.names, optional, ...) {
  components <- names(x)
  D <- x[[1L]]
  colnames <- colnames(D)
  colnames(D) <- make.names(paste0(components[1L], ".", colnames),
                            unique = TRUE)
  if (length(components) == 1L)
    return(as.data.frame(D))
  for (i in 1L:length(components)) {
    D1 <- x[[i]]
    colnames <- colnames(D1)
    colnames(D1) <- make.names(paste0(components[i], ".", colnames),
                               unique = TRUE)
    D <- Merge(D, D1)
  }
  D
}

print.cvDataFrame <- function(x,
                              digits = getOption("digits") - 2L,
                              ...) {
  NextMethod(digits = digits)
}

summary.cvDataFrame <- function(object,
                                formula,
                                fun = mean,
                                include=c("cv", "folds", "all"),
                                ...) {
  # object: inheriting from "cvDataFrame"
  # formula: of the form some.criterion ~ classifying.variable(s)
  # fun: summary funcction to apply
  # include: what to include
  #          "cv" = overall CV results
  #          "folds" = results for individual folds
  #          "all" = everything (usually not sensible)
  include <- match.arg(include)
  if (include == "cv") {
    object <- object[object$fold == 0, ]
  } else if (include == "folds") {
    object <- object[object$fold != 0, ]
  }
  helper <- function(formula, data) {
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mf
  }
  data <- helper(formula, data=object)
  car::Tapply(formula, fun=fun, data = data)
}
