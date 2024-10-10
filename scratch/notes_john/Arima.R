diff.data.frame <- function(x, index, d = 1, ...) {
  y <- x[, index]
  if (!is.numeric(y))
    stop("can't difference a non-numeric variable")
  y <- diff(y, difference = d)
  x <- x[-(1:d), ]
  x[, index] <- y
  x
}

Arima <- function(formula,
                  data,
                  subset = NULL,
                  na.action = na.pass,
                  order = c(1L, 0L, 0L),
                  arima.method = c("arima", "gls"),
                  method = c("REML", "ML"),
                  ...) {

  Arima_arima <- function() {
    which.int <- which("(Intercept)" == colnames(x))
    if (length(which.int > 0))
      x <- x[, -which.int, drop = FALSE]
    result <- list(
      formula = formula,
      data = data,
      subset = subset,
      na.action = na.action,
      order = order,
      call = cl,
      dots = dots
    )
    if (length(formula) == 2) {
      if (!(ncol(x) == 1) && is.numeric(x[, 1]))
        stop("formula must specify a single response")
      result$arima <- stats::arima(x[, 1], order = order, ...)
      result$response <- x[, 1]
    } else {
      y <- model.response(mf, "numeric")
      if (!((is.ts(y) || is.vector(y)) && is.numeric(y)))
        stop("formula must specify a single response")
      result$arima <- stats::arima(y, order = order, xreg = x, ...)
      result$response <- y
      result$model.matrix <- x
    }
    class(result) <- c("ARIMA_arima", "ARIMA")
    result
  }

  Arima_gls <- function() {

    if (!require("nlme", quietly=TRUE))
      stop("nlme package not available")
    y <- model.response(mf, "numeric")
    p <- order[1]
    d <- order[2]
    q <- order[3]
    if (!((is.ts(y) || is.vector(y)) && is.numeric(y)))
      stop("formula must specify a single response")

    if (d > 0) {
      y <- diff(y, difference = d)
      x <- x[-(1:d), , drop=FALSE]
    }
    res <- nlme::gls(
      y ~ x - 1,
      correlation = corARMA(form = ~ 1, p = p, q = q),
      subset = subset,
      na.action = na.action,
      method = method,
      ...
    )
    b <- coefficients(res)
    vcov <- vcov(res)
    rownames(vcov) <- colnames(vcov) <- names(b) <- colnames(x)

    result <- list(
      formula = formula,
      data = data,
      subset = subset,
      na.action = na.action,
      contrasts = res$contrasts,
      order = c(p, d, q),
      call = cl,
      dots = dots,
      response = y,
      model.matrix = x,
      coefficients = b,
      vcov = vcov,
      terms = mt,
      sigma = res$sigma,
      corStruct = res$modelStruct$corStruct,
      fitted = res$fitted,
      residuals = res$residuals,
      logLik = res$logLik,
      method = res$method,
      gls = res
    )
    class(result) <- c("ARIMA_gls", "ARIMA")
    result
  }

  arima.method <- match.arg(arima.method)
  method <- match.arg(method)
  dots <- list(...)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)

  if (arima.method == "arima") {
    Arima_arima()
  } else {
    Arima_gls()
  }

}

print.ARIMA_arima <- function(x, ...){
  xx <- x$arima
  xx$call <- getCall(x)
  print(xx)
  invisible(x)
}

update.ARIMA_arima <- function(object, ...){
  cl0 <- object$call
  cl <- match.call()
  args <- object$dots
  args[c("formula", "data", "subset", "na.action", "order")] <-
    object[c("formula", "data", "subset", "na.action", "order")]
  args$arima.method <- "arima"
  dots <- list(...)
  names <- names(dots)
  # for (name in names){
  #   args[[name]]  <- dots[[name]]
  # }
  args[names] <- dots[names]
  result <- do.call(Arima, args)
  names <- names(cl)
  for (name in names){
    if (name == "" || name == "object") next
    cl0[[name]] <- cl[[name]]
  }
  result$call <- cl0
  result
}

update.ARIMA_gls <- function(object, ...){
  args <- list(...)
  arg.names <- names(args)
  if (!("formula" %in% arg.names)) args$formula <- object$formula
  if (!("data" %in% arg.names)) {
    args$data <- object$data
    data.name <- object$call[["data"]]
    } else {
      data.name <- NULL
    }
  if (!("subset" %in% arg.names)) args$subset <- object$subset
  if (!("na.action" %in% arg.names)) args$na.action <- object$na.action
  if (!("order" %in% arg.names)) args$order <- object$order
  if (!("method" %in% arg.names)) args$method <- object$method
  args$arima.method <- "gls"
  result <- do.call("Arima", args)
  if (!is.null(data.name)){
    result$call[["data"]] <- data.name
  }
    result
}

coef.ARIMA_arima <- function(object, ...) coef(object$arima)
coef.ARIMA_gls <- function(object, ...) {
  b <- object$coefficients
  cor <- capture.output(print(object$corStruct))
  names.b <- names(b)
  names.c <- strsplit(trimws(cor[2]), "[[:space:]]+")[[1]]
  c <- as.numeric(strsplit(trimws(cor[3]), "[[:space:]]+")[[1]])
  result <- c(c, b)
  names(result) <- c(names.c, names.b)
  result
}

model.matrix.ARIMA <- function(object, ...) object$model.matrix

predict.ARIMA_arima <- function(object, n.ahead, newdata = NULL,
                          se.fit = FALSE, ...){
  if (missing(n.ahead) && is.null(newdata)) return(NULL)
  x <- model.matrix(object)
  new.x <- if (!is.null(newdata) && !is.null(model.matrix(object))) {
    n.ahead <- nrow(newdata)
    model.frame(object$formula[-2], data=newdata)
  } else {
    NULL
  }
  predict(object$arima, n.ahead=n.ahead, newxreg=new.x,
          se.fit=se.fit, ...)
}

print.ARIMA_gls <- function(object, ...) {
  result <- object$gls
  result$coefficients <- object$coefficients
  result$varBeta <- object$vcov
  call <- object$call
  names(call)[which(names(call) == "formula")] <- "model"
  result$call <- call
  print(result)
}

summary.ARIMA_gls <- function(object, ...) {
  result <- object$gls
  result$coefficients <- object$coefficients
  result$varBeta <- object$vcov
  call <- object$call
  names(call)[which(names(call) == "formula")] <- "model"
  result$call <- call
  summary(result)
}

if (FALSE) {

  library(cv)

  Lake <- data.frame(level = LakeHuron, year = time(LakeHuron))
  L1 <- diff(Lake, "level", d = 1)
  L2 <- diff(Lake, "level", d = 2)
  head(Lake)
  head(L1)
  head(L2)

  a.a <- Arima(log(level) ~ I(year - 1920),
        order = c(1, 1, 1),
        data = Lake)
  a.g <- Arima(log(level) ~ I(year - 1920) - 1,
        order = c(1, 1, 1),
        data = Lake,
        method = "ML",
        arima.method="gls")
  a.a
  a.g
  coef(a.a)
  coef(a.g)
  predict(a.a, newdata=data.frame(year=c(1973, 1974)))

  a.a <- Arima(level ~ I(year - 1920),
        order = c(2, 0, 0),
        data = Lake)
  a.g <- Arima(level ~ I(year - 1920),
        order = c(2, 0, 0),
        data = Lake,
        method = "ML",
        arima.method = "gls")
  a.a
  a.g
  coef(a.a)
  coef(a.g)
  head(model.matrix(a.a))
  head(model.matrix(a.g))
  update(a.a, order=c(1, 0, 1))
  update(a.g, order=c(1, 0, 1))

  LH <- data.frame(lh = lh)
  a.a <- Arima( ~ lh, data = LH)
  a.g <- Arima(lh ~ 1, data = LH, method = "ML", arima.method="gls")
  a.a
  summary(a.g)
  coef(a.a)
  coef(a.g)
  update(a.a, subset = 1:40)
  update(a.g, subset = 1:40)
  predict(a.a, n.ahead=3)

}
