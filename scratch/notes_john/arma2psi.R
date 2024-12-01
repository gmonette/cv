# Ripley, B. D. (2002). “Time series in R 1.5.0”. R News, 2(2), 2–7.
# https://www.r-project.org/doc/Rnews/Rnews_2002-2.pdf

arma2psi <- function(ar=0, ma=0, ar.seasonal=0, ma.seasonal=0,
                     period, lag.max=100, trunc.psi=TRUE,
                     trunc.at=0.001, ...){
  lag.max.tot <- if (!(missing(period) || is.na(period))) lag.max*period else lag.max
  psi <- ARMAtoMA(ar = ar, ma = ma, lag.max=lag.max.tot)
  if (!(missing(period) || is.na(period))) {
    psi.seasonal <- ARMAtoMA(ar = ar.seasonal, ma = ma.seasonal, lag.max=lag.max)
    psi <- psi + as.vector(rbind(matrix(0, period - 1, lag.max),
                                                     psi.seasonal))
  }
  if (trunc.psi){
    which.psi <- which(abs(psi) >= trunc.at)
    if (length(which.psi) == 0) {
      return(numeric(0))
    }
    if (max(which.psi) == lag.max.tot) {
      warning("all ", lag.max.tot, " psi weights retained")
    } else {
      psi <- psi[1:max(which.psi)]
    }
  }
  psi
}

if (FALSE){
  arma2psi(ar=c(0.6, 0.2))
  ARMAtoMA(ar=c(0.6, 0.2), lag.max=10)
  arma2psi(ma=c(0.6, 0.2))
  arma2psi(ar=c(0.6, 0.2), ma=-c(0.6, 0.2))
  ARMAtoMA(ar=c(0.6, 0.2), ma=-c(0.6, 0.2), lag.max=10)
  arma2psi(ma.seasonal=c(0.6, 0.2), period=4)
  arma2psi(ar.seasonal=c(0.6, 0.2), ma.seasonal=-c(0.6, 0.2), period=4)
  arma2psi(ma=c(0.6, 0.2),
    ar.seasonal=c(0.6, 0.2), ma.seasonal=-c(0.6, 0.2), period=4)
  arma2psi(ar=c(0.6, 0.2),
            ar.seasonal=c(0.6, 0.2), ma.seasonal=-c(0.6, 0.2), period=4)
  arma2psi(ar=c(0.6, 0.2), ma=-c(0.6, 0.2),
            ar.seasonal=c(0.6, 0.2), ma.seasonal=-c(0.6, 0.2), period=4)
  arma2psi(ar=0.999)
  arma2psi(ar.seasonal=0.999, period=4)
}

Predict <- function(model, ...){
  UseMethod("Predict")
}

Predict.ARIMA <- function(model, all.cases=FALSE, n.ahead,
                          newdata, ...){

  predict_i <- function(i){
    start <- max(i - max.lag, 1)
    end <- i - 1
    len <- end - start + 1
    if (len < 1) return(c(res.i=0, yhat.i=NA))
    res.i <- sum(psi[1:len]*residuals[end:start])
    yhat.i <- if (length(b) > 0) res.i + sum(b*X[i, ]) else res.i
    c(res.i=res.i, yhat.i=yhat.i)
  }

  times <- time(model$data)
  ttsp <- tsp(times)
  X <- model.matrix(model)
  residuals <- as.vector(residuals(model))
  response <- as.vector(model$response)
  tsp(residuals) <- ttsp
  class(residuals) <- "ts"
  n <- length(residuals)
  i.na <- which(is.na(residuals))
  residuals[i.na] <- 0
  order <- model$order
  seasonal <- model$seasonal
  if (is.list(seasonal)) {
    period <- seasonal$period
    seasonal <- seasonal$order
  } else {
    period <- NA
  }
  coef <- coef(model)
  coef.nms <- names(coef)
  if ("(Intercept)" %in% coef.nms){
    X <- if (!is.null(X)){
    cbind(1, X)
    } else {
      matrix(1, nrow=n, ncol=1)
    }
  }
  ar <- if (order[1] > 0){
    coef[grepl("^ar[[:digit:]]+$", coef.nms)]
  } else {
    0
  }
  ma <- if (order[3] > 0){
    coef[grepl("^ma[[:digit:]]+$", coef.nms)]
  } else {
    0
  }
  ar.seasonal <- if (seasonal[1] > 0){
    coef[grepl("^sar[[:digit:]]+$", coef.nms)]
  } else {
    0
  }
  ma.seasonal <- if (seasonal[3] > 0){
    coef[grepl("^sma[[:digit:]]+$", coef.nms)]
  } else {
    0
  }

  if (is.na(period) && any(seasonal != 0)){
    period <- ttsp[3]
  }

  diff <- order[2]
  seasonal.diff <- seasonal[2]
  psi <- arma2psi(ar=ar, ma=ma, ar.seasonal=ar.seasonal,
                  ma.seasonal=ma.seasonal, period=period, ...)
  max.lag <- length(psi)
  arma.coefs <- c(ar, ma, ar.seasonal, ma.seasonal)
  arma.coefs <- arma.coefs[arma.coefs != 0]
  b <- coef[-(1:length(arma.coefs))]
  if (missing(n.ahead)){
    yhat <- residuals.psi <- rep(NA, n)
    tsp(yhat) <- tsp(residuals.psi) <- ttsp
    class(yhat) <- class(residuals.psi) <- "ts"
    is <- if (all.cases) 1:n else i.na
    for (i in is) {
      pre.i <- predict_i(i)
      yhat[i] <- pre.i["yhat.i"]
      residuals.psi[i] <- pre.i["res.i"]
      }
  } else {
    new.x <- if (!missing(newdata)){
      if (n.ahead != nrow(newdata))
        stop("'newdata' doesn't have ", n.ahead, " rows")
      model.frame(model$formula[-2L], data=newdata)
    } else if ("(Intercept)" %in% names(coef(model))) {
      matrix(1, nrow=n.ahead, ncol=1)
    } else {
      NULL
    }
    if (!is.null(new.x)){
      X <- rbind(X, new.x)
    }
    residuals.psi <- yhat <- ts(rep(NA, n.ahead),
                                start=ttsp[2L] + deltat(times),
                                frequency=ttsp[3L])
    residuals <- c(residuals, rep(0, n.ahead))
    for (i in 1:n.ahead){
      pre.i <- predict_i(n + i)
      residuals.psi[i] <- pre.i["res.i"]
      yhat[i] <- pre.i["yhat.i"]
    }
  }
  if (diff > 0) {
    start <- if (!missing(n.ahead)){
      length(response) - diff + 1
    } else {
      ## FIXME: should handle each predicted case separately
      1
    }
    yhat[is.na(yhat)] <- 0
#    if (seasonal.diff == 0) {
      y.previous <- response[start:(start + diff - 1)]
      y.previous[is.na(y.previous)] <- mean(response, na.rm=TRUE)
    # } else {
    #   y.previous <- rep(0, diff)
    # }
    yhat <- diffinv(yhat, lag=1, differences=diff, xi=y.previous)[-(1:diff)]
    tsp(yhat) <- tsp(residuals.psi)
    class(yhat) <- c("ts", class(yhat))
  }

  if (seasonal.diff > 0) {
    start <- if (!missing(n.ahead)){
      length(response) - seasonal.diff*period + 1
    } else {
      1
    }
    yhat[is.na(yhat)] <- 0
#    if (diff == 0){
      y.previous <- response[start:(start + seasonal.diff*period - 1)]
      y.previous[is.na(y.previous)] <- mean(response, na.rm=TRUE)
    # } else {
    #   y.previous <- rep(0, seasonal.diff*period)
    # }
    yhat <- diffinv(yhat, lag=period, differences=seasonal.diff,
                    xi=y.previous)[-(1:(seasonal.diff*period))]
    tsp(yhat) <- tsp(residuals.psi)
    class(yhat) <- c("ts", class(yhat))
  }

  list(yhat=yhat, residuals=residuals.psi, psi=psi)

}

if (FALSE){
  Presidents <- ts_data_frame(approval=presidents,
                              time=time(presidents))

  m.pres.ar <- Arima(approval ~ time, data=Presidents,
                  order=c(1, 0, 0), seasonal=c(1, 0, 0))
  Predict(m.pres.ar)

  m.pres.ma <- Arima(approval ~ time, data=Presidents,
                     order=c(0, 0, 1), seasonal=c(0, 0, 1))
  Predict(m.pres.ma)

  m.pres.diff <- Arima(~ approval, data=Presidents,
                     order=c(0, 1, 1))
  summary(m.pres.diff)
  Predict(m.pres.diff)
  Predict(m.pres.diff, n.ahead=5)$yhat
  predict(m.pres.diff$arima, n.ahead=5)$pred

  m.pres.diff.2 <- Arima(~ approval, data=Presidents,
                       order=c(0, 1, 2))
  summary(m.pres.diff.2)
  Predict(m.pres.diff.2, n.ahead=5)$yhat
  predict(m.pres.diff.2$arima, n.ahead=5)$pred

  m.pres.diff.3 <- Arima(~ approval, data=Presidents,
                         order=c(1, 1, 0))
  summary(m.pres.diff.3)
  Predict(m.pres.diff.3, n.ahead=5)$yhat
  predict(m.pres.diff.3$arima, n.ahead=5)$pred

  m.pres.diff.4 <- Arima(~ approval, data=Presidents,
                         order=c(0, 2, 1))
  summary(m.pres.diff.4)
  Predict(m.pres.diff.4, n.ahead=5)$yhat
  predict(m.pres.diff.4$arima, n.ahead=5)$pred

  m.pres.diff.5 <- Arima(~ approval, data=Presidents,
                         order=c(1, 2, 0))
  summary(m.pres.diff.5)
  Predict(m.pres.diff.5, n.ahead=5, trunc.at=0.0001)$yhat
  predict(m.pres.diff.5$arima, n.ahead=5)$pred

  m.pres.diff.6 <- Arima(~ approval, data=Presidents,
                         order=c(0, 0, 0),
                         seasonal=c(0, 1, 1))
  summary(m.pres.diff.6)
  Predict(m.pres.diff.6, n.ahead=8)$yhat
  predict(m.pres.diff.6$arima, n.ahead=8)$pred

  m.pres.diff.7 <- Arima(~ approval, data=Presidents,
                         order=c(0, 0, 0),
                         seasonal=c(1, 1, 0))
  summary(m.pres.diff.7)
  Predict(m.pres.diff.7, n.ahead=8)$yhat
  predict(m.pres.diff.7$arima, n.ahead=8)$pred

  m.pres.diff.8 <- Arima(~ approval, data=Presidents,
                         order=c(0, 1, 1),
                         seasonal=c(0, 1, 1))
  summary(m.pres.diff.8)
  Predict(m.pres.diff.8, n.ahead=8)$yhat
  predict(m.pres.diff.8$arima, n.ahead=8)$pred

  m.pres.ma.1 <- Arima(~ approval, data=Presidents,
                     order=c(0, 0, 1))
  plot(residuals(m.pres.ma.1), Predict(m.pres.ma.1, all=TRUE)$residuals)
  plot(fitted(m.pres.ma.1), Predict(m.pres.ma.1, all=TRUE)$yhat)
  abline(0, 1)
  which(abs(fitted(m.pres.ma.1) -
              Predict(m.pres.ma.1, all=TRUE)$yhat) > 0.001)
  which(is.na(Presidents$approval))

  predict(m.pres.ma.1$arima, n.ahead=5)$pred
  Predict(m.pres.ma.1, n.ahead=5)$yhat

  m.pres.ar.1 <- Arima(~ approval, data=Presidents,
                       order=c(1, 0, 0))
  plot(residuals(m.pres.ar.1), Predict(m.pres.ar.1, all=TRUE)$residuals)
  plot(fitted(m.pres.ar.1), Predict(m.pres.ar.1, all=TRUE)$yhat)
  abline(0, 1)

  predict(m.pres.ar.1$arima, n.ahead=5)$pred
  Predict(m.pres.ar.1, n.ahead=5)$yhat
  predict(m.pres.ar.1$arima, n.ahead=5)$pred -
    Predict(m.pres.ar.1, n.ahead=5)$yhat

  m.pres.arima.2.0.1 <- Arima(~ approval, data=Presidents,
                       order=c(2, 0, 1))
  plot(residuals(m.pres.arima.2.0.1), Predict(m.pres.arima.2.0.1, all=TRUE)$residuals)
  plot(fitted(m.pres.arima.2.0.1), Predict(m.pres.arima.2.0.1, all=TRUE)$yhat)
  abline(0, 1)

  predict(m.pres.arima.2.0.1$arima, n.ahead=5)$pred
  Predict(m.pres.arima.2.0.1, n.ahead=5)$yhat
  predict(m.pres.arima.2.0.1$arima, n.ahead=5)$pred -
    Predict(m.pres.arima.2.0.1, n.ahead=5)$yhat


  # inverse differences commute
  set.seed(123)
  x <- rnorm(10, 10, 1)
  x
  diffinv(diffinv(x, lag=1, differences=2),
          lag=4, differences=1)
  diffinv(diffinv(x, lag=4, differences=1),
          lag=1, differences=2)

  # inverse differences with initial values don't commute
  diffinv(diffinv(x, lag=1, differences=2, xi=1:2),
          lag=4, differences=1, xi=3:6)
  diffinv(diffinv(x, lag=4, differences=1, xi=3:6),
          lag=1, differences=2, xi=1:2)

}
