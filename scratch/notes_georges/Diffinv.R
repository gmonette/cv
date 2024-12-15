#'
#' Differencing and Antidifferencing time series: vectors and matrices
#'
#' - diffts now works on vectors, ts objects, matrices
#' - the usual properties of 'ts' objects don't work
#'   because these objects can have multiple seasonal
#'   differencing components
#' - Issues: differencing factors
#'

library(magrittr)
library(cv)
library(latticeExtra)

Diff <- function(x, ...) {
  UseMethod('Diff')
}

Diff.default <- function(
    x,
    order = 1 * (seasonal == 0),
    seasonal = 1 * (period != 0),
    period = 0) {
  if(seasonal > 0) {                 # seasonal first
    for(i in 1:seasonal){
      x <- diffts(x, order = 1, period = period)
    }
  }
  if(order > 0) {
    for(i in 1:order){
      x <- diffts(x, order = 1, period = 1)
    }
  }
  x
}


diffts <- function(x, order, period) {
  #
  # works on matrices and vectors
  #
  # - returns a result that is 'order * period' shorter
  #   than input.
  #

  if(order > 1) return( diffts(diffts(x, order = order - 1, period = period), order = 1 ))
  if(order == 0) return(x)

  # order == 1

  ismat <- is.matrix(x)

  # FIXME: do a more informative error message
  stopifnot( (length(x) %% period) == 0, period > 0, length(x) > 0)
  stopifnot(order %% 1 == 0,  order >= 0)

  r <- 0 * cbind(x)[-(1:period),,drop = FALSE]

  atts <- list(
    x = x,
    period = period,
    periods = c(attr(x,'periods'), period),
    offset = if(is.null(attr(x,'offset'))){
      period
    } else {
      period + attr(x,'offset')
    })

  x <- cbind(x)
  nperiods <- nrow(x) / period
  indices_x <- seq(1, nrow(x), period) - 1
  indices_r <- indices_x[-length(indices_x)]

  for(i in 1:period) {
    for(j in 1:ncol(r)){
      r[i+indices_r, j] <- diff(x[i+indices_x, j])
    }
  }
  if(!ismat) r <- r[,]         # return a vector if x was a vector
  for(i in seq_along(atts)) attr(r, names(atts[i])) <- atts[[i]]
  r
}

rediffts <- function(x, like) {
  #
  # apply differencing to x like differencing performed in 'like'
  #
  periods <- attr(like, 'periods')

  ret <- x

  for(p in periods){
    ret <- diffts(ret, order = 1, period = p)
  }
  ret
}

if(FALSE){   # tests

  #
  # if x is a matrix:
  #
  x <- cbind((1:36)^3, 1)
  colnames(x) <- c('one','two')
  rownames(x) <- 1:36
  Diff(x, period = 12)
  Diff(x, period = 12) %>% class
  #
  # if x is a vector
  #
  x <- (1:36)^3
  names(x) <- paste0('x',1:36)
  Diff(x, period = 12)
  Diff(x, period = 12) %>% Diff(period = 3) %>% Diff()
  #
  #
  xx <- Diff(x, order = 2, seasonal = 2, period = 12)
  str(xx)
  Diff(x, seasonal = 3, period = 12) %>% attributes
  # Diff(x, seasonal = 4, period = 12)

  x <- (1:24)^3
  Diff(x) %>% attributes
  diffts(x, order = 1, period = 1) %>% attributes


  Diff(Diff(x))
  Diff(Diff(Diff(x)))
  xx <- Diff(Diff(Diff(x)))


  xx2 <- Diff(x, order = 3)

  all.equal(xx, xx2)

  # Seasonal

  Diff(x, period = 3)

  x1 <- Diff(x, seasonal = 2, period = 3)
  x2 <- Diff(Diff(x, period = 3), period = 3)

  all.equal(x1, x2)

  #
  # Does 'rediffts' do the right thing?
  #
  xm <- cbind(cube =(1:36)^3, square = (1:36)^2)

  xmres <- Diff(xm, order = 2, seasonal = 1, period = 12)

  xmres2 <- rediffts(xm, like = xmres)
  all.equal(xmres, xmres2)


  #
  # Periods must divide into previous periods
  #

  Diff(x, period = 12)
  Diff(x, period = 12) %>% Diff(period = 3) # monthly and quarterly

}
#
# but
#

# Diff(x, period = 12) %>% Diff(period = 5) # doesn't work

# Diff(x, period = 12) %>% Diff(order = 3) # always works if there's anything left

#
# Antiderivative
#

Diffinv <- function(x, ...) {
  UseMethod('Diffinv')
}
Diffinv.default <- function(x, all = TRUE, at = NA, val = NA, ...){
  # idea:
  # Take a differenced object and
  #
  # - just antidifference 'into' the original object by using
  #   the corresponding values derived from the original object
  #   when antidifferencing.  This has the result of returning
  #   the original object.
  # - antidifference on step
  # - antidifference into a different object as one would
  #   when using a model to predict with a different
  #   trigger series
  #
  #
  if(is.na(at)) {
    if(!all) {
      diffinvts(x)  # antidifference one step
    } else {
      depth <- length(attr(x, 'periods'))
      ret <- x
      for(i in 1:depth) {
        ret <- diffinvts(ret)
      }
    }
  }
  ret
}


diffinvts <- function(x, at = 1, value = NA, period = attr(x, 'period')) {
  #
  # period:  only used if antidifferencing a raw vector
  #          i.e. a vector not created by Diff
  #          For vectors created by Diff, the period attribute of x is used

  if(is.null(period)) period <- 1

  if((length(x) %% period) != 0) stop('period must be a divisor of length(x)')

  xp <- attr(x,'x')
  if(is.null(xp)) xp <- rep(0, length(x) + period)

  nperiods <- length(x) / period
  indices_r <- seq(1, length(x) + period, period) - 1
  indices_x <- indices_r[-length(indices_r)]
  r <- numeric(length(x) + period)
  for(i in 1:period) {
    r[i+indices_r] <- cumsum(c(0,x[i+indices_x]))
  }
  indices_xp <- seq(0, period - 1)
  if(!is.na(at)) {
    if(!is.na(value)) {
      r <- r - r[at - attr(r,'offset') + indices_xp] + value
    } else {
      r <- r - r[at + indices_xp] + xp[at + indices_xp]
    }
  }
  attr(r,'x') <- attr(attr(x, 'x'), 'x')
  attr(r,'period') <- attr(attr(x, 'x'), 'period')
  if(is.null(attr(r, 'period'))) attr(r, 'period') <- period
  attr(r,'periods') <- attr(attr(x,'x'), 'periods')
  r
}

if(FALSE){

  xx <- Diff((1:12)^3, period = 3)
  xx2 <- Diff((1:12)^3, order = 2, period = 3)

  xx2 %>% diffinvts
  xx2 %>% diffinvts %>% diffinvts
  xx2 %>% diffinvts %>% diffinvts %>% diffinvts # FIX has period but no periods attr
  #
}
#
# Approach
#
# - Statistical analysis of
#   differenced y ~ differenced xreg
# - produces an ARMA model
# - For prediction using ARMA model
#   given:
#   - a position from which to predict ('at'), and n.ahead (a multiple of largest period)
#   - xregpred: xreg in original x's with n.ahead rows
#   - xreghist and yhist history (a multiple of largest period for past up to at)
# - Algorithm:
#   - difference  cbind(xreghist, xregpred) and yhist like fit
#   - predict using ARMA model
#   - undifference ARMA model prediction
#
#
#

#
# Modelling and predicting
#
# The goal is to be able to apply an arima model
# to predict with new data
#
# To do this we need to:
#
# - difference new Y and newxreg
# - Use ARMA prediction on the result
# - Antidifference to get prediction
# - Note:
#   Only differenced newxreg matters!
#   Antidifferencing uses antidifferenced Ys
#
# Example 1: simple differencing
#

# Need additional tools:
#
# - Use a Ts objects structure to apply to a new object
#   - currently easy to undo but not to redo
#   - design: do we make the object more complicated
#     or just write an extractor function
#
# Let's try this without xreg
#

# Simple predict that can be applied to a new sequence of Y's

arma2psi <- function(ar=0, ma=0, ar.seasonal=0, ma.seasonal=0,
                     period, lag.max=100, trunc.psi=TRUE,
                     trunc.at=0.001, ...){
  #
  # Copied from scratch/notes_john/arma2psi.R
  #
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

arma2pi <- function(ar=0, ma=0, ar.seasonal=0, ma.seasonal=0,
                    period, lag.max=100, trunc.pi=TRUE,
                    trunc.at=0.001, ...){
  #
  # Returns 'pi' weights to predict Y(t+h) recursively from Y(t-k),...,Y(t)
  # where Ys are ARMA(ar,ma)
  #
  # Dual to psi:
  #
  # use -ARMAtoMA(ar = -ma, ma = -ar)  with B&D signs for ma parameters
  #
  #

  lag.max.tot <- if (!(missing(period) || is.na(period))) lag.max*period else lag.max
  Pi <- - ARMAtoMA(ar = -ma, ma = -ar, lag.max=lag.max.tot)
  if (!(missing(period) || is.na(period))) {
    Pi.seasonal <- - ARMAtoMA(ar = - ma.seasonal, ar = - ma.seasonal, lag.max=lag.max)
    Pi <- Pi + as.vector(rbind(matrix(0, period - 1, lag.max),
                               Pi.seasonal))
  }
  if (trunc.pi){
    which.pi <- which(abs(Pi) >= trunc.at)
    if (length(which.pi) == 0) {
      return(numeric(0))
    }
    if (max(which.pi) == lag.max.tot) {
      warning("all ", lag.max.tot, " pi weights retained")
    } else {
      Pi <- Pi[1:max(which.pi)]
    }
  }
  Pi
}

if(FALSE) {  # test inversion
  arma2pi(ma = arma2psi(ar=c(.2,.2)))
  arma2psi(ar = arma2pi(ma=c(-.2,-.2)))

  arma2pi(ma = arma2psi(ar=c(.45,.45)))
  arma2psi(ar = arma2pi(ma=c(-.45,-.45)))

  arma2pi(ma = arma2psi(ar=c(.499,.499)))  # nearly non-stationary
  arma2psi(ar = arma2pi(ma=c(-.499,-.499)))  # nearly non-invertible

}

{

  # sample model data

  set.seed(233123)
  N <- 9999
  eps <- arima.sim(list(ar=.8,ma=.8), n = N)
  x <- rnorm(N)
  dd <- data.frame(y = cumsum(eps + x), x = x)

  # sample prediction data

  Ns <- 42
  n.ahead <- 9
  eps_pred <- c(arima.sim(list(ar=.8,ma=.8), n = Ns), rep(NA, n.ahead))
  x_pred <- rnorm(Ns + n.ahead)
  newdata <- data.frame(y = cumsum(eps_pred + x_pred), x = x_pred)

  modelns <- cv::Arima(y ~ x, data = dd, order = c(1,1,1))
  model <- cv::Arima(y~ 1 +x, data = dd, order = c(1,0,1), seasonal = list(order = c(1,1,1), period = 3))
  model <- cv::Arima(y~ x, data = dd, order = c(1,1,1), seasonal = list(order = c(1,1,1), period = 3))

  # str(model)
  # str(modelns)
}


##
## Rough version of prediction function
##


pred <- function(model, newdata, n.ahead = 1) {
  #
  # Rough version function to see if the idea works
  #
  # Prediction with new predictors, Xs, and a sequence of
  # 'lead-up' response values, Ys, preceding the predicted values
  # requires predictor values for the lead-up sequence as
  # well as for the values to be predicted.
  #
  # 'newdata' is a data.frame similar to the data frame
  # used to fit the 'model'.  It consists of
  # n.leadup row with values for Xs and Y
  # followed by n.ahead rows with values for Xs
  # and NA's for Y.
  #
  # n.leadup must be long enough to allow ordinary
  # and seasonal differencing to the order in 'model'.
  #
  # The function computes predicted values for the Ys
  # that are missing by:
  #
  # 1. Getting the model matrix and the response
  #    vector corresponding to 'newdata' using
  #    the model.
  # 2. Differencing the new model matrix and response
  #    as in 'model'
  # 3. Using the resulting ARMA sequence to predict
  #    n.ahead ARMA Ys using pi weights on the
  #    differenced Ys and Xs.
  # 4. Undifferencing the result, including the
  #    leadup data to get predictions on the original
  #    scale for Ys.
  # The function returns the vector of predicted
  # values and the new data frame with leadup and
  # predicted values.
  #
  # Note that, for cross-validation, this approach can use
  # a model that omits an internal fold but the 'leadup'
  # data precedes the predicted values.  The method
  # does not use distant future data to predict the near future
  # but it can use distant future data to estimate the
  # coefficients used to predict the near future.
  #

  # Get y and xreg matrix from new Ys and xreg
  #

  model_new <- update(model, data = newdata)
  xreg_new <- model.matrix(model_new)
  y_new <- model_new$response

  diff_order <- model$order[2]
  seasonal_diff_order <- model$seasonal$order[2]

  ##FIXME: following should be use truncation of leadup data with a warning
  ##       and extension of prediction xreg followed by truncation

  if(seasonal_diff_order > 0) {
    if((sum(!is.na(y_new)) %% seasonal_diff_order) != 0 ){
      stop('Number of lead-up data rows should be a multiple of seasonal differencing order (',
           seasonal_diff_order,')')
    }
    if((sum(is.na(y_new)) %% seasonal_diff_order) != 0 ) {
      stop('Number of predicted data rows should be a multiple of seasonal differencing order (',
           seasonal_diff_order,')')
    }
  }

  coefs <- coef(model)

  if ("(Intercept)" %in% names(coefs)){   # from cv::Predict.ARIMA
    xreg_new <- if (!is.null(xreg_new)){
      cbind(1, xreg_new)
    } else {
      matrix(1, nrow=length(y_new), ncol=1)
    }
  }
  ar <- coefs[grepl('^ar[0-9]+$',names(coefs))]
  ma <- coefs[grepl('^ma[0-9]+$',names(coefs))]
  sar <- coefs[grepl('^sar[0-9]+$',names(coefs))]
  sma <-  coefs[grepl('^sma[0-9]+$',names(coefs))]
  beta <- coefs[-seq_len(length(ar) + length(ma) + length(sar) + length(sma))]

  # Get differenced Ys and xreg

  y_new_diff <- y_new
  xreg_new_diff <- xreg_new

  if(model$seasonal$order[2] > 0) {
    stopifnot(length(y_new) %% model$seasonal$period == 0)
    stopifnot(sum(is.na(y_new)) %% model$seasonal$period == 0)
    y_new_diff <- Diff(y_new_diff,
                       seasonal = model$seasonal$order[2],
                       period = model$seasonal$period)
    xreg_new_diff <- Diff(xreg_new_diff,
                          seasonal = model$seasonal$order[2],
                          period = model$seasonal$period)
  }
  if(model$order[2] > 0) {
    y_new_diff <- Diff(y_new_diff,
                       order = model$order[2])
    xreg_new_diff <- Diff(xreg_new_diff,
                          order = model$order[2])
  }

  # Demean differenced Y

  mean_y_new_diff <- mean(y_new_diff, na.rm = T)
  y_new_diff_cent <- y_new_diff - mean_y_new_diff

  # Pi weights

  Pi <- arma2pi(ar = ar, ma = ma)

  convolve <- function(y, w) {
    yrev <- rev(y)
    len <- min(length(y),length(w))
    sum(yrev[1:len]*w[1:len])
  }

  topred <- which(is.na(y_new_diff_cent))

  for(i in topred) y_new_diff_cent[i] <- convolve(y_new_diff_cent[1:(i-1)], Pi)

  new_pred <- y_new_diff_cent + mean_y_new_diff + if(!is.null(beta)) xreg_new_diff %*% beta else 0
  y_new_diff[topred] <- new_pred[topred]

  if(diff_order > 0 || seasonal_diff_order > 0){
    y_pred <- Diffinv(y_new_diff)
  } else {
    y_pred <- y_new_diff
  }

  newdata[[as.character(model$call$formula[[2]])]] <- y_pred
  newdata$.predicted <- rep(c(FALSE,TRUE), c(nrow(newdata) - sum(is.na(y_new)), sum(is.na(y_new))))
  attr(newdata,'Pi') <- Pi
  attr(newdata,'Mod') <- Mod(polyroot(c(1, -Pi)))
  newdata
}

if(FALSE){   # test pred

  # sample model data

  N <- 9999
  eps <- arima.sim(list(ar=.9,ma=.8), n = N)
  x <- rnorm(N)
  dd <- data.frame(y = cumsum(eps + x), x = x)

  # sample prediction data

  Ns <- 42
  n.ahead <- 9
  eps_pred <- c(arima.sim(list(ar=.8,ma=.8), n = Ns), rep(NA, n.ahead))
  x_pred <- rnorm(Ns + n.ahead)
  newdata <- data.frame(y = cumsum(eps_pred + x_pred), x = x_pred)

  model <- cv::Arima(y~ 1 +x, data = dd, order = c(1,0,1), seasonal = list(order = c(1,1,1), period = 3))
  model <- cv::Arima(y~ x, data = dd, order = c(1,1,1), seasonal = list(order = c(1,1,1), period = 3))
  model <- cv::Arima(y ~ x, data = dd, order = c(1,1,1))
  model <- cv::Arima(y ~ x, data = dd, order = c(1,0,1))

  # str(model)
  #
  #
  ret <- pred(model, newdata)
  library(latticeExtra)
  head(ret)
  tail(ret)
  ret$time <- 1:nrow(ret)
  xyplot(y ~ time, ret, groups = .predicted) %>% print

  # no differencing

  # model data
  N <- 9999
  eps <- arima.sim(list(ar=.9,ma=.8), n = N)
  x <- rnorm(N)
  dd <- data.frame(y = eps + x, x = x)

  model <- cv::Arima(y ~ x, data = dd, order = c(1,0,1))

  # leadup data and prediction data
  n.ahead <- 4
  eps_pred <- c(arima.sim(list(ar=.8,ma=.8), n = Ns), rep(NA, n.ahead))
  x_pred <- rnorm(Ns + n.ahead)
  newdata <- data.frame(y = cumsum(eps_pred + x_pred), x = x_pred)

  model
  ret <- pred(model, newdata)
  ret$time <- 1:nrow(ret)
  xyplot(y ~ time, ret, groups = .predicted) %>% print

  attr(ret, 'Mod')

  # debug(pred)
  # undebug(pred)
}

