# Ripley, B. D. (2002). “Time series in R 1.5.0”. R News, 2(2), 2–7.
# https://www.r-project.org/doc/Rnews/Rnews_2002-2.pdf


arma2psi <- function(ar=0, ma=0, ar.seasonal=0, ma.seasonal=0,
                     period, lag.max=100, trunc.psi=TRUE,
                     trunc.at=0.001){
  lag.max.tot <- if (!missing(period)) lag.max*period else lag.max
  psi <- ARMAtoMA(ar = ar, ma = ma, lag.max=lag.max.tot)
  if (!missing(period)) {
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
