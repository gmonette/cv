#
#  Testing tspred
#
#  Discrepancies with predict.ARIMA seem consistent with
#  regularization using a Kalman filter
#
library(cv)
library(latticeExtra)
source('tspred.R')
predplot <- function(ar = numeric(0), ma = numeric(0) ,int = 0, xcoef=0, N = 9999,
                     main = '', seed = 321,
                     ylim = NULL,
                     sub = paste("ar=",
                                 paste(ar,collapse = ' '), " ma=",
                                 paste(ma, collapse = ' '), "int=", int, collapse = '')

)  {

  {
    print(polyroot(c(1,-ar)))
    print(Mod(polyroot(c(1,-ar))))
    print(2*pi/Arg(polyroot(c(1,-ar))))
    set.seed(seed)
    dd <- rts(N, ar = ar, ma = ma, int = int, xcoef = xcoef)
  }

  system.time({
    models <- within(
      list(),
      {
        `101` <- cv::Arima(y ~ x, data = dd, order = c(1,0,1))
        `202` <- cv::Arima(y ~ x, data = dd, order = c(2,0,2))
        `111` <- cv::Arima(y ~ x, data = dd, order = c(1,1,1))
        `211` <- cv::Arima(y ~ x, data = dd, order = c(2,1,1))
        `312` <- cv::Arima(y ~ x, data = dd, order = c(3,1,2))
        `101/111(3)` <- cv::Arima(y ~ 1 +x, data = dd, order = c(1,0,1), seasonal = list(order = c(1,1,1), period = 3))
        `111/111(3)` <- cv::Arima(y ~ x, data = dd, order = c(1,1,1), seasonal = list(order = c(1,1,1), period = 3))

      }
    ) %>% rev
  })

  Show_pred(models, dd,  last = 30, show = 50,main = main, sub = sub, ylim = ylim)
  print(lapply(models, summary))
}

  if(require(spida2)) spida2::td(pch=1:10)

predplot(ar=c(1.95,-.99), int = 0)
predplot(ar=c(1.95,-.99), int = 0, seed = 431)
predplot(ar=c(1.95,-.99), int = 0, seed = 531)
predplot(ar=c(1.95,-.99), ma = .95, int = 0, seed = 731)
predplot(ar=c(1.95,-.99), int = 1)
predplot(ar=c(1.95,-.99), int = 1, ylim = 2*c(-2000, 1000))
predplot(ar=c(1.95,-.99), int = 1, seed = 431)
predplot(ar=c(1.95,-.99), int = 1, seed = 631)
predplot(ar=c(1.95,-.99), int = 1, seed = 542)




predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 0)
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1)
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1, ylim = 2*c(-2000, 1000))
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1, seed = 431)
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1, seed = 631)
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1, seed = 542)

