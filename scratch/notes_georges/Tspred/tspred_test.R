#'
#'  Testing tspred
#'
#'  Discrepancies with predict.ARIMA seem consistent with
#'  regularization using a Kalman filter
#'
#' We will test regularization by shortening or tapering
#' the pi sequence.  This should work fine with AR series
#' but not with an MA sequence that is close to non-invertible.
#'
#' Variables:
#'
#' - type of series:
#'   - stationary and invertible
#'   - nearly non-stationary 2nd order
#'   - nearly non-invertible 2nd order
#'   - both
#'   X
#'   - integrated or not
#'   X
#'   - seasonal
#' - models
#'   - 202 101 212 111 + same seasonal 12
#' - prediction
#'   - arima
#'   - tspred
#'   - tspred with taper
#'
library(knitr)
opts_chunk$set(error=TRUE,fig.dim=c(8,7))
# hook_plot = knit_hooks$get('plot')
# knit_hooks$set(plot = function(x, options) paste('\n', hook_plot(x, options), sep = ''))
#'
#'
library(cv)

library(latticeExtra)
source('tspred.R')
predplot <- function(N = 9999, ar = numeric(0), ma = numeric(0) ,int = 0, xcoef = 0,
                     half_taper = 50,
                     period = numeric(0), sar = numeric(0), sma = numeric(0),
                     sint = 0, sigma = 1,
                     main = '', seed = NA,
                     ylim = NULL,
                     sub = paste0(
                       "ar=",paste(ar,collapse = ' '),
                       " ma=", paste(ma, collapse = ' '),
                       " int=", int,
                       " period=", sint,
                       "sar=",paste(sar,collapse = ' '),
                       "sma=", paste(sma, collapse = ' '),
                       "sint=", sint))  {
  {
    print(polyroot(c(1,-ar)))
    print(Mod(polyroot(c(1,-ar))))
    print(2*pi/Arg(polyroot(c(1,-ar))))
    if(!is.na(seed)) set.seed(seed)
    dd <- rts(N, ar = ar, ma = ma, int = int,
              period = period, sar = sar, sma = sma, sint = sint, xcoef= xcoef)
    # print(xyplot(y ~ time, dd, type = 'l'))
  }
  models <- within(
    list(),
    {
      `101` <- print(cv::Arima(y ~ x, data = dd, order = c(1,0,1)))
      `202` <- print(cv::Arima(y ~ x, data = dd, order = c(2,0,2)))
      `111` <- print(cv::Arima(y ~ x, data = dd, order = c(1,1,1)))
      `212` <- print(cv::Arima(y ~ x, data = dd, order = c(2,1,2)))
      `101/101(3)` <- print(cv::Arima(y ~ 1 +x, data = dd,
                                order = c(1,0,1),
                                seasonal = list(order = c(1,0,1), period = 3)))
      `202/202(3)` <- print(cv::Arima(y ~ x, data = dd,
                                order = c(2,0,2),
                                seasonal = list(order = c(2,0,2), period = 3)))
      `111/111(3)` <- print(cv::Arima(y ~ 1 +x, data = dd,
                                order = c(1,1,1),
                                seasonal = list(order = c(1,1,1), period = 3)))
      `212/212(3)` <- print(cv::Arima(y ~ x, data = dd,
                                order = c(2,1,2),
                                seasonal = list(order = c(2,1,2), period = 3)))

    }
  ) %>% rev

  Show_pred(models, dd,  last = 30, show = 50,main = main, sub = sub, ylim = ylim,
            half_taper = half_taper)
  print(lapply(models, summary))
}
#'
#'
#'

cols <- c('#88000088','#00880088','#00008888')

trellis.par.set('superpose.symbol',
                list(pch=c(0,2,6,1),
                col = cols ))
trellis.par.set('superpose.line', list(lty=1:4, col = cols ))

#'
#' ## AR .8, MA .8
#'

replicate(3, predplot(ar=.8, ma =.8, main = '101'))

#'
#'
#' ## AR -1.95,-, MA .8
#'

replicate(3,predplot(ar=c(1.95,-.97), ma =.8, main = '201'))


#'
#' ## AR .8, MA .8  PERIOD 3 MAR .8  SMA .8
#'

replicate(3, predplot(
  ar=.8, ma =.8,
  period = 3, sar = .8, sma = .8,
  main = '101/101(3)'))

#'
#' ## AR -1.95,-.97, MA .8, PERIOD 3 SAR -1.95,-.97  SMA .8
#'


replicate(3,
          predplot(
            ar=c(1.95,-.97), ma =.8,
            period = 3, sar=c(1.95,-.97), sma =.8,
            main ='201/201(3)') )


#'
#' ## AR -1.95,-.97, MA .8, INT: 1, PERIOD 3 SAR -1.95,-.97  SMA .8
#'

replicate(3,
          predplot(
            ar=c(1.95,-.97), int = 1, ma =.8,
            period = 3, sar=c(1.95,-.97), sma =.8,
            main='211/201(3)' ))

#'
#' ## AR -1.95,-.97, MA .8, INT: 1, PERIOD 3 SAR -1.95,-.97  SMA .8
#'

replicate(3,predplot(
  ar=c(1.95,-.97),  ma =.8, period = 3,
  sar=c(1.95,-.97), sma =.8, sint = 1,
  main = '201/211(3)'))



#'
knitr::knit_exit()
#'
#' ## AR 1.95, MA .8
#'
#+ error=TRUE

predplot(ar=c(1.95,-.97), ma =.8)
predplot(ar=c(1.95,-.97), ma =.8)
predplot(ar=c(1.95,-.97), ma =.8)




#+ error=TRUE
predplot(ma=c(.9), int = 0)
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
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1, ylim = 2 * c(-2000, 1000))
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1, seed = 431)
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1, seed = 631)
predplot(ar=c(1.95,-.99),ma=c(1.95,-.97), int = 1, seed = 542)

