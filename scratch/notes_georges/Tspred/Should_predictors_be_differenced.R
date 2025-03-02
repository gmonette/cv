#'
#' Should the predictors be differenced?
#'
#' Answer: No. The Arima model fits a model in which the
#' the residuals of the regression have an ARIMA
#' structure.
#'
library(cv)
dd <- within(
  data.frame(time = 1:10000),
  {
    x <- rnorm(10000)
    ygen <- arima.sim(model = list(order = c(1,0,1), ar = .9, ma = .7, sd =.1), n = 10000)
    y1 <- cumsum(ygen + x)

    y2 <- cumsum(ygen) + x

  }
)


head(dd)

fit1 <- Arima(y1 ~ x, dd, order = c(1,1,1))
fit2 <- Arima(y2 ~ x, dd, order = c(1,1,1))
str(fit1)
fit1$arima
summary(fit1)
summary(fit2)
dd$r1 <- residuals(fit1)
dd$r2 <- residuals(fit2)

fitr1 <- Arima(~ r1, order = c(1,1,1), dd)
fitr2 <- Arima(~ r2, order = c(1,1,1), dd)
summary(fitr1)
summary(fitr2)

