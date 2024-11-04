# test that Arima() and predict.ARIMA() work properly

Lake <- data.frame(level=LakeHuron, year=time(LakeHuron))
lake.arima <- Arima(level ~ I(year - 1920), data=Lake,
                    order=c(2, 0, 0))

  # test model fit

lake.arima.2 <- arima(x = Lake$level, xreg = Lake$year - 1920,
                      order = c(2, 0, 0))
test_that('Arima() fits model correctly', {
  expect_equal(as.vector(coef(lake.arima)),
               as.vector(coef(lake.arima.2)))
})

  # test fitted values

test_that('Arima() fitted values are correct', {
  expect_equal(as.vector(predict(lake.arima)),
               as.vector(Lake$level - lake.arima.2$residuals))
})

  # test predictions

test_that('Arima() predictions are correct', {
  expect_equal(as.vector(predict(lake.arima, n.ahead=3,
                                 newdata=data.frame(year=1973:1975))),
               as.vector(predict(lake.arima.2, n.ahead=3,
                                 newxreg = (1973:1975) - 1920)$pred))
})

# test that correct models are fit for CV of ARIMA models

cv.lake.cum <- cv(lake.arima, fold.type="cumulative")
m.first <- update(lake.arima, data=Lake[1:25, ])
m.last <- update(lake.arima, data=Lake[1:97, ])
test_that('CV of ARIMA fold.type="cumulative" first fold', {
  expect_equal(coef(m.first),
               cv.lake.cum$details$coefficients[[1]])
})
test_that('CV of ARIMA fold.type="cumulative" last fold', {
  expect_equal(coef(m.last),
               cv.lake.cum$details$coefficients[[73]])
})

cv.lake.win <- cv(lake.arima, fold.type="window")
m.first <- update(lake.arima, data=Lake[1:25, ])
m.last <- update(lake.arima, data=Lake[73:97, ])
test_that('CV of ARIMA fold.type="window" first fold', {
  expect_equal(coef(m.first),
               cv.lake.win$details$coefficients[[1]])
})
test_that('CV of ARIMA fold.type="window" last fold', {
  expect_equal(coef(m.last),
               cv.lake.win$details$coefficients[[73]])
})

cv.lake.pre <- cv(lake.arima, fold.type="preceding", k=5)
m.first <- update(lake.arima, data=Lake[1:20, ])
m.last <- update(lake.arima, data=Lake[61:79, ])
test_that('CV of ARIMA fold.type="window" first fold', {
  expect_equal(coef(m.first),
               cv.lake.pre$details$coefficients[[1]])
})
test_that('CV of ARIMA fold.type="window" last fold', {
  expect_equal(coef(m.last),
               cv.lake.pre$details$coefficients[[4]])
})

# test that parallel computations are correct

test_that('CV of ARIMA fold.type="cumulative" parallel', {
  expect_equal(cv.lake.cum,
               cv(lake.arima, fold.type="cumulative", ncores=2))
})

test_that('CV of ARIMA fold.type="window" parallel', {
  expect_equal(cv.lake.win,
               cv(lake.arima, fold.type="window", ncores=2))
})

test_that('CV of ARIMA fold.type="preceding" parallel', {
  expect_equal(cv.lake.pre,
               cv(lake.arima, fold.type="preceding", k=5, ncores=2))
})

# test computation of MSE

  # for "preceding" folds

m.2nd <- update(lake.arima, data=Lake[21:40, ])
m.3rd <- update(lake.arima, data=Lake[41:60, ])
yhat <- matrix(NA, 98, 3)
levels <- matrix(Lake$level, 98, 3)
yhat[matrix(c(21:23, 1:3), 3, 2)] <- predict(m.first, newdata=Lake[21:23, , drop=FALSE])
yhat[matrix(c(41:43, 1:3), 3, 2)] <- predict(m.2nd, newdata=Lake[41:43, , drop=FALSE])
yhat[matrix(c(61:63, 1:3), 3, 2)] <- predict(m.3rd, newdata=Lake[61:63, , drop=FALSE])
yhat[matrix(c(80:82, 1:3), 3, 2)] <- predict(m.last, newdata=Lake[80:82, , drop=FALSE])
cv.lake.pre.3 <- cv(lake.arima, fold.type="preceding", k=5,
                    lead=1:3)
test_that('MSE for CV of ARIMA fold.type="preceding"', {
  expect_equal(colMeans((levels - yhat)^2, na.rm=TRUE),
               as.vector(cv.lake.pre.3$"CV crit"))
})

  # for "cumulative" folds

yhat <- matrix(NA, 98, 3)
for (i in 25:97){
  m <- update(lake.arima, data=Lake[1:i, ])
  last <- if (i < 96) 3 else if (i == 96) 2 else 1
  index <- matrix(c((i + 1):(i + last), 1:last), last, 2)
  yhat[index] <-  predict(m, n.ahead=last,
                          newdata=Lake[(i + 1):(i + last), ])
}
cv.lake.cum.3 <- cv(lake.arima, fold.type="cumulative",
                    lead=1:3)
test_that('MSE for CV of ARIMA fold.type="cumulative"', {
  expect_equal(colMeans((yhat - levels)^2, na.rm=TRUE),
               as.vector(cv.lake.cum.3$"CV crit"))
})

  # for "window" folds

yhat <- matrix(NA, 98, 3)
for (i in 25:97){
  m <- update(lake.arima, data=Lake[(i - 24):i, ])
  last <- if (i < 96) 3 else if (i == 96) 2 else 1
  index <- matrix(c((i + 1):(i + last), 1:last), last, 2)
  yhat[index] <-  predict(m, n.ahead=last,
                          newdata=Lake[(i + 1):(i + last), ])
}
cv.lake.win.3 <- cv(lake.arima, fold.type="window",
                    lead=1:3)
test_that('MSE for CV of ARIMA fold.type="window"', {
  expect_equal(colMeans((yhat - levels)^2, na.rm=TRUE),
               as.vector(cv.lake.win.3$"CV crit"))
})
