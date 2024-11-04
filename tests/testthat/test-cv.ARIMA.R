# test that correct models are fit for CV of ARIMA models

Lake <- data.frame(level=LakeHuron, year=time(LakeHuron))
lake.arima <- Arima(level ~ I(year - 1920), data=Lake,
                    order=c(2, 0, 0))

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
