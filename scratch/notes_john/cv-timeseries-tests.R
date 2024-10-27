library(cv)

D <- data.frame(lh = lh)

res <- Arima(~lh, data=D)
res
summary(res)
predict(res, n.ahead=1)
coef(res)
update(res, data=D[-(1:10), , drop=FALSE])
update(res, order=c(1, 1, 1))
plot(res)

cv.a <- cv(res, fold.type="cumulative")
summary(cv.a)
cv.a$details

summary(cv(res, k=5, fold.type="preceding"))
summary(cv(res, fold.type="window"))

cv.a.p <- cv(res, ncores=2, fold.type="cumulative")
all.equal(cv.a, cv.a.p)

cv.b <- cv(res, lead=5, fold.type="cumulative")
summary(cv.b)
cv.a$details

cv.b.p <- cv(res, lead=5, ncores=2, fold.type="cumulative")
all.equal(cv.b, cv.b.p)

cv.c <- cv(res, lead=1:5, fold.type="cumulative")
summary(cv.c, fold.type="cumulative")
plot(cv.c, fold.type="cumulative")

cv.d <- cv(res, lead=c(2, 4, 5), fold.type="cumulative")
summary(cv.d)
plot(cv.d)
all.equal(cvInfo(cv.c, "CV criterion")[c(2, 4, 5)],
          cvInfo(cv.d, "CV criterion"))

cv.c.p <- cv(res, lead=1:5, ncores=TRUE, fold.type="cumulative")
all.equal(cv.c, cv.c.p, fold.type="cumulative")

cv.a.i <- cv(res, k=5, fold.type="preceding")
summary(cv.a.i)
cv.a.i$details

cv.a.i.p <- cv(res, k=5, fold.type="preceding", ncores=2)
all.equal(cv.a.i, cv.a.i.p)

cv.a.w <- cv(res, fold.type="window")
summary(cv.a.w)
cv.a.w.p <- cv(res, fold.type="window", ncores=2)
all.equal(cv.a.w, cv.a.w.p)

cv.b.w <- cv(res, fold.type="window", lead=5)
summary(cv.b.w)

cv.c.w <- cv(res, fold.type="window", lead=1:5)
summary(cv.c.w)

cv.c.w.p <- cv(res, fold.type="window", lead=1:5, ncores=2)
all.equal(cv.c.w, cv.c.w.p)

DD <- data.frame(level=LakeHuron, year=time(LakeHuron))
res.x <- Arima(level ~ I(year - 1920), data=DD, order=c(2, 0, 0))
res.x
summary(res.x)
summary(Arima(level - mean(level) ~ I(year - 1920) - 1, data=DD, order=c(2, 0, 0)))
predict(res.x, newdata=data.frame(year=c(1973, 1974)))
predict(res.x, n.ahead=2, newdata=data.frame(year=c(1973, 1974)))
predict(res.x)
update(res.x, data=DD[-(1:10), ])
plot(res.x)
plot(Effect("year", res.x, residuals=TRUE))

cv.ax <- cv(res.x, lead=1:5, fold.type="cumulative")
summary(cv.ax)
cv.ax$details
plot(cv.ax)
cv.ax.p <- cv(res.x, lead=1:5, ncores=2, fold.type="cumulative")
all.equal(cv.ax, cv.ax.p)

cv.ax.i <- cv(res.x, k=5, fold.type="preceding")
summary(cv.ax.i)
cv.ax.i$details

cv.ax.i.p <- cv(res.x, k=5, fold.type="preceding", ncores=2)
all.equal(cv.ax.i, cv.ax.i.p)

cv.bx.i <- cv(res.x, fold.type="window", lead=1:5)
cv.bx.i.p <- cv(res.x, fold.type="window", lead=1:5, ncores=2)
all.equal(cv.bx.i, cv.bx.i.p)
summary(cv.bx.i)
plot(cv.bx.i)

res.x2 <- update(res.x, . ~ . + I((year - 1920)^2))
summary(res.x2)
cv(res.x2, lead=1:5)

## ------- folds --------


ff <- folds(100, 10)
ff
fold(ff, 2)


ff <- folds(100, 10, fold.type="cumulative", begin.with=20)
ff
fold(ff, 2)
fold(ff, 2, predicted=TRUE)
fold(ff, 2, predicted=TRUE, lead=5)
fold(ff, 2, predicted=TRUE, lead=1:5)
fold(ff, 9)
fold(ff, 9, predicted=TRUE)
fold(ff, 9, predicted=TRUE, lead=20)

ff <- folds(100, 10, fold.type="preceding")
ff
fold(ff, 2)
fold(ff, 2, predicted=TRUE, lead=1:5)
fold(ff, 9)
fold(ff, 9, predicted=TRUE, lead=1:5)


ff <- folds(100, fold.type="cumulative")
ff
fold(ff, 2)
fold(ff, 2, predicted=TRUE)
fold(ff, 2, predicted=TRUE, lead=5)
fold(ff, 2, predicted=TRUE, lead=1:5)
fold(ff, 74)
fold(ff, 74, predicted=TRUE, lead=1:5)
fold(ff, 75, predicted=TRUE, lead=1:5)
fold(ff, 76, predicted=TRUE, lead=1:5)
fold(ff, 75, predicted=TRUE, lead=2:5)

ff <- folds(100, fold.type="window")
ff
fold(ff, 3)
fold(ff, 3, predicted=TRUE, lead=1:5)
