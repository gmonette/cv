library(cv)

debug(selectModelList)
undebug(selectModelList)

debug(cvSelect)
undebug(cvSelect)

data("Auto", package="ISLR2")
for (p in 1:10){
  command <- paste0("m.", p, "<- lm(mpg ~ poly(horsepower, ", p,
                    "), data=Auto)")
  eval(parse(text=command))
}


cv.auto.10 <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 data=Auto, seed=2120)
cv.auto.10
cv.auto.10[which.min(sapply(cv.auto.10, function(x) x$"CV crit"))]

res <- selectModelList(Auto, 1:40,
                       models(m.1, m.2, m.3, m.4, m.5,
                              m.6, m.7, m.8, m.9, m.10))
res

selectModelList(Auto,
                model=models(m.1, m.2, m.3, m.4, m.5,
                             m.6, m.7, m.8, m.9, m.10),
                seed=2120)

res <- cv(selectModelList, Auto,
          working.model=models(m.1, m.2, m.3, m.4, m.5,
                               m.6, m.7, m.8, m.9, m.10),
          save.model=TRUE,
          seed=2120)
res
unclass(res)
cv(m.7, seed=2120)

res.alt <- cv(models(m.1, m.2, m.3, m.4, m.5,
          m.6, m.7, m.8, m.9, m.10),
   data=Auto,
   seed=2120,
   recursive=TRUE,
   save.model=TRUE)
res.alt$selected.model
res.alt
all.equal(res, res.alt)

res.par <- cv(selectModelList, Auto,
              working.model=models(m.1, m.2, m.3, m.4, m.5,
                                   m.6, m.7, m.8, m.9, m.10),
              save.model=TRUE,
              seed=2120,
              ncores=2)

all.equal(res, res.par)
res
res.par

res.repl <- cv(selectModelList, Auto,
              working.model=models(m.1, m.2, m.3, m.4, m.5,
                                   m.6, m.7, m.8, m.9, m.10),
              save.model=TRUE,
              seed=2120,
              reps=5)
res.repl

res <- cv(selectModelList, Auto,
          working.model=models(m.1, m.2, m.3, m.4, m.5,
                               m.6, m.7, m.8, m.9, m.10),
          save.model=TRUE,
          seed=2120)
res

res10 <- cv(selectModelList, Auto, k=10,
          working.model=models(m.1, m.2, m.3, m.4, m.5,
                               m.6, m.7, m.8, m.9, m.10),
          save.model=TRUE,
          seed=2120)
res10
all.equal(res, res10)

res10.10 <- cv(selectModelList, Auto, k=10, k.recurse=10,
            working.model=models(m.1, m.2, m.3, m.4, m.5,
                                 m.6, m.7, m.8, m.9, m.10),
            save.model=TRUE,
            seed=2120)
res10.10
all.equal(res, res10.10)

res.cv.10 <- cv(models(m.1, m.2, m.3, m.4, m.5,
                       m.6, m.7, m.8, m.9, m.10),
                Auto, k=10,
                recursive=TRUE,
                save.model=TRUE,
                seed=2120)
all.equal(res, res.cv.10)

res.loo <- cv(selectModelList, Auto, k="loo",
               working.model=models(m.1, m.2, m.3, m.4, m.5,
                                    m.6, m.7, m.8, m.9, m.10),
               save.model=TRUE)
res.loo

res.loo.loo <- cv(selectModelList, Auto, k="loo", k.recurse="loo",
              working.model=models(m.1, m.2, m.3, m.4, m.5,
                                   m.6, m.7, m.8, m.9, m.10),
              save.model=TRUE)
all.equal(res.loo, res.loo.loo)

res.cv.loo <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 recursive=TRUE,
                 k="loo",
                 save.model=TRUE)
all.equal(res.loo, res.cv.loo)


res.10.loo <- cv(selectModelList, Auto, k=10, k.recurse="loo",
              working.model=models(m.1, m.2, m.3, m.4, m.5,
                                   m.6, m.7, m.8, m.9, m.10),
              save.model=TRUE,
              seed=123)
res.10.loo

res.cv.10.loo <- cv(models(m.1, m.2, m.3, m.4, m.5,
                        m.6, m.7, m.8, m.9, m.10),
                 recursive=TRUE,
                 k=10, k.recurse="loo",
                 save.model=TRUE,
                 seed=123)
all.equal(res.10.loo, res.cv.10.loo)


system.time(res <- cv(selectModelList, Auto,
          working.model=models(m.1, m.2, m.3, m.4, m.5,
                               m.6, m.7, m.8, m.9, m.10),
          save.model=TRUE,
          seed=2120))
res

system.time(res.par <- cv(selectModelList, Auto,
              working.model=models(m.1, m.2, m.3, m.4, m.5,
                                   m.6, m.7, m.8, m.9, m.10),
              save.model=TRUE,
              seed=2120,
              ncores=2))
all.equal(res, res.par)


system.time(res.loo <- cv(selectModelList, Auto, k="loo",
          working.model=models(m.1, m.2, m.3, m.4, m.5,
                               m.6, m.7, m.8, m.9, m.10),
          save.model=TRUE))
res.loo

system.time(res.par.loo <- cv(selectModelList, Auto, k="loo",
              working.model=models(m.1, m.2, m.3, m.4, m.5,
                                   m.6, m.7, m.8, m.9, m.10),
              save.model=TRUE,
              ncores=2))
all.equal(res.loo, res.par.loo)
