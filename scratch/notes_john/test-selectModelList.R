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

