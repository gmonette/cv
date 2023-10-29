library(cv)
library(MASS)
library(ISLR)

m <- lm(mpg ~ . - name - origin, data=Auto)

stepAIC(m, trace=FALSE)
stepAIC(m, k=log(nrow(Auto)), trace=FALSE) # BIC

stepAIC(lm(mpg ~ . - name - origin, data=Auto[1:50, ]),
           trace=FALSE)
stepAIC(lm(mpg ~ . - name - origin, data=Auto[1:50, ]),
        k=log(nrow(Auto)), trace=FALSE) # BIC

selectStepAIC(Auto, 1:50, model=m)
selectStepAIC(Auto, 51:100, model=m)
selectStepAIC(Auto, model=m)

selectStepAIC(Auto, 1:50, model=m, k.=log(nrow(Auto))) # BIC

cvSelect(selectStepAIC, Auto, k=5, seed=123, model=m)
cvSelect(selectStepAIC, Auto, k=5, seed=321, model=m)
cvSelect(selectStepAIC, Auto, k=5, seed=123, model=m,
         k.=log(nrow(Auto))) # BIC

cvSelect(selectStepAIC, Auto, k=5, seed=321, model=m, reps=5)


system.time(print(cvSelect(selectStepAIC, Auto, seed=123, model=m)))
system.time(print(cvSelect(selectStepAIC, Auto, seed=123, model=m,
                           ncores=5))) # slower
