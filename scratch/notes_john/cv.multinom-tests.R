library(nnet)
data("BEPS", package="carData")
m.beps <- multinom(vote ~ age + gender + economic.cond.national +
                     economic.cond.household + Blair + Hague + Kennedy +
                     Europe*political.knowledge, data=BEPS)
BayesRuleMulti <- function(y, yhat){
  result <- mean(y != yhat)
  attr(result, "casewise loss") <- "y != yhat"
  result
}

GetResponse.multinom <- function(model, ...) {
  insight::get_response(model)
}

m.beps <- update(m.beps, trace=FALSE)

(cv1 <- cv(m.beps, seed=3465))
(cv2 <- cv(m.beps, seed=3465, ncores=2))
all.equal(cv1, cv2)

system.time(print(cv1 <- cv(m.beps, k="loo")))
system.time(print(cv2 <- cv(m.beps, k="loo", ncores=2)))
all.equal(cv1, cv2)
