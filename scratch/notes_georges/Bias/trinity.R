#
#
# Using numDeriv to obtain numerical Rao(score), Wald and Wilks(LRT) tests
# from the likelihood function for the binomial
# with two parametrizations: probability and log odds
#
#


logLikp <- function(n,Y) {
  #
  # creates the log likelihood for a binomial with n observations and Y successes
  # as a function of probability of a success
  #
  ret <- function(p) {
    Y*log(p) + (n-Y)*log(1-p)
  }
  phat <- Y/n
  attr(ret, 'mle') <- phat
  ret
}
logLiklo <- function(n,Y) {
  #
  # creates the log likelihood for a binomial with n observations and Y successes
  # as a function of the log odds of success
  #
  ret <- function(lo) {
    Y*log(1/(1+exp(-lo))) + (n-Y)*log(1/(1+exp(lo)))
  }
  phat <- Y/n
  attr(ret, 'mle') <- log(phat/(1-phat))
  ret
}

tests <- function(llik, h0) {
  mle <- attr(llik,'mle')
  c(
    wald = -(mle - h0)^2 * hessian(llik, mle),
    wald2 = -(mle - h0)^2 * hessian(llik, h0),
    rao = -(grad(llik, h0)^2) / hessian(llik, h0),
    wilks = 2 * (llik(mle) - llik(h0))
  )
}



llp <- logLikp(20,8)
lllo <- logLiklo(20,8)

tests(llp, .4)                   # should be 0
tests(lllo, log(.4/((1 -.4))))   # should be 0

# Rao and wilks should be the same except for numerical errors in differentiation
tests(llp, .5)                   # should be 0
tests(lllo, log(.5/((1 -.5))))   # should be 0

tests(llp, .3)                   # should be 0
tests(lllo, log(.3/((1 -.3))))   # should be 0

tests(llp, .5)
tests(llp, .3)
tests(lllo, log(.4/(1-.4)))
tests(lllo, log(.5/(1-.5)))
tests(lllo, log(.3/(1-.3)))

library(latticeExtra)
plot(seq(0,1,.01), sapply(seq(0,1,.01), function(x) llp(x)))

h0s <- seq(.001,.999,.001)

ts <- lapply(h0s, function(h0) tests(llp, h0))



help(p=numDeriv)
