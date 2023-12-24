#'
#' ## General issues
#'
#'
#' - [Information criteria vs cross validation - data analysis / model validation - Datamethods Discussion Forum](https://discourse.datamethods.org/t/information-criteria-vs-cross-validation/2375)
#'   - Interesting comments by Frank Harrell
#' - Blog post by Rob Hyndman 2010: [Rob J Hyndman - Why every statistician should know about cross-validation](https://robjhyndman.com/hyndsight/crossvalidation/)
#' - 2010 survey: [A survey of cross-validation procedures for model selection](https://projecteuclid.org/journals/statistics-surveys/volume-4/issue-none/A-survey-of-cross-validation-procedures-for-model-selection/10.1214/09-SS054.full)
#'


#
# Davison Hinkley
#
# foldWise vs caseWise computation
#
# Comparing DH equation 6.47
# and weighted average of 6.43
#
#

Dfun <- function(Fa, Ft, fit, prederr = function(y,ypred) {(y-ypred)^2}) {    #  (6.39)
  #
  # Aggregate prediction error: Empirical version of 6.39 of DH
  #
  ypred <- predict(update(fit, data = Ft), newdata = Fa)
  yname <- formula(fit)[[2]]
  y <- Fa[[yname]]
  prederrs <- mapply(prederr, y, ypred)
  mean(prederrs)
}

# test

head(mtcars)

fit <- lm(mpg ~ disp + wt, mtcars)
summary(fit)
mean(resid(fit)^2)

Dfun(mtcars, mtcars, fit)   # produces 'apparent' error                           (6.41)

#
# expected excess error                                                           (6.42)
#
# e(F) = E{ D(F,Ffit) - D(Ffit, Ffit) }
#      = Delta(F) - E{ D(Ffit, Ffit) }
#
#
# Cross-validation
#
# Using training set, Ft and 'assessment set' Fa
#
# Dfun(Fa, Ft, fit)                                                               (6.43)
#
# Note that Dfun is an average of case-wise values of 'prederr'
#
# K-fold CV
#
cvk <- function(fit,
                prederr = function(y, ypred) (y - ypred)^2,
                K = min(10, round(sqrt(N))))
  {

  # Caution: Folds are not randomized

  # generate list of indices for K-fold exclusion
  {
    data <- insight::get_data(fit)
    N <- nrow(data)
    if(K == 'loo') K <-N
    n <- floor(N/K)
    remain <- N - n*K
    sizes <- rep(n,K)
    sizes[seq_len(remain)] <- n+1
    stopifnot(sum(sizes) == N)
    ends <- cumsum(sizes)
    starts <- c(1,ends+1)[-(length(ends) + 1)]
    folds <- lapply(1:length(ends), function(i) starts[i]:ends[i])
  }

  #
  # Apply 6.47 of DH  caseWise computation
  #

  ypred <- rep(NA, N)

  for(fold in folds) {
    ypred[fold] <- predict(update(fit, data = data[-fold, ]), newdata = data[fold, ])
  }

  yname <- formula(fit)[[2]]
  y <- data[[yname]]
  prederrs <- mapply(prederr, y, ypred)

  ret <- c('Delhat_CV,K (6.47)' = mean(prederrs))                             # (6.47)

  #
  # weighted foldWise computation
  #

  foldWiseErrs <-sapply(folds, function(fold) Dfun(data[fold, ], data[-fold, ], fit, prederr))
  foldWiseNs <- sapply(folds, length)
  foldWiseErr <- weighted.mean(foldWiseErrs, foldWiseNs)

  ret <- c(ret, 'foldWiseWeightedErr' = foldWiseErr)

  #
  # Apparent error
  #

  ret <- c(ret, 'Apparent Error (6.41)' = Dfun(data, data, fit, prederr))

  #
  # Weighted cross validation predictions on all responses
  #

  foldWiseOnAllErrs <- sapply(folds, function(fold) Dfun(data, data[-fold, ], fit, prederr))

  ret <- c(ret, 'foldWiseOnAllErrs' = weighted.mean(foldWiseOnAllErrs, foldWiseNs))


#  ret <- c(ret, 'Adjustment' = ret[['Apparent Error (6.41)']] - ret[['foldWiseOnAllErrs']])

  ret <- c(ret, 'Delhat_ACV,K (6.48)' = ret[['Delhat_CV,K (6.47)']] + ret[['Apparent Error (6.41)']] - ret[['foldWiseOnAllErrs']])

  class(ret) <- 'cvk'
  ret

}

print.cvk <- function(x) {
  invisible(print(cbind(Estimate = x)))
}

cvallk <- function(fit,
                   prederr = function(y, ypred) (y - ypred)^2,
                   K = 2:N) {
  N <- nrow(insight::get_data(fit))
  ret <- t(sapply(K, function(k) cvk(fit, prederr = prederr, K =k )))
  class(ret) <- 'cvallk'
  ret
}


plot.cvallk <- function(x, lty = 1:5, col = 1:6, lwd = 2, perm = c(1,2,5,4,3) ) {
  matplot(x[,perm], type = 'l', lwd = lwd, lty = lty, col = col,
          ylab = 'Estimate')
  legend('topright', legend = colnames(ret)[perm], lty = lty, col = col, lwd = lwd )
}


plot(cvallk(fit))


plot(cvallk(fit, prederr = function(y, ypred) (y - ypred)^4))
plot(cvallk(fit, prederr = function(y, ypred) abs(y - ypred)))
plot(cvallk(fit, prederr = function(y, ypred) exp(abs(y - ypred))))

head(mtcars)
dd <- mtcars


fit <- glm(I(mpg > 20) ~ wt + disp, mtcars, family = binomial)
summary(fit)
debug(cvk)
cvk(fit, K =6)

plot(cvallk(fit, K = 10:32))




