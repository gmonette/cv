#'
#' Experimenting with CV bias and variance
#'

{
  # Initialize

  library(glmnet)
  library(cv)
  library(insight)
  library(latticeExtra)
  library(spida2)
  library(car)

  getResponse <- insight::get_response
  spida2::setwd_here()
}


{
  #
  # Supernova data from Efron and Hastie (2013) Computer Age Statistical Inference
  #

  if(!file.exists('supernova.csv')){
    # import data
    supernova.in <- read.table("http://hastie.su.domains/CASI_files/DATA/supernova.txt", header=TRUE)
    write.csv(supernova.in, file = 'supernova.csv')
  }
  supernova <- read.csv('supernova.csv', row.names = 1)
}

fit1 <- lm(Magnitude ~ E1, supernova)
fit2 <- lm(Magnitude ~ E1 + E2, supernova)
full <- lm(Magnitude ~ ., supernova)    # E1 to E10

anova(fit1, fit2, full)

summary(fit1)
summary(fit2)
summary(full)

xyplot(Magnitude~predict(full), supernova, pch = 16) +
  layer(panel.dell(...)) +
  layer(panel.lmline(..., lty = 1)) +
  layer(panel.abline(a=0,b=1, lty = 2, lwd = 3, col = 'red'))

cvk <- cv(models(fit1, fit2, full), data = supernova, k = 10)
cvloo <- cv(models(fit1, fit2, full), data = supernova, k = 'loo')
cvkreps <- cv(models(fit1, fit2, full), data = supernova, k = 10, reps = 100)
plot(cvk)
plot(cvloo)
plot(cvkreps)


{

  #
  # Simulate N samples using the supernova data
  #
  N <- 1000
  sims <- lapply(1:N, \(i) supernova[sample(1:39, replace = T),])

  mods <- lapply(
    sims,
    function(d) {
      list(
        fit1 = lm(Magnitude ~ E1, d),
        fit2 = lm(Magnitude ~ E1 + E2, d),
        full = lm(Magnitude ~ ., d)
      )
    }
  )

  kfold.l <- lapply(mods, lapply, cv)

  loo.l <- lapply(mods, lapply, cv, k = 'loo')

  kfoldrep.l <- lapply(mods, lapply,
                       function(mod) {
                         cv:::summarizeReps(cv(mod, k = 10, reps = 10))
                       }
  )

  kfoldcomp.l <- lapply(
    seq_along(mods),
    function(i){
      cv(do.call(models, mods[[i]]), k = 10, data = sims[[i]])
    }
  )

  #
  # Out-of-sample squared-error loss on finite population
  #
  Err_XY <- lapply(
    mods,
    lapply,
    function(mod) sum((supernova$Magnitude - predict(mod, newdata = supernova))^2)/nrow(supernova))

  popfull <- lm(Magnitude ~ ., supernova)
  popfit1 <- lm(Magnitude ~ E1, supernova)
  popfit2 <- lm(Magnitude ~ E1 + E2, supernova)

  sig2full <- summary(popfull)$sigma^2
  sig2fit1 <- summary(popfit1)$sigma^2
  sig2fit2 <- summary(popfit2)$sigma^2


  summ <- within(
    list(),
    {

      kfold <-    unlist(lapply(kfold.l, lapply, `[[`, 'CV crit'))
      kfold_ba <- unlist(lapply(kfold.l, lapply, `[[`, 'adj CV crit'))
      loo <-      unlist(lapply(loo.l,   lapply, `[[`,  'CV crit'))
      kfoldrep <-      unlist(lapply(kfoldrep.l,   lapply, `[[`,  'CV crit'))
      kfoldrep_ba <-   unlist(lapply(kfoldrep.l,   lapply, `[[`,  'adj CV crit'))
      kfoldcomp <-      unlist(lapply(kfoldcomp.l,   lapply, `[[`,  'CV crit'))
      kfoldcomp_ba <-   unlist(lapply(kfoldcomp.l,   lapply, `[[`,  'adj CV crit'))

      Err_XY <- unlist(Err_XY)

      model <- rep(c('fit1','fit2','full'), length(sims))
      sim <- rep(1:length(sims), each = 3)

      sig2 <- rep(c(sig2fit1, sig2fit2, sig2full), length(sims))
      n_sampled <- rep( sapply(sims, \(d) length(unique(d$Magnitude))), each = 3)


      # model comparisons

      kfold_mc <- capply(kfold, sim, \(v) v - v[3])
      loo_mc <- capply(loo, sim, \(v) v - v[3])
      kfoldrep_mc <- capply(kfoldrep, sim, \(v) v - v[3])
      kfoldcomp_mc <- capply(kfoldcomp, sim, \(v) v - v[3])

    }
  )

  sapply(summ, length)
  summ <- as.data.frame(rev(summ))

  summ$Err       <- with(summ, capply(Err_XY, model, mean))
  summ$kfold_var <- with(summ, capply(kfold, model, var))
  summ$loo_var   <- with(summ, capply(loo, model, var))

  xyplot(kfold_mc ~ factor(model), summ, groups = sim, type = 'l')
  xyplot(kfold ~ factor(model), summ, groups = sim, type = 'l')
  tab(summ, ~ model)

  xyplot(Err_XY ~ kfold | model, summ)
  xyplot(Err_XY ~ loo + kfold + kfoldrep + kfoldcomp, summ, subset = model == 'full', outer = T)
  xyplot(Err_XY ~ loo | model, summ)
  scatterplotMatrix(~Err_XY + kfold + loo + kfoldrep + kfoldcomp, data = subset(summ, model == 'full'))

  ##
  ## Re BHT Corollary 3: E (Errhat - Err_XY)^2 - E( Errhat - Err)^2 ~ Omega(1/n)
  ##
  ## BHT p. 6:
  ##
  ## If Err_hat and Err_XY are asymptotically uncorrelated (not in this example):
  ##
  ## MSE( Err_hat - Err_XY) = {E(Err_hat - Err)}^2 + var(Err_hat) + var(Err_XY)
  ##                              bias^2
  ##

  ms <- function(x) sum(x^2)/length(x)

  # LHS: MSE( Err_hat - Err_XY)

  (MSE_kfold_XY <- with(summ, by(summ, ~ model, \(d) with(d, ms(kfold -  Err_XY )) )))
  (MSE_kfold_ba_XY <- with(summ, by(summ, ~ model, \(d) with(d, ms(kfold_ba -  Err_XY )) )))
  (MSE_loo_XY <- with(summ, by(summ, ~ model, \(d) with(d, ms(loo -  Err_XY ) ))))
  (MSE_kfoldrep_XY <- with(summ, by(summ, ~ model, \(d) with(d, ms(kfoldrep -  Err_XY ) ))))

  # Bias

  (bias2_kfold <- with(summ, by(summ, ~ model, \(d) with(d, mean(kfold -  Err ) )))^2)
  (bias2_kfold_ba <- with(summ, by(summ, ~ model, \(d) with(d, mean(kfold_ba -  Err ) )))^2)
  (bias2_loo <- with(summ, by(summ, ~ model, \(d) with(d, mean(loo -  Err ) )))^2)
  (bias2_kfoldrep <- with(summ, by(summ, ~ model, \(d) with(d, mean(kfoldrep -  Err ) )))^2)

  # var(Err_hat)

  (var_kfold <- with(summ, by(summ, ~ model, \(d) with(d, var(kfold) ))))
  (var_kfold_ba <- with(summ, by(summ, ~ model, \(d) with(d, var(kfold_ba) ))))
  (var_loo <- with(summ, by(summ, ~ model, \(d) with(d, var(loo) ))))
  (var_kfoldrep <- with(summ, by(summ, ~ model, \(d) with(d, var(kfoldrep) ))))

  # var Err_XY
  (var_Err_XY <- with(summ, by(summ, ~ model, \(d) with(d, var(Err_XY ) ))))


  summ_kfold <- cbind(MSE_kfold_XY, bias2_kfold, var_kfold, var_Err_XY)
  summ_kfold_ba <- cbind(MSE_kfold_ba_XY, bias2_kfold_ba, var_kfold_ba, var_Err_XY)
  summ_loo <- cbind(MSE_loo_XY, bias2_loo, var_loo, var_Err_XY)
  summ_kfoldrep <- cbind(MSE_kfoldrep_XY, bias2_kfoldrep, var_kfoldrep, var_Err_XY)


  round(summ_kfold,6)
  round(summ_kfold_ba,6)
  round(summ_loo,6)
  round(summ_kfoldrep,6)

  summ2 <- as.data.frame(rbind(summ_kfold, summ_kfold_ba, summ_loo, summ_kfoldrep))
  names(summ2) <- c('MSE_CV_XY', 'bias2_CV','var_CV', 'var_Err_XY')
  summ2$model <- factor(rep(c('fit1','fit2','full'), 4))
  summ2$CV <- factor(rep(c('kfold','kfold_ba','loo','kfoldrep'), each = 3))


  summ2 <- sortdf(summ2, ~CV)

  round(summ2,5)
  options(scipen=20, digits = 4)

  (xyplot(Err_XY ~ loo + kfold + kfold_ba + kfoldrep, summ, subset = model == 'full', outer = T) +
    layer(panel.spline(...)))

  xyplot(var_CV + bias2_CV ~ CV, summ2, groups = model, type = 'b',
         auto.key = T, scales = list(y=list(log=T, relation = 'free'))) %>% print


  scatterplotMatrix(~ Err_XY + loo + kfold + kfold_ba, data = subset(summ, model == 'full'), regline = TRUE) %>% print
  scatterplotMatrix(~ Err_XY + loo + kfold + kfold_ba, data = subset(summ, model == 'fit2'), regline = TRUE) %>% print


} # end of simulation



#'
#' BH
#'
#'



library(p3d)
library(viridis)
Init3d(cex = 2)
ns <- length(unique(summ$n_sampled))

Err_full <- mean(subset(summ, model == 'full')$Err_XY)

cols <- viridis_pal()(ns)
Plot3d(Err_XY ~ kfold + loo | n_sampled, subset(summ, model == 'full'), col = cols)
Plot3d(Err_XY ~ loo + n_sampled, subset(summ, model == 'full'), col = cols)
Plot3d(Err_XY ~ kfold + loo | n_sampled, subset(summ, model == 'full'), col = cols)
Plot3d(log(Err_XY) ~ kfold + n_sampled , subset(summ, model == 'full'), col = cols)

# Fit3d(function(...) Err_full)

Axes3d()
Plot3d( log(Err_XY) ~ n_sampled+ loo , subset(summ, model == 'full'), col = heat.colors(6))


fit1 <- lm(log(Err_XY) ~ loo + I(loo^2),  subset(summ, model == 'full'))

fit <- lm( log(Err_XY) ~ loo * n_sampled + I(loo^2) + I(n_sampled^2) , subset(summ, model == 'full'))
summary(fit)
Fit3d(fit)
Fit3d(fit1)



Err <- mean(r$Err_XY)
Fit3d

f <- function(loo, kfold) Err
Fit3d(f)

# Considerations:
#
# To get Err_X by simulation, we need to be able to generate Y conditional on X
#
# - Note in passing that although X and Y are generated with 0 marginal means
#   the models have intercepts, so we don't assume 0 marginal means in models
#
#
#
#
#
# 1000 simulations from a normal with the same covariance as the
# supernova data
#
# In order to estimate Err_X by simulation we need to simulate
# many Ys for the same X which can be done economically
# using Z matrices and a Choleski variance factor
#

pop_cov <- var(as.matrix(supernova))
eigen(pop_cov)$values
rchol <- chol(pop_cov)
max(abs(pop_cov - t(rchol) %*% rchol))

max(abs(pop_cov - crossprod(rchol)))

sim <- function(R, n, ny = 0) {
  # R is the right choleski factor of the variance of Xs and Y
  #   whose marginal mean is assumed to be 0
  # n is the sample size
  # ny is the number of additional conditionally independent columns of Ys
  #
  d <- ncol(R)
  Z <- matrix(rnorm(d * n),n)
  Data <- Z %*% R      # X matrix and Y vector
  EY <- Z[,-d] %*% R[-d, d]    # conditional expectation of Y
  Ymat <- matrix(rnorm(n * ny), n) * c(R[d,d]) + c(EY)    # matrix of additional Ys
  dimnames(Data) <- NULL
  dimnames(Ymat) <- NULL
  ret <- data.frame(Y = Data[,d])
  ret$X <- Data[,-d]
  ret$Ymat <- Ymat
  ret
}

sims$X


ret$X
sims <- sim(rchol, 39, 3)
sims %>% names
lm(Y ~ X, sims)
lm(Ymat ~ X, sims)


sims$X
lapply(sims, \(d) {
  lm(Magnitude ~ )
}
sims$X

condvy <- pop_cov[11,11] - pop_cov[11,-11]%*% solve( pop_cov[1:10,1:10], pop_cov[-11,11])

# should be equal
condvy
rchol[11,11]^2


popfit <- lm(Magnitude ~ ., as.data.frame(matrix(rnorm(10000000*11), 10000000) %*% rchol))
summary(popfit)

truey <- function(d) {
  ypred <- predict(popfit, newdata = d)
  ypred + rnorm(nrow(d)) * summary(popfit)$sigma
}


sims <- lapply(1:1000,
               function(i) {
                 z <- matrix(rnorm(39*11), 39)
                 rcholp <- rchol[,39]
                 rcholp[39] <- 0
                 ey <- z %*% rcholp
                 data <- z %*% rchol
                 list( data = as.data.frame(data), ey = ey, )
                   mat <-
                   data = as.data.frame(matrix(rnorm(39*11), 39) %*% rchol
               })
sims[[1]]
modlist <- lapply(sims, \(d) lm(Magnitude ~ ., d))



cvlistkfold <- lapply(modlist, cv)
cvlistloo <- lapply(modlist, cv, k = 'loo')

Err_XY <- lapply(seq_along(modlist),
                   \(i) sum((predict(modlist[[i]]) - truey(sims[[i]]))^2)/nrow(sims[[i]]))
Err_X <- lapply(seq_along(modlist),
                 \(i) sum((predict(modlist[[i]]) - truey(sims[[i]]))^2)/nrow(sims[[i]]))

r <- within(list(),
            {
              kfold <- sapply(cvlistkfold, `[[`, 'CV crit')
              kfold_ba <- sapply(cvlistkfold, `[[`  ,'adj CV crit')
              loo <- sapply(cvlistloo, `[[`, 'CV crit')
              Err_XY <- unlist(Err_XY)
            }
)
r <- as.data.frame(r)

r <- cbind(r, )
#
# Estimating Err_X as defined in BHT
#



library(p3d)
Init3d(cex = 2)
Plot3d( Err_XY ~ kfold + loo , r)
Axes3d()


summary(lm(Err_XY ~ loo, r))
summary(lm(Err_XY ~ kfold, r))

sapply(r, mean)






results %>% lapply(length)
results %>% lapply(mode)
which.max(results$Err_XY)

sims[[312]]

# number of unique values caught




plomodlist <- do.call(models,modlist)
cvlist <- cv(modlist)
?sample
methods(cv)


