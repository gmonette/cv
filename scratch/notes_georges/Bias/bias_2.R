#'
#' Experimenting with bias and variance
#'

spida2::setwd_here()

if(FALSE){
  supernova.in <- read.table("http://hastie.su.domains/CASI_files/DATA/supernova.txt", header=TRUE)
  supernova.in
  write.csv(supernova.in, file = 'supernova.csv')
}
supernova <- read.csv('supernova.csv', row.names = 1)

# source nested CV

library(glmnet)
library(cv)
library(insight)
library(latticeExtra)
library(spida2)

getResponse <- insight::get_response


fit <- lm(Magnitude ~ ., supernova)
summary(fit)

fit2 <- lm(Magnitude ~ E1+ E2, supernova)
summary(fit2)

xyplot(Magnitude~predict(fit), supernova, pch = 16) +
  layer(panel.dell(...)) +
  layer(panel.lmline(..., lty = 1)) +
  layer(panel.abline(a=0,b=1, lty = 2, lwd = 3, col = 'red'))


ocv <- cv(fit)
ncv <- cv(fit, reps = 10)

models(fit, fit2) %>% cv(k=10, data = supernova) %>% plot




?cv
plot(ocv)
class(ncv)
ncv
cv(fit2)
cv(models(full=fit, rest=fit2), k=2)
?models

cv(list(fit, fit2), data = supernova)
fit2
methods(class = 'cvList')
cv(fit, reps = 10)
%>% as.data.frame -> zncv

zncv


methods(class='cv')
z
?cv
ls(3)
search()
?cvSelect


#
# 1000 simulations using the supernova data
#
sims <- lapply(1:1000,
               function(i) {
                 supernova[sample(1:39, replace = T),]
               })


modlist <- lapply(sims, function(d) lm(Magnitude ~ ., d))

cvlistkfold <- lapply(modlist, cv)
cvlistloo <- lapply(modlist, cv, k = 'loo')

#
# Out-of-sample squared-error loss on finite population
#
Err_XY <- lapply(modlist, function(mod) sum((supernova$Magnitude - predict(mod, newdata = supernova))^2)/nrow(supernova))


r <- within(list(),
  {
  kfold <- sapply(cvlistkfold, `[[`,      'CV crit')
  kfold_ba <- sapply(cvlistkfold, `[[`  , 'adj CV crit')
  loo <- sapply(cvlistloo, `[[`,          'CV crit')

  sample_size <- sapply(sims, \(d) length(unique(d$Magnitude)))  # number of elements in sample
  # uniq_cat <- local(
  #   {
  #     breaks <- c(-Inf, sort(uniq_num)[c(5,20,32)], Inf)
  #     cut(uniq_num, breaks, c(names(breaks)[-1]))
  #   }
  #   )
  Err_XY <- unlist(Err_XY)
  }
)
r <- as.data.frame(r)

library(p3d)
Init3d(cex = 2)
Plot3d( Err_XY ~ kfold + loo | sample_size, r, col = heat.colors(length(unique(r$uniq_cat))))
Axes3d()
Plot3d( Err_XY ~ uniq_num + loo , r, col = heat.colors(6))

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

simdata <- function(R, n, ny = 0) {
  # R is the right choleski factor of the variance of Xs and Y
  #   whose marginal mean is assumed to be 0
  # n is the sample size
  # ny is the number of additional conditionally independent columns of Ys
  #

  d <- ncol(R)
  Z <- matrix(rnorm(d * n),n)
  D <- Z %*% R      # X matrix and Y vector
  EY <- Z[,-d] %*% R[-d, d]    # conditional expectation of Y
  Ymat <- matrix(rnorm(n * ny), n) * c(R[d,d]) + c(EY)    # matrix of additional Ys
  ret <- data.frame(Y = D[,d])
  ret$X <- D[,-d]
  ret$Ymat <- Ymat
  ret
}

simanal <- function(R, n, ny = 100) {

  # R: right choleski factor of variance of X and Y
  # n: sample size
  # ny: number of Ys to simulate for Err_X simulated estimate

  dat <- simdata(R = R, n = n, ny = ny)
  fit <- lm(Y ~ X, data = dat)
  fitmat <- lm(Ymat ~ X, data = dat)
  cvk <- cv(fit, k = 10)
  cvloo <- cv(fit, k = 'loo')
  cvnested <- cv:::summarizeReps(cv(fit, reps = 10))

  # Err_X: See Bias_notes.xopp

  theta <- coef(fit)
  Err_XY <- theta[1]^2 + sum((R %*% c(theta[-1],-1) )^2)

  # simulated estimate of Err_X

  Theta <- coef(fitmat)
  Err_X <- (sum(Theta[1,]^2) + sum((R %*% rbind(Theta[-1,],-1))^2))/ny

  ret <- within(
    list(),
    {
                kfold <- cvk[['CV crit']]
                kfold_ba <- cvk[['adj CV crit']]
                loo <- cvk[['CV crit']]
                loo_ba <- cvk[['adj CV crit']]
                nested <- cvnested[["CV crit"]]
                nested_ba <- cvnested[["adj CV crit"]]
                nested_sd <- cvnested[["SD CV crit"]]
                nested_ba_sd <- cvnested[["SD adj CV crit"]]
                Err_X <- Err_X
                Err_XY <- as.vector(Err_XY)
                sigma2 <- R[nrow(R),ncol(R)]^2
    })
  rev(unlist(ret))

}







system.time(
simanal(rchol, 30, 1000)

)

sims <- replicate(1000, simanal(rchol,30,500))
sims <- as.data.frame(t(sims))

library(p3d)
Init3d(cex = 1.5)

head(sims)
Plot3d(Err_X ~ nested + Err_XY, sims, sphere = 2, surface =T, residuals = F)


sims <- sim(rchol, 20, 4)



sims$Ymat


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


