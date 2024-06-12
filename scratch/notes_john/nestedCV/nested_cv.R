nested_cv <-
function (X, Y, funcs, reps = 50, n_folds = 10, alpha = 0.1,
          bias_reps = NA, funcs_params = NULL, n_cores = 1, verbose = F)
{
  if (verbose) {
    t1 <- Sys.time()
    temp <- nested_cv_helper(X, Y, funcs, n_folds,
                                        funcs_params = funcs_params)
    t2 <- Sys.time()
    print(paste0("Estimated time required: ", (t2 - t1) *
                   reps))
  }
  var_pivots <- c()
  gp_errs <- c()
  ho_errs <- c()
  if (n_cores == 1) {
    raw <- lapply(1:reps, function(i) {
      nested_cv_helper(X, Y, funcs, n_folds,
                                  funcs_params = funcs_params)
    })
  }
  else {
    raw <- parallel::mclapply(1:reps, function(i) {
      nested_cv_helper(X, Y, funcs, n_folds,
                                  funcs_params = funcs_params)
    }, mc.cores = n_cores)
  }
  for (i in 1:reps) {
    temp <- raw[[i]]
    var_pivots <- rbind(var_pivots, temp$pivots)
    ho_errs <- c(ho_errs, temp$errs)
  }
  assign("ab0", cbind(var_pivots[, 1]^2, var_pivots[, 2]), envir=.GlobalEnv)
  n_sub <- floor(length(Y) * (n_folds - 1)/n_folds)
  ugp_infl <- sqrt(max(0, mean(var_pivots[, 1]^2 - var_pivots[,
                                                              2])))/(sd(ho_errs)/sqrt(n_sub))
  ugp_infl <- max(1, min(ugp_infl, sqrt(n_folds)))
  infl_est2 <- sqrt(pmax(0, sapply(1:reps, function(i) {
    mean(var_pivots[1:(i * n_folds), 1]^2 - var_pivots[1:(i *
                                                            n_folds), 2])
  })))/(sd(ho_errs)/sqrt(n_sub))
  cv_means <- c()
  bias_est <- 0
  if (is.na(bias_reps)) {
    bias_reps <- ceiling(reps/5)
  }
  if (bias_reps == 0) {
    bias_est <- 0
  }
  else {
    for (i in 1:bias_reps) {
      temp <- nestedcv:::naive_cv(X, Y, funcs, n_folds,
                                  funcs_params = funcs_params)
      cv_means <- c(cv_means, temp$err_hat)
    }
    bias_est <- (mean(ho_errs) - mean(cv_means)) * (1 + ((n_folds -
                                                            2)/(n_folds))^(1.5))
  }
  assign("mean.ho_errs", mean(ho_errs), envir=.GlobalEnv)
  assign("mean.cv_means", mean(cv_means), envir=.GlobalEnv)
  pred_est <- mean(ho_errs) - bias_est
  list(sd_infl = ugp_infl, err_hat = pred_est, ci_lo = pred_est -
         qnorm(1 - alpha/2) * sd(ho_errs)/sqrt(length(Y)) * ugp_infl,
       ci_hi = pred_est + qnorm(1 - alpha/2) * sd(ho_errs)/sqrt(length(Y)) *
         ugp_infl, raw_mean = mean(ho_errs), bias_est = bias_est,
       sd = sd(ho_errs), running_sd_infl = infl_est2)
}
