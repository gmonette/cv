nested_cv_helper <-
function (X, Y, funcs, n_folds = 10, funcs_params = NULL)
{
  fold_id <- 1:nrow(X)%%n_folds + 1
  fold_id <- sample(fold_id[1:(nrow(X)%/%n_folds * n_folds)])
  fold_id <- c(fold_id, rep(0, nrow(X)%%n_folds))
  ho_errors <- array(0, dim = c(n_folds, n_folds, nrow(X)%/%n_folds))
  for (f1 in 1:(n_folds - 1)) {
    for (f2 in (f1 + 1):n_folds) {
      test_idx <- c(which(fold_id == f1), which(fold_id ==
                                                  f2))
      fit <- funcs$fitter(X[-test_idx, ], Y[-test_idx],
                          funcs_params = funcs_params)
      preds <- funcs$predictor(fit, X, funcs_params = funcs_params)
      ho_errors[f1, f2, ] <- funcs$loss(preds[fold_id ==
                                                f1], Y[fold_id == f1], funcs_params = funcs_params)
      ho_errors[f2, f1, ] <- funcs$loss(preds[fold_id ==
                                                f2], Y[fold_id == f2], funcs_params = funcs_params)
    }
  }
  out_mat <- matrix(0, n_folds, 2)
  for (f1 in 1:(n_folds)) {
    test_idx <- which(fold_id == f1)
    fit <- funcs$fitter(X[-test_idx, ], Y[-test_idx], funcs_params = funcs_params)
    preds <- funcs$predictor(fit, X[test_idx, ], funcs_params = funcs_params)
    e_out <- funcs$loss(preds, Y[test_idx])
    e_bar_t <- c()
    for (f2 in 1:n_folds) {
      if (f2 == f1) {
        next
      }
      e_bar_t <- c(e_bar_t, ho_errors[f2, f1, ])
    }
    out_mat[f1, 1] <- mean(e_bar_t) - mean(e_out)
    out_mat[f1, 2] <- var(e_out)/length(test_idx)
  }
  all_ho_errs <- c()
  for (f1 in 1:(n_folds - 1)) {
    for (f2 in (f1 + 1):n_folds) {
      all_ho_errs <- c(all_ho_errs, ho_errors[f1, f2, ],
                       ho_errors[f2, f1, ])
    }
  }
  return(list(pivots = out_mat, errs = all_ho_errs))
}
