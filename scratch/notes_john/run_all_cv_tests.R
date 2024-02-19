run_all_cv_tests <- function(run=TRUE){
  if (run) Sys.setenv(RUN_ALL_CV_TESTS = "true")
  else Sys.setenv(RUN_ALL_CV_TESTS = "")
}

run_all_cv_tests()
