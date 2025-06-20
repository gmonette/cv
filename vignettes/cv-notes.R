## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = TRUE,
  warning = TRUE,
  fig.align = "center",
  fig.height = 6,
  fig.width = 7,
  fig.path = "fig/",
  dev = "png",
  comment = "#>" #,
)

# save some typing
knitr::set_alias(w = "fig.width",
                 h = "fig.height",
                 cap = "fig.cap")

# colorize text: use inline as `r colorize(text, color)`
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
      x)
  } else x
}


.opts <- options(digits = 5)

## ----Auto-data-quaadratic-regression------------------------------------------
data("Auto", package="ISLR2")
m.auto <- lm(mpg ~ poly(horsepower, 2), data = Auto)
summary(m.auto)

library("cv")
summary(cv(m.auto, k = "loo") ) # default method = "hatvalues"
summary(cv(m.auto, k = "loo", method = "naive"))
summary(cv(m.auto, k = "loo", method = "Woodbury"))

## ----cv.lm.timings, cache=TRUE------------------------------------------------
m.auto.glm <- glm(mpg ~ poly(horsepower, 2), data = Auto)
boot::cv.glm(Auto, m.auto.glm)$delta

microbenchmark::microbenchmark(
  hatvalues = cv(m.auto, k = "loo"),
  Woodbury = cv(m.auto, k = "loo", method = "Woodbury"),
  naive = cv(m.auto, k = "loo", method = "naive"),
  cv.glm = boot::cv.glm(Auto, m.auto.glm),
  times = 10
)

## ----Mroz-logistic-regression-------------------------------------------------
data("Mroz", package="carData")
m.mroz <- glm(lfp ~ ., data = Mroz, family = binomial)
summary(m.mroz)

summary(cv(m.mroz, # default method = "exact"
   k = "loo", 
   criterion = BayesRule))
summary(cv(m.mroz,
   k = "loo",
   criterion = BayesRule,
   method = "Woodbury"))
summary(cv(m.mroz,
   k = "loo",
   criterion = BayesRule,
   method = "hatvalues"))

## ----glm.timings, cache=TRUE--------------------------------------------------
microbenchmark::microbenchmark(
  hatvalues = cv(
    m.mroz,
    k = "loo",
    criterion = BayesRule,
    method = "hatvalues"
  ),
  Woodbury = cv(
    m.mroz,
    k = "loo",
    criterion = BayesRule,
    method = "Woodbury"
  ),
  exact = cv(m.mroz, k = "loo", criterion = BayesRule),
  cv.glm = boot::cv.glm(Mroz, m.mroz,
                        cost = BayesRule),
  times = 10
)

## -----------------------------------------------------------------------------
AUC <- function(y, yhat = seq_along(y)) {
  s <- sum(y)
  if (s == 0)
    return(0)
  if (s == length(y))
    return(1)
  Metrics::auc(y, yhat)
}

## -----------------------------------------------------------------------------
Ymat <- function(n_v, exclude_identical = FALSE) {
  stopifnot(n_v > 0 &&
              round(n_v) == n_v)    # n_v must be a positive integer
  ret <- sapply(0:(2 ^ n_v - 1),
                function(x)
                  as.integer(intToBits(x)))[1:n_v,]
  ret <- if (is.matrix(ret))
    t(ret)
  else
    matrix(ret)
  colnames(ret) <- paste0("y", 1:ncol(ret))
  if (exclude_identical)
    ret[-c(1, nrow(ret)),]
  else
    ret
}

## -----------------------------------------------------------------------------
Ymat(3)

## -----------------------------------------------------------------------------
Ymat(3, exclude_identical = TRUE)

## -----------------------------------------------------------------------------
cbind(Ymat(3), AUC = apply(Ymat(3), 1, AUC))

## -----------------------------------------------------------------------------
resids <- function(n_v,
                   exclude_identical = FALSE,
                   tol = sqrt(.Machine$double.eps)) {
  Y <- Ymat(n_v, exclude_identical = exclude_identical)
  AUC <- apply(Y, 1, AUC)
  X <- cbind(1 - Y, Y)
  opts <- options(warn = -1)
  on.exit(options(opts))
  fit <- lsfit(X, AUC, intercept = FALSE)
  ret <- max(abs(residuals(fit)))
  if (ret < tol) {
    ret <- 0
    solution <- coef(fit)
    names(solution) <- paste0("c(", c(1:n_v, 1:n_v), ",",
                              rep(0:1, each = n_v), ")")
    attr(ret, "solution") <- zapsmall(solution)
  }
  ret
}

## -----------------------------------------------------------------------------
resids(3, exclude_identical = TRUE)

## -----------------------------------------------------------------------------
resids(3, exclude_identical = FALSE)

## -----------------------------------------------------------------------------
resids(4, exclude_identical = TRUE)
resids(4, exclude_identical = FALSE)

## ----coda, include = FALSE----------------------------------------------------
options(.opts)

