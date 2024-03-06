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
  # eval = nzchar(Sys.getenv("REBUILD_VIGNETTES"))
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

## ----AUCcomp------------------------------------------------------------------
AUCcomp <- function(y, yhat) 1 - Metrics::auc(y, yhat)

## ----Mroz-regression----------------------------------------------------------
data("Mroz", package="carData")
m.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
summary(m.mroz)

AUCcomp(with(Mroz, as.numeric(lfp == "yes")), fitted(m.mroz))

## ----Mroz-CV-ROC--------------------------------------------------------------
library("cv")
cv(m.mroz, criterion=AUCcomp, seed=3639)

## -----------------------------------------------------------------------------
mse

cv:::getLossFn(mse(rnorm(100), rnorm(100)))

## ----BEPS-data----------------------------------------------------------------
data("BEPS", package="carData")
head(BEPS)

## ----BEPS-model---------------------------------------------------------------
library("nnet")
m.beps <- multinom(vote ~ age + gender + economic.cond.national +
                       economic.cond.household + Blair + Hague + Kennedy +
                       Europe*political.knowledge, data=BEPS)

car::Anova(m.beps)

## ----BEPS-plot, fig.width=9, fig.height=5-------------------------------------
plot(effects::Effect(c("Europe", "political.knowledge"), m.beps,
            xlevels=list(Europe=1:11, political.knowledge=0:3),
            fixed.predictors=list(given.values=c(gendermale=0.5))),
     lines=list(col=c("blue", "red", "orange")),
     axes=list(x=list(rug=FALSE), y=list(style="stacked")))

## ----BayesRuleMulti-----------------------------------------------------------
head(BEPS$vote)
yhat <- predict(m.beps, type="class")
head(yhat)

BayesRuleMulti <- function(y, yhat){
  result <- mean(y != yhat)
  attr(result, "casewise loss") <- "y != yhat"
  result
}

BayesRuleMulti(BEPS$vote, yhat)

## ----BEPS-response-distribution-----------------------------------------------
xtabs(~ vote, data=BEPS)/nrow(BEPS)

## ----BEPS-test-default, error=TRUE--------------------------------------------
cv(m.beps, seed=3465, criterion=BayesRuleMulti)

## ----GetResponse.multinom-----------------------------------------------------
GetResponse.multinom <- function(model, ...) {
  insight::get_response(model)
}

head(GetResponse(m.beps))

## ----BEPS-test-default-2, error=TRUE------------------------------------------
cv(m.beps, seed=3465, criterion=BayesRuleMulti)

## ----cv.nultinom--------------------------------------------------------------
cv.multinom <- function (model, data, criterion=BayesRuleMulti, k, reps,
                         seed, ...){
  NextMethod(type="class", criterion=criterion,
             criterion.name=deparse(substitute(criterion)))
}

## ----BEPS-cv------------------------------------------------------------------
m.beps <- update(m.beps, trace=FALSE)
cv(m.beps, seed=3465)

## ----cv.lme-------------------------------------------------------------------
cv:::cv.lme

## ----GetResponse.glmmPQL------------------------------------------------------
GetResponse.glmmPQL <- function(model, ...){
  f <- formula(model)
  f[[3]] <- 1 # regression constant only on RHS
  model <- suppressWarnings(glm(f, data=model$data, family=model$family,
                                control=list(maxit=1)))
  cv::GetResponse(model)
}

## ----cv.glmmPQL---------------------------------------------------------------
cv.glmmPQL <- function(model, data = model$data, criterion = mse,
                     k, reps = 1, seed, ncores = 1, clusterVariables, 
                     blups=coef, fixed.effects=nlme::fixef, ...){
  cvMixed(
    model,
    package="MASS",
    data=data,
    criterion=criterion,
    k=k,
    reps=reps,
    seed=seed,
    ncores=ncores,
    clusterVariables=clusterVariables,
    predict.clusters.args=list(object=model,
                               newdata=data,
                               level=0,
                               type="response"),
    predict.cases.args=list(object=model,
                            newdata=data,
                            level=1,
                            type="response"),
    blups=blups,
    fixed.effects=fixed.effects,
    verbose=FALSE,
    ...)
}

## ----glmmPQL-example----------------------------------------------------------
library("MASS")
m.pql <- glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
             family = binomial, data = bacteria)
summary(m.pql)

## ----compare-to-lme4----------------------------------------------------------
library("lme4")
m.glmer <- glmer(y ~ trt + I(week > 2) + (1 | ID),
               family = binomial, data = bacteria)
summary(m.glmer)

  # comparison of fixed effects:
car::compareCoefs(m.pql, m.glmer) 

## ----swiss--------------------------------------------------------------------
library("leaps")
head(swiss)
nrow(swiss)

## ----swiss-lm-----------------------------------------------------------------
m.swiss <- lm(Fertility ~ ., data=swiss)
summary(m.swiss)

cv(m.swiss, seed=8433)

## ----subset-selection---------------------------------------------------------
swiss.sub <- regsubsets(Fertility ~ ., data=swiss)
summary(swiss.sub)

(bics <- summary(swiss.sub)$bic)
which.min(bics)

car::subsets(swiss.sub, legend="topright")

## ----best-model---------------------------------------------------------------
m.best <- update(m.swiss, . ~ . - Examination)
summary(m.best)

cv(m.best, seed=8433) # use same folds as before

## ----selectSubsets------------------------------------------------------------
selectSubsets <- function(data=insight::get_data(model), 
                          model,
                          indices,
                          criterion=mse,
                          details=TRUE, 
                          seed,
                          save.model=FALSE, 
                          ...){
  
  if (inherits(model, "lm", which=TRUE) != 1)
    stop("selectSubsets is appropriate only for 'lm' models")
  
  y <- GetResponse(model)
  formula <- formula(model)
  X <- model.matrix(model)

  if (missing(indices)) {
    if (missing(seed) || is.null(seed)) seed <- sample(1e6, 1L)
    # select the best model from the full data by BIC
    sel <- leaps::regsubsets(formula, data=data, ...)
    bics <- summary(sel)$bic
    best <- coef(sel, 1:length(bics))[[which.min(bics)]]
    x.names <- names(best)
    # fit the best model; intercept is already in X, hence - 1:
    m.best <- lm(y ~ X[, x.names] - 1) 
    fit.all <- predict(m.best, newdata=data)
    return(list(
      criterion = criterion(y, fit.all),
      model = if (save.model) m.best else NULL # return best model
      ))
  }

  # select the best model omitting the i-th fold (given by indices)
  sel.i <- leaps::regsubsets(formula, data[-indices, ], ...)
  bics.i <- summary(sel.i)$bic
  best.i <- coef(sel.i, 1:length(bics.i))[[which.min(bics.i)]]
  x.names.i <- names(best.i)
  m.best.i <- lm(y[-indices] ~ X[-indices, x.names.i] - 1)
              # predict() doesn't work here:
  fit.all.i <- as.vector(X[, x.names.i] %*% coef(m.best.i))
  fit.i <- fit.all.i[indices]
  # return the fitted values for i-th fold, CV criterion for all cases, 
  #   and the regression coefficients
  list(fit.i=fit.i, # fitted values for i-th fold
       crit.all.i=criterion(y, fit.all.i), # CV crit for all cases
       coefficients = if (details){ # regression coefficients
         coefs <- coef(m.best.i)
         
         # fix coefficient names
         names(coefs) <- sub("X\\[-indices, x.names.i\\]", "",
                             names(coefs))
         
         coefs
       }  else {
         NULL
       }
  )
}

## ----test-selectSubsets-------------------------------------------------------
selectSubsets(model=m.swiss)

## ----test-selectSubsets-fold--------------------------------------------------
selectSubsets(model=m.swiss, indices=seq(5, 45, by=10))

## ----cvSelect-selectSubsets---------------------------------------------------
(cv.swiss <- cvSelect(selectSubsets, model=m.swiss,
                      data=swiss, seed=8433)) # use same folds

## ----best-models-by-folds-----------------------------------------------------
compareFolds(cv.swiss)

## -----------------------------------------------------------------------------
AUC <- function(y, yhat = seq_along(y)) {
  s <- sum(y)
  if (s == 0) return(0)
  if (s == length(y)) return(1)
  Metrics::auc(y, yhat)
}

## -----------------------------------------------------------------------------
Ymat <- function(n_v, exclude_identical = FALSE) {
  stopifnot(n_v > 0 && round(n_v) == n_v)    # n_v must be a positive integer
  ret <- sapply(0:(2^n_v - 1),
                function(x) as.integer(intToBits(x)) )[1:n_v, ]
  ret <- if (is.matrix(ret)) t(ret) else matrix(ret)
  colnames(ret) <- paste0("y", 1:ncol(ret))
  if (exclude_identical) ret[-c(1, nrow(ret)), ] else ret
}

## -----------------------------------------------------------------------------
Ymat(3)

## -----------------------------------------------------------------------------
Ymat(3, exclude_identical = TRUE)

## -----------------------------------------------------------------------------
cbind(Ymat(3), AUC = apply(Ymat(3), 1, AUC))

## -----------------------------------------------------------------------------
resids <- function(n_v, exclude_identical = FALSE, 
                   tol = sqrt(.Machine$double.eps)) {
  Y <- Ymat(n_v, exclude_identical = exclude_identical)
  AUC <- apply(Y, 1, AUC)
  X <- cbind(1-Y, Y)
  opts <- options(warn = -1)
  on.exit(options(opts))
  fit <- lsfit(X, AUC, intercept = FALSE)
  ret <- max(abs(residuals(fit)))
  if(ret < tol){
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

