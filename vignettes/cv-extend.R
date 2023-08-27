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
  comment = "#>"
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
  mean(y != yhat)
}

BayesRuleMulti(BEPS$vote, yhat)

## ----BEPS-response-distribution-----------------------------------------------
xtabs(~ vote, data=BEPS)/nrow(BEPS)

## ----BEPS-test-default, error=TRUE--------------------------------------------
cv(m.beps, seed=3465, criterion=BayesRuleMulti)

## ----getRespons.multinom------------------------------------------------------
getResponse.multinom <- function(model, ...) {
  insight::get_response(model)
}

head(getResponse(m.beps))

## ----BEPS-test-default-2, error=TRUE------------------------------------------
cv(m.beps, seed=3465, criterion=BayesRuleMulti)

## ----cv.nultinom--------------------------------------------------------------
cv.multinom <- function (model, data, criterion=BayesRuleMulti, k, reps,
                         seed, ...){
  NextMethod(type="class", criterion=criterion)
}

## ----BEPS-cv------------------------------------------------------------------
m.beps <- update(m.beps, trace=FALSE)
cv(m.beps, seed=3465)

