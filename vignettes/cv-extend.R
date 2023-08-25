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

