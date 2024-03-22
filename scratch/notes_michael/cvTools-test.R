#' ---
#' title: cvTools plotting methods
#' output: html_document
#' ---
#'

#' As guidance for what kinds of display methods might be useful for the `cv` package, here are
#' some examples from the `cvTools` package.
#'
#' I'm doing this to compare what is available in what cvTools returns in relation to graphs
#' vs. what we have already in `cv`.
#'

#+ echo=FALSE
knitr::opts_chunk$set(warning=FALSE,
                      message=FALSE,
                      R.options=list(digits=4))

#' This is from `example(cvTools::plot.cv)`

library(cvTools)
library("robustbase")
data("coleman")
set.seed(1234)  # set seed for reproducibility

# set up folds for cross-validation, K=5 folds, R=10 reps
folds <- cvFolds(nrow(coleman), K = 5, R = 10)


#' ## compare LS, MM and LTS regression

#' perform cross-validation for an LS regression model
fitLm <- lm(Y ~ ., data = coleman)
cvFitLm <- cvLm(fitLm, cost = rtmspe,
                folds = folds, trim = 0.1)

#' perform cross-validation for an MM regression model
fitLmrob <- lmrob(Y ~ ., data = coleman, k.max = 500)
cvFitLmrob <- cvLmrob(fitLmrob, cost = rtmspe,
                      folds = folds, trim = 0.1)

#' perform cross-validation for an LTS regression model
fitLts <- ltsReg(Y ~ ., data = coleman)
cvFitLts <- cvLts(fitLts, cost = rtmspe,
                  folds = folds,  trim = 0.1)

#' combine results into one object
cvFits <- cvSelect(LS = cvFitLm, MM = cvFitLmrob, LTS = cvFitLts)
cvFits

#' what does it contain?
names(cvFits)

cvFits$cv
cvFits$se

cvFits$reps


#' ## plot results for the MM regression model
#' Of these the density plot is more useful
plot(cvFitLmrob, method = "bw")
plot(cvFitLmrob, method = "density")

#' ## plot combined results
plot(cvFits, method = "bw")
plot(cvFits, method = "density")
plot(cvFits, method = "xy")
plot(cvFits, method = "dot")

#' ## Try to better with ggplot, ggdist

library(dplyr)
library(ggplot2)
library(ggdist)

reps <- cvFits$reps |>
  rename(method = Fit)

ggplot(reps, aes(x = CV, color=method, fill=method)) +
  geom_density(alpha = 0.3) +
  geom_rug() +
  labs(x = "Root mean squared prediction error (rtmsepe)") +
  theme_bw(base_size = 14)


#' ## Raincloud-ish plot

ggplot(reps, aes(x = method, y = CV, fill = method)) +
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.5,
    # move to the right
    justification = -0.2,
    alpha = 0.4,
    .width = c(0.68, 0.9, 0.95),
    linewidth = 3,
    # point_colour = NA
  ) +
  geom_boxplot(
    width = 0.12,
    alpha = 0.3,
    outlier.size = 3) +
  geom_jitter(aes(color = method), width = 0.1) +
  labs(y = "Root mean squared prediction error (rtmsepe)",
       x = "Estimation method") +
  theme_bw(base_size = 14)


