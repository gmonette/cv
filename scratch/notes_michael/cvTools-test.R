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

# set up folds for cross-validation, K=5 folds, R=25 reps
folds <- cvFolds(nrow(coleman), K = 5, R = 25)


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

cvFits$reps |> head()


#' ## plot results for the MM regression model
#' Of these the density plot is more useful
#+ out.width="50%", fig.show="hold"
plot(cvFitLmrob, method = "bw")
plot(cvFitLmrob, method = "density")

#' ## plot combined results
#+ out.width="50%", fig.show="hold"
plot(cvFits, method = "bw")
plot(cvFits, method = "density")

#+ out.width="50%", fig.show="hold"
plot(cvFits, method = "xy")
plot(cvFits, method = "dot")

#' ## Try to better with ggplot, ggdist

library(dplyr)
library(ggplot2)
library(ggdist)

#' Wrangle the `cvFits` to give data frames for plotting in a density plot. Need to get the data for
#' the overall CV measures in a separate data.frame
reps <- cvFits$reps |> rename(method = Fit)
cv <- cvFits$cv |> rename(method = Fit)
se <- cvFits$se |> rename(method = Fit, se = CV)
cv_all <- data.frame(cv, se = se$se,
                     y = 0.25 + .2* as.numeric(cv$method))  # scale relative to Density (y) axis

ggplot(reps, aes(x = CV, color=method, fill=method, shape=method)) +
  geom_density(alpha = 0.3) +
#  geom_rug() +
  geom_jitter(aes(y = -0.05), width=0, height = 0.05) +
  geom_point(data =cv_all, aes(x = CV, y = y, color = method), size = 5) +
  # plot error bars twice
  geom_errorbarh(data =cv_all,
                 aes(xmin = CV - se, xmax = CV + se, y = y, color = method),
                 height = .2, linewidth = 1.6) +
  geom_errorbarh(data =cv_all,
                 aes(xmin = CV - 2*se, xmax = CV + 2*se, y = y, color = method),
                 height = .1, linewidth = 1) +
  # add labels, rather than a legend
  geom_label(data =cv_all, aes(x = CV, y = y, label=method, fill = NULL),
             nudge_y = .2
             ) +
  labs(x = "Root mean squared prediction error (rtmsepe)",
       y = "Density") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")


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
  geom_point(data =cv_all, aes(x = method, y = CV, color = method), size = 5) +
  labs(y = "Root mean squared prediction error (rtmsepe)",
       x = "Estimation method") +
  theme_bw(base_size = 14)


