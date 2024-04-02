#' ---
#' title: testing cv plotting ideas
#' ---

library(cv)
library(car)
library(ggplot2)
library(dplyr)

#' ## Simplest example
#' Single model, no replications, folds only

data("Auto", package="ISLR2")
m.auto <- lm(mpg ~ horsepower, data=Auto)
(cv.auto <- cv(m.auto, seed=1234, confint=TRUE))

cv.auto.df <- as.data.frame(cv.auto, columns="criteria")


#' ## First attempt: Plot showing overall mse and the folds on different lines.
#' should show also the adj.mse with confi interval ...
cv.auto.df |>
  mutate(Fold = ifelse(fold==0, "Overall", "Folds")) |>

  ggplot(aes(x=mse, y=Fold, color=Fold, shape=Fold)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = adj.mse.lower, xmax = adj.mse.upper), height=.1)

# but seems I should also show the adjusted.mse and full.mse on the Overall line ????

#' ## Try again, separating folds from overall
cv.auto.folds <- as.data.frame(cv.auto, rows="folds", columns="criteria") |>
  select(fold, mse)
cv.auto.all <- as.data.frame(cv.auto, rows="cv", columns="criteria")

# using car::densityPlot
densityPlot(~ mse, data=cv.auto.folds,
            xlab = "Cross-validation MSE")
y <- 0.01
with(cv.auto.all, {
  arrows(x0 = adj.mse.lower, x1 = adj.mse.upper, y0=y, y1=y,
         lwd = 2, angle=90, code=3, length = 0.1)
  points(x = adjusted.mse, y = y, pch = 16, cex=2)
})


ggplot(data=cv.auto.folds, aes(x=mse)) +
  geom_density(alpha = 0.3, fill="pink") +
  geom_jitter(aes(y = 0), size = 2, width=0, height = 0.002) +
  geom_point(data=cv.auto.all,
             aes(x = adjusted.mse, y = 0.01), size=4) +
  geom_errorbarh(data=cv.auto.all,
                 aes(xmin = adjusted.mse - SE.adj.mse,
                     xmax = adjusted.mse + SE.adj.mse, y = 0.01),
                 height = .005, linewidth = 1.4) +
  xlim(5, 45) +
  labs(x = "Cross-validation MSE",
       y = "Density") +
  theme_bw(base_size = 14)

#' ## Plotting coefficients

cv.auto.coef <- as.data.frame(cv.auto, columns="coefficients") |>
  rename(Intercept = coef.Intercept,
         horsepower = coef.horsepower)

dataEllipse(horsepower ~ Intercept, data=cv.auto.coef,
            levels = 0.68, pch = 15)
# show the overall estimate
points(cv.auto.coef[1, "Intercept"], cv.auto.coef[1, "horsepower"],
       pch = "+", cex = 2)
# add the confidence ellipse for the model
confidenceEllipse(m.auto, levels=0.68, fill=TRUE, col="lightgreen", add=TRUE)


#' ## Replications
#'

cv.auto.reps <- cv(m.auto, seed=1234, confint=TRUE, reps = 20)

cv.auto.reps.df <- as.data.frame(cv.auto.reps)
names(cv.auto.reps.df)
head(cv.auto.reps.df)

cv.auto.reps.folds <- as.data.frame(cv.auto.reps, rows="folds", columns="criteria") |>
  select(rep, fold, mse)

cv.auto.reps.all <- as.data.frame(cv.auto.reps, rows="cv", columns="criteria")




