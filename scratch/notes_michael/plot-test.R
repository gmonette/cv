library(cv)

data("Auto", package="ISLR2")
m.auto <- lm(mpg ~ horsepower, data=Auto)
(cv.auto <- cv(m.auto, seed=1234, confint=TRUE))

cv.auto.df <- as.data.frame(cv.auto)

library(ggplot2)
library(dplyr)

# First attempt: Plot showing overall mse and the folds on different lines.
# should show also the adj.mse with confi interval ...
cv.auto.df |>
  mutate(Fold = ifelse(fold==0, "Overall", "Folds")) |>

  ggplot(aes(x=mse, y=Fold, color=Fold, shape=Fold)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = adj.mse.lower, xmax = adj.mse.upper), height=.1)

# but seems I should also show the adjusted.mse and full.mse on the Overall line ????
