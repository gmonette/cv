data("Duncan", package="carData")
m.lm <- lm(prestige ~ income + education, data=Duncan)

cv(m.lm, criterion=mse, seed=123)
cv(m.lm, criterion=rmse, seed=123)
