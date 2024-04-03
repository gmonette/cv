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

## ----HSB-data-----------------------------------------------------------------
data("MathAchieve", package = "nlme")
dim(MathAchieve)
head(MathAchieve, 3)
tail(MathAchieve, 3)

data("MathAchSchool", package = "nlme")
dim(MathAchSchool)
head(MathAchSchool, 2)
tail(MathAchSchool, 2)

## ----data---------------------------------------------------------------------
# Parameters:
set.seed(9693)
Nb <- 100     # number of groups
Nw <- 5       # number of individuals within groups
Bb <- 0       # between-group regression coefficient on group mean
SDre <-
  2.0   # between-group SD of random level relative to group mean of x
SDwithin <- 0.5  # within group SD
Bw <- 1          # within group effect of x
Ay <- 10         # intercept for response
Ax <- 20         # starting level of x
Nx <- Nw * 10    # number of distinct x values

Data <- data.frame(group = factor(rep(1:Nb, each = Nw)),
                   x = Ax + rep(1:Nx, length.out = Nw * Nb)) |>
  within({
    xm  <- ave(x, group, FUN = mean) # within-group mean
    y <- Ay +
      Bb * xm +                      # contextual effect
      Bw * (x - xm) +                # within-group effect
      rnorm(Nb, sd = SDre)[group] +  # random level by group
      rnorm(Nb * Nw, sd = SDwithin)  # random error within groups
  })

## ----plot1--------------------------------------------------------------------
library("lattice")
library("latticeExtra")
plot <- xyplot(
  y ~ x,
  data = Data[1:Nx,],
  group = group,
  ylim = c(4, 16),
  par.settings = list(superpose.symbol = list(pch = 1, cex =
                                                0.7))
) +
  layer(panel.ellipse(..., center.cex = 0))
plot # display graph

## -----------------------------------------------------------------------------
summary(lm(y ~ x, data=Data))

## ----include=FALSE, echo=FALSE------------------------------------------------
library(lme4) # necessary for some reason to knit vignette in RStudio, harmless otherwise

## -----------------------------------------------------------------------------
# random intercept only:
mod.0 <- lmer(y ~ 1 + (1 | group), Data)
summary(mod.0)

## -----------------------------------------------------------------------------
# effect of x and random intercept:
mod.1 <- lmer(y ~ x + (1 | group), Data)

# effect of x, contextual (student) mean of x, and random intercept:
mod.2 <- lmer(y ~ x + xm + (1 | group), Data)
        # equivalent to y ~ I(x - xm) + xm + (1 | group)

# model generating the data (where Bb = 0)
mod.3 <- lmer(y ~ I(x - xm) + (1 | group), Data)

## -----------------------------------------------------------------------------
Data <- within(Data, {
  fit_mod0.fe <- predict(mod.0, re.form = ~ 0) # fixed effects only
  fit_mod0.re <- predict(mod.0) # fixed and random effects (BLUPs)
  fit_mod1.fe <- predict(mod.1, re.form = ~ 0)
  fit_mod1.re <- predict(mod.1)
  fit_mod2.fe <- predict(mod.2, re.form = ~ 0)
  fit_mod2.re <- predict(mod.2)
  fit_mod3.fe <- predict(mod.3, re.form = ~ 0)
  fit_mod3.re <- predict(mod.3)
})

## -----------------------------------------------------------------------------
Data_long <- reshape(Data[1:Nx, ], direction = "long", sep = ".", 
              timevar = "effect", varying = grep("\\.", names(Data[1:Nx, ])))
Data_long$id <- 1:nrow(Data_long)
Data_long <- reshape(Data_long, direction = "long", sep = "_", 
              timevar = "modelcode",  varying = grep("_", names(Data_long)))
Data_long$model <- factor(
  c("~ 1", "~ 1 + x", "~ 1 + x + xm", "~ 1 + I(x - xm)")
  [match(Data_long$modelcode, c("mod0", "mod1", "mod2", "mod3"))]
)

## ----plot-fits-mod0-----------------------------------------------------------
(
  plot +
    xyplot(
      fit ~ x,
      subset(Data_long, modelcode == "mod0" & effect == "fe"),
      groups = group,
      type = "l",
      lwd = 2
    ) +
    xyplot(
      fit ~ x,
      subset(Data_long, modelcode == "mod0" &  effect == "re"),
      groups = group,
      type = "l",
      lwd = 2,
      lty = 3
    )
) |> update(
  main="Model: y ~ 1 + (1 | group)",
  key=list(
    corner=c(0.05, 0.05),
    text=list(c("fixed effects only","fixed and random")),
    lines=list(lty=c(1, 3))))

## -----------------------------------------------------------------------------
summary(mod.1)

## ----plot-fits-mod1-----------------------------------------------------------
(
  plot +
    xyplot(
      fit ~ x,
      subset(Data_long, modelcode == "mod1" & effect == "fe"),
      groups = group,
      type = "l",
      lwd = 2
    ) +
    xyplot(
      fit ~ x,
      subset(Data_long, modelcode == "mod1" & effect == "re"),
      groups = group,
      type = "l",
      lwd = 2,
      lty = 3
    )
) |> update(
  main="Model: y ~ 1 + x + (1 | group)",
  ylim=c(-15, 35),
  key=list(
    corner=c(0.95, 0.05),
    text=list(c("fixed effects only","fixed and random")),
    lines=list(lty=c(1, 3))))

## -----------------------------------------------------------------------------
summary(mod.2)

## ----plot-fits-mod2-----------------------------------------------------------
(
  plot +
    xyplot(
      fit ~ x,
      subset(Data_long, modelcode == "mod2" & effect == "fe"),
      groups = group,
      type = "l",
      lwd = 2
    ) +
    xyplot(
      fit ~ x,
      subset(Data_long, modelcode == "mod2" & effect == "re"),
      groups = group,
      type = "l",
      lwd = 2,
      lty = 3
    )
) |> update(
  main="Model: y ~ 1 + x + xm + (1 | group)",
  ylim=c(4, 16),
  key=list(
    corner=c(0.05, 0.05),
    text=list(c("fixed effects only","fixed and random")),
    lines=list(lty=c(1, 3))))

## -----------------------------------------------------------------------------
summary(mod.3)

## ----plot-fits-mod3-----------------------------------------------------------
(
  plot +
    xyplot(
      fit ~ x,
      subset(Data_long, modelcode == "mod3" & effect == "fe"),
      groups = group,
      type = "l",
      lwd = 2
    ) +
    xyplot(
      fit ~ x,
      subset(Data_long, modelcode == "mod3" & effect == "re"),
      groups = group,
      type = "l",
      lwd = 2,
      lty = 3
    )
) |> update(
  main="Model: y ~ 1 + I(x - xm) + (1 | group)",
  ylim=c(4, 16),
  key=list(
    corner=c(0.05, 0.05),
    text=list(c("fixed effects only","fixed and random")),
    lines=list(lty=c(1, 3))))

## ----echo=FALSE,include=FALSE-------------------------------------------------
library(cv) # unclear why it's necessary to reload cv

## ----cross-validation-clusters------------------------------------------------
modlist <- models(
  "~ 1" = mod.0,
  "~ 1 + x" = mod.1,
  "~ 1 + x + xm" = mod.2,
  "~ 1 + I(x - xm)" = mod.3
)
cvs_clusters <-
  cv(
    modlist,
    data = Data,
    cluster = "group",
    k = 10,
    seed = 6449
  )
plot(cvs_clusters, main = "Model Comparison, Cluster-Based CV")

## ----cross-validation-cases---------------------------------------------------
cvs_cases <- cv(modlist, data = Data, seed = 9693)
plot(cvs_cases, main = "Model Comparison, Case-Based CV")

## ----pigs---------------------------------------------------------------------
head(Pigs, 9)
head(xtabs( ~ id + week, data = Pigs), 3)
tail(xtabs( ~ id + week, data = Pigs), 3)

## ----pigs-graph---------------------------------------------------------------
plot(weight ~ week, data = Pigs, type = "n")
for (i in unique(Pigs$id)) {
  with(Pigs, lines(
    x = 1:9,
    y = Pigs[id == i, "weight"],
    col = "gray"
  ))
}
abline(lm(weight ~ week, data = Pigs),
       col = "blue",
       lwd = 2)
lines(
  with(Pigs, loess.smooth(week, weight, span = 0.5)),
  col = "magenta",
  lty = 2,
  lwd = 2
)

## ----pigs-lmer----------------------------------------------------------------
m.p <- lmer(
  weight ~ week + (1 | id) + (1 | week),
  data = Pigs,
  REML = FALSE, # i.e., ML
  control = lmerControl(optimizer = "bobyqa")
)
summary(m.p)

## ----pigs-cv------------------------------------------------------------------
cv(m.p, clusterVariables = "id")

cv(m.p, clusterVariables = "week")

cv(
  m.p,
  clusterVariables = c("id", "week"),
  k = 10,
  seed = 8469
)

## ----pigs-cv-cases------------------------------------------------------------
cv(m.p, k = 10, seed = 8469)

## ----coda, include = FALSE----------------------------------------------------
options(.opts)

