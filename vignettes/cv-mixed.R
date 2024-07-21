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
library(cv)
library(lme4)
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

## ----include=FALSE, echo=FALSE------------------------------------------------
library("glmmTMB") # necessary for some reason to knit vignette in RStudio, harmless otherwise

## ----parameters---------------------------------------------------------------

# Parameters:

Nb <- 20     # number of patients
Nw <- 5      # number of occasions for each patient

Bb <- 1.0    # between-patient regression coefficient on patient means
Bw <- -0.5   # within-patient effect of x

SD_between <- c(0, 5, 6, 8)               # SD between patients
SD_within <- rep(2.5, length(SD_between)) # SD within patients

Nv <- length(SD_within)       # number of variance profiles
SD_ratio <- paste0('SD ratio = ', SD_between,' / ',SD_within)
SD_ratio <- factor(SD_ratio, levels = SD_ratio)

set.seed(833885) 

Data_template <- expand.grid(patient = 1:Nb, obs = 1:Nw) |>
  within({
    xw <- seq(-2, 2, length.out = Nw)[obs]
    x <- patient + xw
    xm  <- ave(x, patient)   # within-patient mean

    # Scaled random error within each SD_ratio_i group

    re_std <- scale(resid(lm(rnorm(Nb*Nw) ~ x)))
    re_between <- ave(re_std, patient)
    re_within <- re_std - re_between
    re_between <- scale(re_between)/sqrt(Nw)
    re_within <- scale(re_within)
  })

Data <- do.call(
  rbind,
  lapply(
    1:Nv,
    function(i) {
      cbind(Data_template, SD_ratio_i = i)
    }
  )
)

Data <- within(
  Data,
  {
    SD_within_ <- SD_within[SD_ratio_i]
    SD_between_ <- SD_between[SD_ratio_i]
    SD_ratio <- SD_ratio[SD_ratio_i]
    y <- 10 +
      Bb * xm +                  # contextual effect
      Bw * (x - xm) +            # within-patient effect
      SD_within_ * re_within +   # within patient random effect
      SD_between_ * re_between   # adjustment to between patient random effect
  }
)

## ----plot1--------------------------------------------------------------------
library("lattice")
library("latticeExtra")
plot <- xyplot(y ~ x | SD_ratio, data = Data, group = patient,
           layout = c(Nv, 1),
           par.settings = list(superpose.symbol = list(pch = 1, cex = 0.7))) +
      layer(panel.ellipse(..., center.pch = 16, center.cex = 1.5,  
                          level = 0.5),
            panel.abline(a = 10, b = 1))
plot # display graph

## ----model-fits---------------------------------------------------------------
model.formulas <- c(
  ' ~ 1'             =  y ~ 1 + (1 | patient),
  '~ 1 + x'          =  y ~ 1 + x + (1 | patient),
  '~ 1 + x + xm'     =  y ~ 1 + x + xm + (1 | patient),
  '~ 1 + I(x - xm)'  =  y ~ 1 + I(x - xm) + (1 | patient)
)

fits <- lapply(split(Data, ~ SD_ratio),
               function(d) {
                 lapply(model.formulas, function(form) {
                   glmmTMB(form, data = d)
                 })
               })

## ----predict------------------------------------------------------------------
# predicted fixed and random effects:
pred.BLUPs <- lapply(fits, lapply, predict)
# predicted fixed effects:
pred.fixed <- lapply(fits, lapply, predict, re.form = ~0)  

## ----data-predictions---------------------------------------------------------
Dataf <- lapply(split(Data, ~ SD_ratio),
                    function(d) {
                      lapply(names(model.formulas), 
                             function(form) cbind(d, formula = form))
                    }) |> 
             lapply(function(dlist) do.call(rbind, dlist)) |> 
             do.call(rbind, args = _)
    
Dataf <- within(
  Dataf,
  {
    pred.fixed <- unlist(pred.fixed)
    pred.BLUPs <- unlist(pred.BLUPs)
    panel <- factor(formula, levels = c(names(model.formulas), 'data'))
  }
)

Data$panel <- factor('data', levels = c(names(model.formulas), 'data'))

## ----plot-fits-fixed----------------------------------------------------------
{
  xyplot(y ~ x |SD_ratio * panel, Data,
         groups = patient, type = 'n', 
         par.strip.text = list(cex = 0.7),
         drop.unused.levels = FALSE) +
    glayer(panel.ellipse(..., center.pch = 16, center.cex = 0.5,  
                         level = 0.5),
           panel.abline(a = 10, b = 1)) +
    xyplot(pred.fixed  ~ x |SD_ratio * panel, Dataf, type = 'l', 
           groups = patient,
           drop.unused.levels = F,
           ylab = 'fixed-effect predictions') 
}|> 
  useOuterStrips() |> print()

## ----plot-fits-blups----------------------------------------------------------
{
  xyplot(y ~ x | SD_ratio * panel, Data,
         groups = patient, type = 'n',
         drop.unused.levels = F, 
         par.strip.text = list(cex = 0.7),
         ylab = 'fixed- and random-effect predictions (BLUPS)') +
    glayer(panel.ellipse(..., center.pch = 16, center.cex = 0.5,  
                         level = 0.5),
           panel.abline(a = 10, b = 1)) +
    xyplot(pred.BLUPs  ~ x | SD_ratio * panel, Dataf, type = 'l', 
           groups = patient,
           drop.unused.levels = F) 
}|> 
  useOuterStrips() |> print()

## ----echo=FALSE,include=FALSE-------------------------------------------------
library(cv) # unclear why it's necessary to reload cv
.opts <- options(warn = -2) 

## ----cross-validation,cache=TRUE----------------------------------------------

model_lists <- lapply(fits, function(fitlist) do.call(models, fitlist))

cvs_cases <-
  lapply(1:Nv,
         function(i){
           cv(model_lists[[i]], k = 10, 
              data = split(Data, ~ SD_ratio)[[i]])
         })

cvs_clusters <-
  lapply(1:Nv,
         function(i){
           cv(model_lists[[i]],  k = 10, 
              data = split(Data, ~SD_ratio)[[i]], 
              clusterVariables = 'patient')
         })

## ----plot-cv-example----------------------------------------------------------
plot(cvs_clusters[[2]], main="Comparison of Fixed Effects")

## ----cross-validation-data----------------------------------------------------
names(cvs_clusters) <- names(cvs_cases) <- SD_ratio 
dsummary <- expand.grid(SD_ratio_i = names(cvs_cases), model = names(cvs_cases[[1]]))

dsummary$cases <-
  sapply(1:nrow(dsummary), function(i){
    with(dsummary[i,], cvs_cases[[SD_ratio_i]][[model]][['CV crit']])
  })

dsummary$clusters <-
  sapply(1:nrow(dsummary), function(i){
    with(dsummary[i,], cvs_clusters[[SD_ratio_i]][[model]][['CV crit']])
  })

## ----cross-validation-data-plot-----------------------------------------------
xyplot(cases + clusters ~ model|SD_ratio_i, dsummary,
       auto.key = list(space = 'top', reverse.rows = T, columns = 2), type = 'b',
       xlab = "Fixed Effects",
       ylab = 'CV criterion (MSE)',
       layout= c(Nv,1),
       par.settings =
         list(superpose.line=list(lty = c(2, 3), lwd = 3),
              superpose.symbol=list(pch = 15:16, cex = 1.5)),
       scales = list(y = list(log = TRUE), x = list(alternating = F, rot = 60))) |> print()

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

