#' @name Pigs
#' @docType data
#'
#' @title Body Weights of 48 Pigs in 9 Successive Weeks
#'
#' @description This data set appears in Table 3.1 of Diggle, Liang,
#' and Zeger (1994).
#'
#' @usage data("Pigs", package = "cv")
#'
#' @format A data frame with 432 rows and 3 columns.
#' \describe{
#'   \item{id}{Pig id number, 1--48.}
#'   \item{week}{Week number, 1--9.}
#'   \item{weight}{Weight in kg.}
#'   }
#'
#' @source P. J. Diggle, K.-Y. Liang, and S. L. Zeger,
#' \emph{Analysis of Longitudinal Data} (Oxford, 1994).
#'
#' @examples
#' library("lme4")
#' m.p <- lmer(weight ~ week + (1 | id) + (1 | week),
#'             data=Pigs, REML=FALSE,
#'             control=lmerControl(optimizer="bobyqa"))
#' summary(m.p)
#' cv(m.p, clusterVariables=c("id", "week"), k=10, seed=8469)
"Pigs"
