#' SCR: Sequential Covariate-Adjusted Randomization
#'
#' Implements the sequential covariate-adjusted randomization (SCR) of
#' Yang, Qin, Li and Hu (2024), which assigns participants one at a time by
#' hierarchically minimizing the marginal imbalance and a modified
#' Mahalanobis distance.
#'
#' Main functions:
#' \itemize{
#'   \item [SCR()]: allocate a sequence of units.
#'   \item [scr_next()]: allocate one newly arrived unit.
#'   \item [scr_imbalance()]: covariate and marginal imbalance measures.
#'   \item [scr_test()]: randomization test for the treatment effect.
#'   \item [smart()]: multi-arm extension.
#' }
#'
#' @references Yang, H., Qin, Y., Li, Y., and Hu, F. (2024). Sequential
#'   covariate-adjusted randomization via hierarchically minimizing
#'   Mahalanobis distance and marginal imbalance. *Biometrics*, 80(2),
#'   ujae047. \doi{10.1093/biomtc/ujae047}
#' @keywords internal
"_PACKAGE"
