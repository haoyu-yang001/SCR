#' Randomization test for the treatment effect under SCR
#'
#' Tests the null hypothesis of no treatment effect after a trial randomized
#' by [SCR()]. The reference distribution of the t statistic is obtained by
#' re-running SCR `B` times on the observed covariates (in the same order of
#' arrival) with the outcomes held fixed, as in the supplementary code of
#' Yang, Qin, Li and Hu (2024).
#'
#' @details
#' The test statistic is the t statistic of the treatment indicator in the
#' linear model \eqn{y_i = \mu_2 + \tau T_i + \beta^\top x_i + \epsilon_i}
#' (`adjust = TRUE`; categorical covariates enter as dummy variables) or
#' \eqn{y_i = \mu_2 + \tau T_i + \epsilon_i} (`adjust = FALSE`, the
#' difference-in-means estimator). The p-value is
#' \eqn{(1 + \#\{b: |t_b| \ge |t_{obs}|\}) / (B + 1)}, which counts the
#' observed allocation as one of the re-randomizations, so it is never
#' smaller than \eqn{1/(B + 1)}.
#'
#' Both treatment effect estimates are reported: the difference-in-means
#' estimator \eqn{\hat\tau} (Equation 2 of the paper) and the
#' regression-adjusted estimator \eqn{\tilde\tau} (Equation 4).
#'
#' @param y Numeric outcome vector.
#' @param assignment Observed 0/1 assignments.
#' @param B Number of re-randomizations.
#' @param adjust Logical; adjust for covariates in the test statistic.
#' @inheritParams SCR
#'
#' @return An object of class `"htest"`, with the extra component `t_null`
#'   holding the `B` re-randomized statistics.
#'
#' @examples
#' set.seed(3)
#' n <- 100; p <- 4
#' x <- matrix(rnorm(n * p), n, p)
#' tr <- SCR(x, d = 5)$assignment
#' y <- 0.5 * tr + drop(x %*% rep(1, p)) + rnorm(n)
#' scr_test(y, x, tr, B = 200)
#' @export
scr_test <- function(y, covariate, assignment, B = 1000, adjust = TRUE,
                     d = 5, q1 = 0.75, q2 = 0.85, categorical = NULL, w = 0.7) {
  dname <- paste(deparse1(substitute(y)), "by", deparse1(substitute(assignment)))
  X <- .as_covariate_matrix(covariate)
  Zc <- .as_categorical_df(categorical)
  if (!is.numeric(y) || anyNA(y)) stop("`y` must be a numeric vector without missing values.", call. = FALSE)
  n <- length(y)
  if (is.null(X) && is.null(Zc)) stop("Supply `covariate`, `categorical`, or both.", call. = FALSE)
  if ((!is.null(X) && nrow(X) != n) || (!is.null(Zc) && nrow(Zc) != n) ||
      length(assignment) != n) {
    stop("`y`, `covariate`, `categorical` and `assignment` must describe the same units.", call. = FALSE)
  }
  tr <- .check_assignment(assignment, n)
  if (!is.numeric(B) || length(B) != 1L || B < 1) stop("`B` must be a positive integer.", call. = FALSE)
  B <- as.integer(B)
  .check_d(d); .check_prob(q1, "q1"); .check_prob(q2, "q2"); .check_w(w)

  # Design matrix of covariates used for adjustment.
  Zadj <- NULL
  if (!is.null(X)) Zadj <- X
  if (!is.null(Zc)) {
    dummies <- stats::model.matrix(~ ., data = as.data.frame(lapply(Zc, factor)))[, -1L, drop = FALSE]
    Zadj <- cbind(Zadj, dummies)
  }

  t_stat <- function(t) .t_stat(y, t, if (adjust) Zadj else NULL)
  t_obs <- t_stat(tr)
  codes <- .strata_codes(Zc)
  t_null <- vapply(seq_len(B), function(b) {
    t_stat(.scr_run(X, codes, n, integer(0), d, q1, q2, w)$assignment)
  }, numeric(1))

  est <- c(
    "difference in means" = mean(y[tr == 1L]) - mean(y[tr == 0L]),
    "regression-adjusted" = .coef_tr(y, tr, Zadj)
  )
  structure(list(
    statistic = c(t = t_obs),
    parameter = c(B = B),
    p.value = (1 + sum(abs(t_null) >= abs(t_obs))) / (B + 1),
    estimate = est,
    null.value = c("treatment effect" = 0),
    alternative = "two.sided",
    method = sprintf("Randomization test under SCR (%s t statistic)",
                     if (adjust) "covariate-adjusted" else "unadjusted"),
    data.name = dname,
    t_null = t_null
  ), class = "htest")
}

.t_stat <- function(y, tr, Z) {
  fit <- if (is.null(Z)) stats::lm(y ~ tr) else stats::lm(y ~ tr + Z)
  cf <- summary(fit)$coefficients
  if (!"tr" %in% rownames(cf)) return(NA_real_)
  cf["tr", "t value"]
}

.coef_tr <- function(y, tr, Z) {
  fit <- if (is.null(Z)) stats::lm(y ~ tr) else stats::lm(y ~ tr + Z)
  unname(stats::coef(fit)["tr"])
}
