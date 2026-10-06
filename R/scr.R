#' Sequential Covariate-Adjusted Randomization (SCR)
#'
#' Allocates units one at a time to two treatment groups by hierarchically
#' minimizing the marginal imbalance (difference in group sizes) and the
#' covariate imbalance (a modified Mahalanobis distance), following
#' Algorithm 1 of Yang, Qin, Li and Hu (2024).
#'
#' @details
#' Let \eqn{T_i \in \{0, 1\}} be the assignment of unit \eqn{i}
#' (1 = treatment, 0 = control) and \eqn{x_i} its continuous covariates.
#' The covariate imbalance after \eqn{i} units is the modified Mahalanobis
#' distance
#' \deqn{\widetilde{M}(i) = \Big(\sum_{j:T_j=1} x_j - \sum_{j:T_j=0} x_j\Big)^\top
#'   \mathrm{cov}(x)^{-1} \Big(\sum_{j:T_j=1} x_j - \sum_{j:T_j=0} x_j\Big),}
#' where \eqn{\mathrm{cov}(x)} is the sample covariance of the covariates of
#' the first \eqn{i} units (a Moore-Penrose inverse is used when it is
#' singular).
#'
#' The procedure is:
#' \enumerate{
#'   \item The first two units are assigned to the two treatments separately,
#'     in random order.
#'   \item For unit \eqn{i \ge 3}, let \eqn{n_1} and \eqn{n_0} be the current
#'     group sizes. If \eqn{|n_1 - n_0| < d}, compute the "potential"
#'     imbalance \eqn{\widetilde{M}_1(i)} (if \eqn{T_i = 1}) and
#'     \eqn{\widetilde{M}_0(i)} (if \eqn{T_i = 0}), and set \eqn{T_i = 1}
#'     with probability \eqn{q_1}, \eqn{1 - q_1} or \eqn{0.5} when
#'     \eqn{\widetilde{M}_1(i)} is smaller than, larger than or equal to
#'     \eqn{\widetilde{M}_0(i)}.
#'   \item Otherwise, set \eqn{T_i = 1} with probability \eqn{q_2} if
#'     \eqn{n_1 < n_0} and \eqn{1 - q_2} if \eqn{n_1 > n_0}.
#' }
#' The marginal constraint uses the strict inequality \eqn{|n_1 - n_0| < d},
#' as in the code accompanying the paper (SCR 1.0.0), so that \eqn{d = 1}
#' corresponds to an ARM-like procedure and \eqn{d = n} balances covariates
#' only. With \eqn{d = 0} the procedure reduces to Efron's biased coin
#' design with bias \eqn{q_2}.
#'
#' **Continuous and categorical covariates (Section 2.2).** If
#' `categorical` is supplied, units are grouped into strata defined by all
#' combinations of the categorical covariates and the imbalance measure
#' becomes
#' \deqn{W(i) = w \widetilde{M}(i) + (1 - w) \|D_i\|_2^2,}
#' where \eqn{D_i} is the vector of within-stratum differences in group
#' sizes. If `covariate` is `NULL`, only the within-stratum imbalance is
#' used.
#'
#' **Sequential use.** `SCR()` allocates all units of `covariate` in row
#' order. In a trial where units arrive one by one, use [scr_next()] to
#' assign each new arrival; with the same random seed, calling
#' [scr_next()] repeatedly gives exactly the same assignments as `SCR()`.
#'
#' @param covariate A numeric matrix or data frame of continuous covariates,
#'   one row per unit, in order of arrival. May be `NULL` if `categorical`
#'   is given.
#' @param assignment Optional 0/1 vector with the assignments of the first
#'   `length(assignment)` units, if they have already been allocated.
#'   `NULL` or `NA` (default) if no unit has been allocated.
#' @param d Marginal imbalance constraint (non-negative). Covariate balance
#'   is targeted only while \eqn{|n_1 - n_0| < d}. The paper recommends a
#'   moderate value such as 5 or 10.
#' @param q1 Biased coin probability of the covariate-adaptive step, in
#'   \eqn{[0.5, 1]}. Default 0.75.
#' @param q2 Biased coin probability of the marginal (sample size) step, in
#'   \eqn{[0.5, 1]}. Default 0.85.
#' @param categorical Optional data frame (or vector) of categorical
#'   covariates, one row per unit, aligned with `covariate`.
#' @param w Weight of the modified Mahalanobis distance in \eqn{W(i)} when
#'   both continuous and categorical covariates are given. Default 0.7, as
#'   in the paper.
#'
#' @return An object of class `"SCR"`, a list with components
#' \item{assignment}{Integer vector of assignments (1 = treatment,
#'   0 = control).}
#' \item{prob}{Probability \eqn{P(T_i = 1)} used for each unit (`NA` for
#'   units whose assignment was supplied).}
#' \item{rule}{The step used for each unit: `"given"`, `"burn-in"`,
#'   `"covariate"` or `"marginal"`.}
#' \item{sample_size}{Group sizes.}
#' \item{Mahalanobis_Distance}{Mahalanobis distance
#'   \eqn{M(n) = n(\bar x_1 - \bar x_0)^\top \mathrm{cov}(x)^{-1}
#'   (\bar x_1 - \bar x_0)} between the two groups.}
#' \item{imbalance}{All imbalance measures, see [scr_imbalance()].}
#' \item{parameters}{The design parameters.}
#'
#' @references Yang, H., Qin, Y., Li, Y., and Hu, F. (2024). Sequential
#'   covariate-adjusted randomization via hierarchically minimizing
#'   Mahalanobis distance and marginal imbalance. *Biometrics*, 80(2),
#'   ujae047. \doi{10.1093/biomtc/ujae047}
#'
#' @seealso [scr_next()] for allocating one new unit, [scr_imbalance()],
#'   [scr_test()], [smart()].
#'
#' @examples
#' set.seed(1)
#' n <- 100; p <- 5
#' x <- matrix(rnorm(n * p), n, p)
#' fit <- SCR(x, d = 5, q1 = 0.75, q2 = 0.85)
#' fit
#' table(fit$assignment)
#'
#' # Continue a trial in which the first three units were already assigned
#' fit2 <- SCR(x, assignment = c(1, 0, 0))
#' head(fit2$assignment)
#'
#' # Continuous and categorical covariates
#' z <- data.frame(sex = sample(c("F", "M"), n, TRUE),
#'                 site = sample(1:3, n, TRUE))
#' SCR(x, categorical = z, w = 0.7)
#' @export
SCR <- function(covariate, assignment = NULL, d = 5, q1 = 0.75, q2 = 0.85,
                categorical = NULL, w = 0.7) {
  X <- .as_covariate_matrix(covariate)
  Z <- .as_categorical_df(categorical)
  if (is.null(X) && is.null(Z)) {
    stop("Supply `covariate`, `categorical`, or both.", call. = FALSE)
  }
  n <- if (!is.null(X)) nrow(X) else nrow(Z)
  if (!is.null(X) && !is.null(Z) && nrow(Z) != n) {
    stop("`covariate` and `categorical` must have the same number of rows.", call. = FALSE)
  }
  given <- .check_assignment(assignment, n)
  .check_d(d); .check_prob(q1, "q1"); .check_prob(q2, "q2"); .check_w(w)

  res <- .scr_run(X, .strata_codes(Z), n, given, d, q1, q2, w)

  out <- list(
    assignment = res$assignment,
    prob = res$prob,
    rule = res$rule,
    sample_size = c(treatment = sum(res$assignment == 1L),
                    control = sum(res$assignment == 0L)),
    Mahalanobis_Distance = NA_real_,
    imbalance = scr_imbalance(X, res$assignment, categorical = Z, w = w),
    parameters = list(d = d, q1 = q1, q2 = q2, w = w)
  )
  out$Mahalanobis_Distance <- out$imbalance[["Mahalanobis_Distance"]]
  class(out) <- "SCR"
  out
}

#' Assign a newly arrived unit by SCR
#'
#' Given the covariates and assignments of the units already in the trial,
#' returns the SCR assignment of the next unit. This is the function to call
#' in a real trial, where each participant must be randomized on arrival.
#'
#' @param covariate_new Continuous covariates of the new unit (a vector, or a
#'   one-row matrix or data frame). `NULL` if only categorical covariates are
#'   used.
#' @param covariate Continuous covariates of the units already assigned
#'   (rows in order of arrival). `NULL` or zero rows if the new unit is the
#'   first one.
#' @param assignment 0/1 assignments of the units already assigned.
#' @param categorical_new,categorical Categorical covariates of the new unit
#'   and of the units already assigned (see [SCR()]).
#' @inheritParams SCR
#'
#' @return A list with components `assignment` (0 or 1), `prob`
#'   (\eqn{P(T = 1)} used for the draw) and `rule` (`"burn-in"`,
#'   `"covariate"` or `"marginal"`).
#'
#' @examples
#' set.seed(2)
#' x <- matrix(rnorm(40 * 3), 40, 3)
#' tr <- integer(0)
#' for (i in seq_len(nrow(x))) {
#'   new <- scr_next(x[i, ], covariate = x[seq_len(i - 1), , drop = FALSE],
#'                   assignment = tr)
#'   tr <- c(tr, new$assignment)
#' }
#' table(tr)
#'
#' # Same assignments as allocating everyone at once with the same seed
#' set.seed(2)
#' x <- matrix(rnorm(40 * 3), 40, 3)
#' identical(SCR(x)$assignment, tr)
#' @export
scr_next <- function(covariate_new, covariate = NULL, assignment = NULL,
                     d = 5, q1 = 0.75, q2 = 0.85,
                     categorical_new = NULL, categorical = NULL, w = 0.7) {
  x_new <- .as_covariate_matrix(covariate_new, "covariate_new")
  z_new <- .as_categorical_df(categorical_new, "categorical_new")
  if (is.null(x_new) && is.null(z_new)) {
    stop("Supply `covariate_new`, `categorical_new`, or both.", call. = FALSE)
  }
  if (!is.null(x_new) && nrow(x_new) != 1L) x_new <- matrix(x_new, nrow = 1L)
  if (!is.null(z_new) && nrow(z_new) != 1L) {
    stop("`categorical_new` must describe a single unit.", call. = FALSE)
  }

  X <- .as_covariate_matrix(covariate)
  Z <- .as_categorical_df(categorical)
  if (!is.null(x_new)) {
    if (is.null(X)) X <- matrix(numeric(0), 0L, ncol(x_new))
    if (ncol(X) != ncol(x_new)) {
      stop("`covariate_new` and `covariate` have different numbers of covariates.", call. = FALSE)
    }
    X <- rbind(X, x_new)
  } else if (!is.null(X) && nrow(X) > 0L) {
    stop("`covariate` is given but `covariate_new` is NULL.", call. = FALSE)
  } else {
    X <- NULL
  }
  if (!is.null(z_new)) {
    if (is.null(Z)) {
      Z <- z_new[0L, , drop = FALSE]
    } else {
      if (ncol(Z) != ncol(z_new)) {
        stop("`categorical_new` and `categorical` have different numbers of covariates.", call. = FALSE)
      }
      names(z_new) <- names(Z)
    }
    Z <- rbind(Z, z_new)
  } else if (!is.null(Z) && nrow(Z) > 0L) {
    stop("`categorical` is given but `categorical_new` is NULL.", call. = FALSE)
  } else {
    Z <- NULL
  }

  n <- if (!is.null(X)) nrow(X) else nrow(Z)
  if (!is.null(X) && !is.null(Z) && nrow(Z) != n) {
    stop("`covariate` and `categorical` must have the same number of rows.", call. = FALSE)
  }
  given <- .check_assignment(assignment, n)
  if (length(given) != n - 1L) {
    stop("`assignment` must contain one value for each unit already in the trial.", call. = FALSE)
  }
  .check_d(d); .check_prob(q1, "q1"); .check_prob(q2, "q2"); .check_w(w)

  res <- .scr_run(X, .strata_codes(Z), n, given, d, q1, q2, w)
  list(assignment = res$assignment[n], prob = res$prob[n], rule = res$rule[n])
}

# Core SCR loop. X: n x p matrix or NULL; codes: stratum codes or NULL;
# given: assignments of the first length(given) units.
.scr_run <- function(X, codes, n, given, d, q1, q2, w) {
  has_cont <- !is.null(X) && ncol(X) > 0L
  has_cat <- !is.null(codes)
  p <- if (has_cont) ncol(X) else 0L
  state <- list(k = 0L, mean = numeric(p), M2 = matrix(0, p, p),
                D = numeric(p), Ds = numeric(if (has_cat) max(codes) else 0L),
                n1 = 0L, n0 = 0L)
  wc <- if (has_cat) w else 1      # weight of the Mahalanobis part
  wd <- if (has_cont) 1 - w else 1 # weight of the within-stratum part

  assignment <- integer(n)
  prob <- rep(NA_real_, n)
  rule <- character(n)
  n_given <- length(given)

  for (i in seq_len(n)) {
    x <- if (has_cont) X[i, ] else NULL
    s <- if (has_cat) codes[i] else NULL
    if (has_cont) state <- .welford_add(state, x)

    if (i <= n_given) {
      t <- given[i]
      rule[i] <- "given"
    } else {
      n1 <- state$n1; n0 <- state$n0
      if (n1 == 0L || n0 == 0L) {
        # Line 1 of Algorithm 1: the first two units go to different arms.
        pr <- if (n1 == 0L && n0 == 0L) 0.5 else if (n1 == 0L) 1 else 0
        rule[i] <- "burn-in"
      } else if (abs(n1 - n0) < d) {
        # Covariate-adaptive step. W_1(i) - W_0(i) = 4 * diff, where
        # diff = wc * x' cov^- D + wd * D_s.
        diff <- 0; scale <- 0
        if (has_cont) {
          AD <- drop(MASS::ginv(.welford_cov(state)) %*% state$D)
          diff <- diff + wc * sum(x * AD)
          scale <- scale + wc * sum(abs(x * AD))
        }
        if (has_cat) {
          diff <- diff + wd * state$Ds[s]
          scale <- scale + wd * abs(state$Ds[s])
        }
        pr <- if (.is_tie(diff, scale)) 0.5 else if (diff < 0) q1 else 1 - q1
        rule[i] <- "covariate"
      } else {
        # Marginal step: favour the smaller group.
        pr <- if (n1 < n0) q2 else if (n1 > n0) 1 - q2 else 0.5
        rule[i] <- "marginal"
      }
      t <- .draw(pr)
      prob[i] <- pr
    }

    assignment[i] <- t
    sgn <- 2L * t - 1L
    if (has_cont) state$D <- state$D + sgn * x
    if (has_cat) state$Ds[s] <- state$Ds[s] + sgn
    if (t == 1L) state$n1 <- state$n1 + 1L else state$n0 <- state$n0 + 1L
  }
  list(assignment = assignment, prob = prob, rule = rule)
}

#' @export
print.SCR <- function(x, digits = 4, ...) {
  cat("Sequential Covariate-Adjusted Randomization (SCR)\n")
  pr <- x$parameters
  cat(sprintf("  d = %s, q1 = %s, q2 = %s", pr$d, pr$q1, pr$q2))
  if (!is.na(x$imbalance[["within_stratum"]])) cat(sprintf(", w = %s", pr$w))
  cat("\n\n")
  cat(sprintf("  Units: %d  (treatment = %d, control = %d, difference = %d)\n",
              length(x$assignment), x$sample_size[["treatment"]],
              x$sample_size[["control"]],
              x$sample_size[["treatment"]] - x$sample_size[["control"]]))
  im <- x$imbalance
  if (!is.na(im[["Mahalanobis_Distance"]])) {
    cat(sprintf("  Mahalanobis distance M(n):            %s\n",
                format(im[["Mahalanobis_Distance"]], digits = digits)))
    cat(sprintf("  Modified Mahalanobis distance:        %s\n",
                format(im[["modified_Mahalanobis_Distance"]], digits = digits)))
  }
  if (!is.na(im[["within_stratum"]])) {
    cat(sprintf("  Within-stratum imbalance ||D_n||^2:   %s\n",
                format(im[["within_stratum"]], digits = digits)))
  }
  if (!is.na(im[["W"]])) {
    cat(sprintf("  Combined imbalance W(n):              %s\n",
                format(im[["W"]], digits = digits)))
  }
  cat("\n  Steps used: ")
  tab <- table(factor(x$rule, levels = c("given", "burn-in", "covariate", "marginal")))
  tab <- tab[tab > 0]
  cat(paste(sprintf("%s = %d", names(tab), as.integer(tab)), collapse = ", "), "\n")
  invisible(x)
}
