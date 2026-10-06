#' Covariate and marginal imbalance between two groups
#'
#' Computes the imbalance measures used in Yang, Qin, Li and Hu (2024) for a
#' given 0/1 allocation.
#'
#' @details
#' With \eqn{\mathrm{cov}(x)} the sample covariance matrix of all units
#' (a Moore-Penrose inverse is used when it is singular):
#' \itemize{
#'   \item `Mahalanobis_Distance`: \eqn{M(n) = n(\bar x_1 - \bar x_0)^\top
#'     \mathrm{cov}(x)^{-1}(\bar x_1 - \bar x_0)}.
#'   \item `modified_Mahalanobis_Distance`: \eqn{\widetilde M(n)}, the same
#'     quadratic form with group sums in place of group means (Equation 1).
#'   \item `within_stratum`: \eqn{\|D_n\|_2^2}, where \eqn{D_n} holds the
#'     treatment-minus-control counts in each stratum formed by the
#'     categorical covariates.
#'   \item `W`: the combined measure \eqn{w \widetilde M(n) + (1 - w)
#'     \|D_n\|_2^2} (or whichever part is available).
#'   \item `marginal`: \eqn{n_1 - n_0}.
#' }
#'
#' @param covariate Numeric matrix or data frame of continuous covariates, or
#'   `NULL`.
#' @param assignment 0/1 vector of assignments.
#' @param categorical Optional data frame (or vector) of categorical
#'   covariates.
#' @param w Weight of the modified Mahalanobis distance in `W`.
#'
#' @return A named numeric vector. Measures that do not apply are `NA`.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200 * 4), 200, 4)
#' scr_imbalance(x, SCR(x)$assignment)
#' scr_imbalance(x, rbinom(200, 1, 0.5))  # complete randomization
#' @export
scr_imbalance <- function(covariate, assignment, categorical = NULL, w = 0.7) {
  X <- .as_covariate_matrix(covariate)
  Z <- .as_categorical_df(categorical)
  n <- length(assignment)
  if (!is.null(X) && nrow(X) != n) {
    stop("`covariate` and `assignment` have different numbers of units.", call. = FALSE)
  }
  if (!is.null(Z) && nrow(Z) != n) {
    stop("`categorical` and `assignment` have different numbers of units.", call. = FALSE)
  }
  tr <- .check_assignment(assignment, n)
  if (length(tr) != n) stop("`assignment` must not be NA.", call. = FALSE)
  .check_w(w)

  M <- Mt <- Ds <- NA_real_
  if (!is.null(X) && ncol(X) > 0L) {
    if (n >= 2L) {
      A <- MASS::ginv(stats::cov(X))
      S <- colSums(X[tr == 1L, , drop = FALSE]) - colSums(X[tr == 0L, , drop = FALSE])
      Mt <- drop(crossprod(S, A %*% S))
      if (any(tr == 1L) && any(tr == 0L)) {
        u <- colMeans(X[tr == 1L, , drop = FALSE]) - colMeans(X[tr == 0L, , drop = FALSE])
        M <- n * drop(crossprod(u, A %*% u))
      }
    }
  }
  if (!is.null(Z)) {
    codes <- .strata_codes(Z)
    Dn <- tabulate(codes[tr == 1L], max(codes)) - tabulate(codes[tr == 0L], max(codes))
    Ds <- sum(Dn^2)
  }
  W <- if (!is.na(Mt) && !is.na(Ds)) w * Mt + (1 - w) * Ds else if (!is.na(Mt)) Mt else Ds

  c(Mahalanobis_Distance = M,
    modified_Mahalanobis_Distance = Mt,
    within_stratum = Ds,
    W = W,
    marginal = sum(tr == 1L) - sum(tr == 0L))
}
