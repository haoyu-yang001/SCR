#' Sequential covariate-adjusted randomization for multi-arm trials
#'
#' A multi-arm version of SCR (as in SCR 1.0.0): units are allocated one at a
#' time to \eqn{K} treatments, balancing the covariates through pairwise
#' modified Mahalanobis distances whenever the group sizes are within the
#' marginal constraint `d`.
#'
#' @details
#' Let \eqn{S_k} be the sum of the covariates of the units in arm \eqn{k}
#' and \eqn{\mathrm{cov}(x)} the sample covariance of the covariates of the
#' units so far. The imbalance between arms \eqn{s} and \eqn{t} is
#' \eqn{(S_s - S_t)^\top \mathrm{cov}(x)^{-1} (S_s - S_t)}, and the overall
#' imbalance is the mean, maximum or median over all pairs (`method`).
#'
#' \enumerate{
#'   \item The first \eqn{K} units are assigned to the \eqn{K} arms in random
#'     order (more generally, while some arm is empty the unit goes to an
#'     empty arm chosen at random).
#'   \item If the largest and smallest group sizes differ by less than `d`,
#'     the arm minimizing the overall potential imbalance is chosen with
#'     probability `q` and each of the other arms with probability
#'     \eqn{(1 - q)/(K - 1)}. Arms tied for the minimum share the
#'     preference equally.
#'   \item Otherwise the unit is assigned at random to one of the arms with
#'     the smallest group size.
#' }
#'
#' @param K Number of treatment arms (at least 2).
#' @param q Biased coin probability, in \eqn{[1/K, 1]}. Default 0.75.
#' @param method How pairwise imbalances are combined: `"mean"`, `"max"` or
#'   `"median"`.
#' @param assignment Optional vector of assignments (coded `1, ..., K`) of
#'   the first units.
#' @inheritParams SCR
#'
#' @return An object of class `"smart"`, a list with components
#' \item{assignment}{Integer vector of arms `1, ..., K`.}
#' \item{prob}{An \eqn{n \times K} matrix of the allocation probabilities
#'   used for each unit (`NA` for units whose assignment was supplied).}
#' \item{rule}{The step used for each unit.}
#' \item{sample_size}{Group sizes.}
#' \item{pairwise_Mahalanobis}{Mahalanobis distance
#'   \eqn{(2n/K^2)(\bar x_s - \bar x_t)^\top \mathrm{cov}(x)^{-1}
#'   (\bar x_s - \bar x_t)} for each pair of arms.}
#' \item{Mahalanobis_Distance}{The pairwise distances combined by
#'   `method`.}
#' \item{parameters}{The design parameters.}
#'
#' @references Yang, H., Qin, Y., Li, Y., and Hu, F. (2024). Sequential
#'   covariate-adjusted randomization via hierarchically minimizing
#'   Mahalanobis distance and marginal imbalance. *Biometrics*, 80(2),
#'   ujae047. \doi{10.1093/biomtc/ujae047}
#'
#'   Yang, H., Qin, Y., Wang, F., Li, Y., and Hu, F. (2023). Balancing
#'   covariates in multi-arm trials via adaptive randomization.
#'   *Computational Statistics & Data Analysis*, 179, 107642.
#'   \doi{10.1016/j.csda.2022.107642}
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(90 * 4), 90, 4)
#' fit <- smart(x, K = 3, d = 5, q = 0.75)
#' fit
#' @export
smart <- function(covariate, assignment = NULL, K, d = 5, q = 0.75,
                  method = c("mean", "max", "median")) {
  method <- match.arg(method)
  X <- .as_covariate_matrix(covariate, allow_null = FALSE)
  n <- nrow(X); p <- ncol(X)
  if (missing(K) || !is.numeric(K) || length(K) != 1L || K < 2 || K != round(K)) {
    stop("`K` must be an integer of at least 2.", call. = FALSE)
  }
  K <- as.integer(K)
  given <- .check_assignment(assignment, n, levels = seq_len(K))
  .check_d(d)
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q < 1 / K || q > 1) {
    stop("`q` must be a single number in [1/K, 1].", call. = FALSE)
  }
  agg <- switch(method, mean = mean, max = max, median = stats::median)
  pairs <- utils::combn(K, 2L)

  state <- list(k = 0L, mean = numeric(p), M2 = matrix(0, p, p))
  S <- matrix(0, K, p)
  counts <- integer(K)
  assignment <- integer(n)
  prob <- matrix(NA_real_, n, K, dimnames = list(NULL, paste0("arm", seq_len(K))))
  rule <- character(n)

  for (i in seq_len(n)) {
    x <- X[i, ]
    state <- .welford_add(state, x)
    if (i <= length(given)) {
      a <- given[i]
      rule[i] <- "given"
    } else {
      if (any(counts == 0L)) {
        pr <- as.numeric(counts == 0L) / sum(counts == 0L)
        rule[i] <- "burn-in"
      } else if (max(counts) - min(counts) < d) {
        A <- MASS::ginv(.welford_cov(state))
        imb <- vapply(seq_len(K), function(k) {
          Sk <- S
          Sk[k, ] <- Sk[k, ] + x
          G <- Sk %*% A %*% t(Sk)
          agg(G[cbind(pairs[1, ], pairs[1, ])] + G[cbind(pairs[2, ], pairs[2, ])] -
                2 * G[t(pairs)])
        }, numeric(1))
        best <- .is_tie(imb - min(imb), abs(min(imb)))
        m <- sum(best)
        other <- (1 - q) / (K - 1)
        pr <- ifelse(best, q / m + (1 - 1 / m) * other, other)
        rule[i] <- "covariate"
      } else {
        small <- counts == min(counts)
        pr <- as.numeric(small) / sum(small)
        rule[i] <- "marginal"
      }
      a <- findInterval(stats::runif(1), cumsum(pr), rightmost.closed = TRUE) + 1L
      a <- min(a, K)
      while (pr[a] == 0) a <- a - 1L  # guard against rounding at the boundary
      prob[i, ] <- pr
    }
    assignment[i] <- a
    S[a, ] <- S[a, ] + x
    counts[a] <- counts[a] + 1L
  }

  pw <- .pairwise_mahalanobis(X, assignment, K)
  out <- list(
    assignment = assignment,
    prob = prob,
    rule = rule,
    sample_size = stats::setNames(counts, paste0("arm", seq_len(K))),
    pairwise_Mahalanobis = pw,
    Mahalanobis_Distance = if (anyNA(pw)) NA_real_ else agg(pw),
    parameters = list(K = K, d = d, q = q, method = method)
  )
  class(out) <- "smart"
  out
}

.pairwise_mahalanobis <- function(X, assignment, K) {
  n <- nrow(X)
  pairs <- utils::combn(K, 2L)
  A <- if (n >= 2L) MASS::ginv(stats::cov(X)) else NULL
  vals <- apply(pairs, 2L, function(st) {
    if (is.null(A) || !any(assignment == st[1]) || !any(assignment == st[2])) return(NA_real_)
    u <- colMeans(X[assignment == st[1], , drop = FALSE]) -
      colMeans(X[assignment == st[2], , drop = FALSE])
    2 / K^2 * n * drop(crossprod(u, A %*% u))
  })
  stats::setNames(vals, apply(pairs, 2L, paste, collapse = "-"))
}

#' @export
print.smart <- function(x, digits = 4, ...) {
  pr <- x$parameters
  cat(sprintf("Multi-arm sequential covariate-adjusted randomization (K = %d)\n", pr$K))
  cat(sprintf("  d = %s, q = %s, method = %s\n\n", pr$d, pr$q, pr$method))
  cat("  Group sizes: ")
  cat(paste(sprintf("%s = %d", names(x$sample_size), x$sample_size), collapse = ", "), "\n")
  cat("  Pairwise Mahalanobis distances:\n")
  print(round(x$pairwise_Mahalanobis, digits))
  cat(sprintf("  Combined (%s): %s\n", pr$method, format(x$Mahalanobis_Distance, digits = digits)))
  invisible(x)
}
