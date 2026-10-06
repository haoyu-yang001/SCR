# Brute-force reference: recompute the potential imbalance W_1(i), W_0(i) of
# Algorithm 1 from scratch (as in SCR 1.0.0) and return P(T_i = 1).
reference_prob <- function(X, codes, tr_prev, d, q1, q2, w) {
  i <- length(tr_prev) + 1L
  n1 <- sum(tr_prev == 1); n0 <- sum(tr_prev == 0)
  if (n1 == 0 || n0 == 0) return(if (n1 == 0 && n0 == 0) 0.5 else if (n1 == 0) 1 else 0)
  if (abs(n1 - n0) >= d) return(if (n1 < n0) q2 else if (n1 > n0) 1 - q2 else 0.5)
  W <- function(t) {
    tt <- c(tr_prev, t)
    val <- 0
    if (!is.null(X)) {
      Xi <- X[seq_len(i), , drop = FALSE]
      u <- colSums(Xi[tt == 1, , drop = FALSE]) - colSums(Xi[tt == 0, , drop = FALSE])
      Mt <- drop(t(u) %*% MASS::ginv(cov(Xi)) %*% u)
      val <- val + (if (is.null(codes)) 1 else w) * Mt
    }
    if (!is.null(codes)) {
      ci <- codes[seq_len(i)]
      m <- max(codes)
      Dn <- tabulate(ci[tt == 1], m) - tabulate(ci[tt == 0], m)
      val <- val + (if (is.null(X)) 1 else 1 - w) * sum(Dn^2)
    }
    val
  }
  W1 <- W(1); W0 <- W(0)
  if (isTRUE(all.equal(W1, W0))) 0.5 else if (W1 < W0) q1 else 1 - q1
}
