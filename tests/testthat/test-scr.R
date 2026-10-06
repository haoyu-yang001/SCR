test_that("SCR returns a valid 0/1 allocation", {
  set.seed(1)
  x <- matrix(rnorm(60 * 3), 60, 3)
  fit <- SCR(x)
  expect_s3_class(fit, "SCR")
  expect_length(fit$assignment, 60)
  expect_true(all(fit$assignment %in% c(0L, 1L)))
  expect_false(anyNA(fit$assignment))
  expect_equal(sort(fit$assignment[1:2]), c(0L, 1L))
  expect_equal(fit$rule[1:2], c("burn-in", "burn-in"))
  expect_equal(unname(fit$sample_size), c(sum(fit$assignment), sum(1 - fit$assignment)))
  expect_output(print(fit), "Sequential Covariate-Adjusted Randomization")
})

test_that("allocation probabilities match a brute-force implementation of Algorithm 1", {
  set.seed(2)
  n <- 80
  x <- cbind(rnorm(n, 50, 10), rexp(n), rnorm(n))
  for (d in c(1, 3, 5, n)) {
    fit <- SCR(x, d = d, q1 = 0.75, q2 = 0.85)
    for (i in 3:n) {
      ref <- reference_prob(x, NULL, fit$assignment[seq_len(i - 1)], d, 0.75, 0.85, 0.7)
      expect_equal(fit$prob[i], ref, info = sprintf("d = %s, unit %d", d, i))
    }
  }
})

test_that("mixed covariates follow W(n) = w M + (1 - w) ||D||^2", {
  set.seed(3)
  n <- 80
  x <- matrix(rnorm(n * 2), n, 2)
  z <- data.frame(a = sample(letters[1:2], n, TRUE), b = sample(1:3, n, TRUE))
  codes <- SCR:::.strata_codes(z)
  fit <- SCR(x, categorical = z, d = 6, w = 0.4)
  for (i in 3:n) {
    ref <- reference_prob(x, codes, fit$assignment[seq_len(i - 1)], 6, 0.75, 0.85, 0.4)
    expect_equal(fit$prob[i], ref, info = sprintf("unit %d", i))
  }
  only <- SCR(NULL, categorical = z, d = 6)
  for (i in 3:n) {
    ref <- reference_prob(NULL, codes, only$assignment[seq_len(i - 1)], 6, 0.75, 0.85, 0.7)
    expect_equal(only$prob[i], ref, info = sprintf("unit %d", i))
  }
})

test_that("SCR() and repeated scr_next() give identical assignments", {
  set.seed(4)
  x <- matrix(rnorm(50 * 4), 50, 4)
  z <- data.frame(g = sample(c("A", "B"), 50, TRUE))
  set.seed(10)
  all_at_once <- SCR(x, categorical = z, d = 4)$assignment
  set.seed(10)
  tr <- integer(0)
  for (i in 1:50) {
    prev <- seq_len(i - 1)
    tr <- c(tr, scr_next(x[i, ], x[prev, , drop = FALSE], tr, d = 4,
                         categorical_new = z[i, , drop = FALSE],
                         categorical = z[prev, , drop = FALSE])$assignment)
  }
  expect_identical(all_at_once, tr)
})

test_that("prior assignments are kept", {
  set.seed(5)
  x <- matrix(rnorm(30 * 2), 30, 2)
  fit <- SCR(x, assignment = c(1, 1, 0, 1))
  expect_equal(fit$assignment[1:4], c(1L, 1L, 0L, 1L))
  expect_equal(fit$rule[1:4], rep("given", 4))
  expect_true(all(is.na(fit$prob[1:4])))
  expect_false(anyNA(fit$prob[-(1:4)]))
  one <- SCR(x, assignment = 1)
  expect_equal(one$assignment[2], 0L)
  expect_identical(SCR(x, assignment = NA)$parameters, SCR(x)$parameters)
})

test_that("marginal step and the constraint d behave as specified", {
  set.seed(6)
  x <- matrix(rnorm(200 * 3), 200, 3)
  # q2 = 1 forces the smaller group, so |n1 - n0| never exceeds d.
  for (d in c(1, 2, 5)) {
    tr <- SCR(x, d = d, q2 = 1)$assignment
    expect_lte(max(abs(cumsum(2 * tr - 1))), d)
  }
  # d = 0: never covariate-adaptive, no missing assignments (bug in 1.0.0).
  fit0 <- SCR(x, d = 0)
  expect_false(anyNA(fit0$assignment))
  expect_true(all(fit0$rule[-(1:2)] == "marginal"))
  bal <- cumsum(2 * fit0$assignment - 1)[2:199]  # n1 - n0 before units 3..200
  expect_equal(fit0$prob[3:200], ifelse(bal == 0, 0.5, ifelse(bal < 0, 0.85, 1 - 0.85)))
  # d = n: covariate step only.
  expect_true(all(SCR(x, d = 200)$rule[-(1:2)] == "covariate"))
})

test_that("ties are randomized with probability 0.5", {
  x <- rbind(c(1, 2), c(1, 2), matrix(rnorm(20), 10, 2))
  fit <- SCR(x)
  expect_equal(fit$prob[3], 0.5)
})

test_that("SCR balances covariates far better than complete randomization", {
  set.seed(7)
  n <- 300; p <- 5
  m_scr <- m_cr <- numeric(20)
  for (r in 1:20) {
    x <- matrix(rnorm(n * p), n, p)
    m_scr[r] <- SCR(x)$Mahalanobis_Distance
    m_cr[r] <- scr_imbalance(x, rbinom(n, 1, 0.5))[["Mahalanobis_Distance"]]
  }
  expect_lt(mean(m_scr), mean(m_cr) / 5)
})

test_that("input errors are informative", {
  x <- matrix(rnorm(20), 10, 2)
  expect_error(SCR(x, assignment = c(1, 2, 2)), "assignment - 1")
  expect_error(SCR(data.frame(a = 1:3, b = c("u", "v", "w"))), "categorical")
  expect_error(SCR(x, q1 = 0.3), "q1")
  expect_error(SCR(x, d = -1), "`d`")
  expect_error(SCR(NULL), "Supply")
  expect_error(scr_next(c(0, 0), x, assignment = c(0, 1)), "one value")
})
