test_that("smart allocates to K arms within the marginal constraint", {
  set.seed(1)
  x <- matrix(rnorm(120 * 3), 120, 3)
  for (method in c("mean", "max", "median")) {
    fit <- smart(x, K = 3, d = 4, method = method)
    expect_s3_class(fit, "smart")
    expect_true(all(fit$assignment %in% 1:3))
    expect_setequal(fit$assignment[1:3], 1:3)
    counts <- sapply(1:3, function(k) cumsum(fit$assignment == k))
    expect_lte(max(apply(counts[-(1:3), ], 1, function(r) max(r) - min(r))), 4)
    expect_equal(unname(rowSums(fit$prob[-(1:3), ])), rep(1, 117))
  }
  expect_length(smart(x, K = 4)$pairwise_Mahalanobis, 6)
  expect_output(print(smart(x, K = 3)), "K = 3")
})

test_that("smart keeps prior assignments and checks input", {
  set.seed(2)
  x <- matrix(rnorm(30 * 2), 30, 2)
  fit <- smart(x, assignment = c(3, 3), K = 3)
  expect_equal(fit$assignment[1:2], c(3L, 3L))
  expect_setequal(fit$assignment[3:4], 1:2)
  expect_error(smart(x, K = 1), "K")
  expect_error(smart(x, assignment = c(0, 1), K = 3), "values")
})
