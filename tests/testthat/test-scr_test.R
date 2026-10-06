test_that("scr_test returns a valid htest", {
  set.seed(1)
  n <- 60
  x <- matrix(rnorm(n * 2), n, 2)
  tr <- SCR(x)$assignment
  y <- 2 * tr + drop(x %*% c(1, 1)) + rnorm(n)
  out <- scr_test(y, x, tr, B = 50)
  expect_s3_class(out, "htest")
  expect_length(out$t_null, 50)
  expect_true(out$p.value >= 0 && out$p.value <= 1)
  expect_lt(out$p.value, 0.05)
  expect_named(out$estimate, c("difference in means", "regression-adjusted"))
  un <- scr_test(y, x, tr, B = 20, adjust = FALSE)
  expect_equal(unname(un$estimate[1]), mean(y[tr == 1]) - mean(y[tr == 0]))
  z <- data.frame(g = sample(c("A", "B"), n, TRUE))
  expect_s3_class(scr_test(y, x, tr, B = 10, categorical = z), "htest")
})

test_that("the randomization p-value is never below 1 / (B + 1)", {
  set.seed(2)
  x <- matrix(rnorm(80), 40, 2)
  tr <- SCR(x)$assignment
  y <- 10 * tr + rnorm(40)
  out <- scr_test(y, x, tr, B = 30)
  expect_equal(out$p.value, (1 + sum(abs(out$t_null) >= abs(out$statistic))) / 31)
  expect_gte(out$p.value, 1 / 31)
})
