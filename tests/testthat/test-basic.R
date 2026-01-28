
test_that("basic DS calibration works and satisfies constraints", {
  Xs <- cbind(
    rep(1, 10),
    c(rep(1, 5), rep(0, 5)),
    1:10
  )
  piks <- rep(0.2, 10)
  d <- 1 / piks
  total <- c(50, 24, 290)

  fit <- calibrate(~ 0 + Xs, w0 = d, const = total,
                   method = "DS", entropy = "SL",
                   solver = "newton",
                   solver_control = list(constraint_tol = 1e-6))

  expect_true(isTRUE(fit$diagnostics$converged))
  expect_true(all(is.finite(weights(fit))))
  gap <- colSums(fit$X * weights(fit)) - fit$const
  expect_lt(max(abs(gap)), 1e-6)
})

test_that("const with NA is rejected", {
  Xs <- cbind(rep(1, 3), 1:3)
  d <- rep(5, 3)
  expect_error(
    calibrate(~ 0 + Xs, w0 = d, const = c(1, NA), method = "DS", entropy = "SL"),
    "must not contain NA"
  )
})

test_that("named const must match model.matrix column names", {
  Xs <- cbind(rep(1, 3), 1:3)
  d <- rep(5, 3)
  expect_error(
    calibrate(~ 0 + Xs, w0 = d, const = c(wrong = 1, also_wrong = 2), method = "DS", entropy = "SL"),
    "missing required names"
  )
})

test_that("factor outcomes are supported via indicator expansion", {
  df <- data.frame(
    x = factor(c("A", "A", "B", "B", "B", "A")),
    y = factor(c("L1", "L2", "L1", "L2", "L2", "L3"), levels = c("L1", "L2", "L3"))
  )

  # Trivial constraints so that w == w0 == 1 is feasible and should converge easily.
  X <- model.matrix(~ x, data = df)
  const <- colSums(X)
  names(const) <- colnames(X)

  fit <- calibrate(y ~ x, w0 = 1, data = df, const = const,
                   method = "BD", entropy = "SL",
                   solver = "newton")

  expect_true(isTRUE(fit$diagnostics$converged))
  expect_true(!is.null(fit$estimate))

  est <- fit$estimate[, "Estimate"]
  expect_equal(unname(est["y:L1"]), 2)
  expect_equal(unname(est["y:L2"]), 3)
  expect_equal(unname(est["y:L3"]), 1)
})

test_that("multivariate numeric outcomes work via cbind()", {
  df <- data.frame(
    x = factor(c("A", "A", "B", "B")),
    y1 = c(1, 2, 3, 4),
    y2 = c(10, 20, 30, 40)
  )

  X <- model.matrix(~ x, data = df)
  const <- colSums(X)
  names(const) <- colnames(X)

  fit <- calibrate(cbind(y1, y2) ~ x, w0 = 1, data = df, const = const,
                   method = "BD", entropy = "SL",
                   solver = "newton")

  expect_true(isTRUE(fit$diagnostics$converged))
  est <- fit$estimate[, "Estimate"]
  expect_equal(unname(est["y1"]), sum(df$y1))
  expect_equal(unname(est["y2"]), sum(df$y2))
})

test_that("~ 0 (no constraints) returns w0 unchanged", {
  df <- data.frame(
    y = c(1, 2, 3, 4),
    w0 = c(2, 3, 4, 5)
  )

  fit <- calibrate(~ 0, w0 = w0, data = df, const = numeric(0),
                   method = "DS", entropy = "SL")

  expect_true(isTRUE(fit$diagnostics$converged))
  expect_equal(weights(fit), df$w0)
  expect_equal(ncol(fit$X), 0)
})

test_that("mean-scale (sum(w)=1) variance does not produce NaN", {
  df <- data.frame(
    y = c(1, 2, 3, 4),
    w0 = c(2, 3, 4, 5)
  )

  fit <- calibrate(y ~ 1, w0 = w0, data = df,
                   const = c("(Intercept)" = 1),
                   method = "DS", entropy = "SL")

  expect_true(isTRUE(fit$diagnostics$converged))
  se <- fit$estimate[, "Std. Error"]
  expect_true(all(is.finite(se)))
  expect_false(any(is.nan(se)))
})
