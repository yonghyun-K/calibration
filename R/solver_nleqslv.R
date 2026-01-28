# nleqslv backend (fast dual root-finding)

#' @keywords internal
.weights_from_lambda <- function(lambda, d, X, spec, intercept, w_scale, G_scale) {
  u <- drop(X %*% lambda) / (G_scale * w_scale)
  w <- d / w_scale * .g_inv_impl(u, spec = spec, intercept = intercept)
  w
}

#' @keywords internal
.residual_from_lambda <- function(lambda, d, X, const, spec, intercept, w_scale, G_scale) {
  w <- .weights_from_lambda(lambda, d = d, X = X, spec = spec, intercept = intercept, w_scale = w_scale, G_scale = G_scale)
  if (any(!is.finite(w))) {
    return(rep(Inf, length(const)))
  }
  colSums(X * w) - const
}

#' @keywords internal
.jacobian_from_lambda <- function(lambda, d, X, spec, intercept, w_scale, G_scale, ridge = 0) {
  u <- drop(X %*% lambda) / (G_scale * w_scale)
  gp <- .g_prime_inv_impl(u, spec = spec, intercept = intercept)
  if (any(!is.finite(gp))) {
    J <- matrix(Inf, nrow = ncol(X), ncol = ncol(X))
    return(J)
  }
  J <- t(X) %*% (X * (d / (w_scale^2 * G_scale) * gp))
  if (is.finite(ridge) && ridge > 0) {
    J <- J + diag(ridge, nrow = ncol(X))
  }
  J
}

#' @keywords internal
.solve_nleqslv <- function(X, const, w0, method, spec,
                          w_scale = 1, G_scale = 1,
                          bounds = NULL,
                          solver_control = list()) {

  if (!is.null(bounds)) {
    stop("Bounds are only supported in the 'cvxr' backend for now.")
  }

  n <- nrow(X)
  p <- ncol(X)

  xtol <- if (!is.null(solver_control$xtol)) solver_control$xtol else 1e-10
  maxit <- if (!is.null(solver_control$maxit)) solver_control$maxit else 200
  allowSingular <- if (!is.null(solver_control$allowSingular)) solver_control$allowSingular else TRUE
  ridge <- if (!is.null(solver_control$ridge)) solver_control$ridge else 0
  tol_gap <- if (!is.null(solver_control$constraint_tol)) solver_control$constraint_tol else 1e-8

  # Method-specific setup
  if (identical(method, "BD")) {
    d <- rep(1, n)
    intercept <- .g_impl(w0 * w_scale, spec)
  } else if (identical(method, "GE")) {
    d <- rep(1, n)
    intercept <- rep(0, n)
  } else if (identical(method, "DS")) {
    if (spec$family %in% c("CE", "PH")) {
      stop("method = 'DS' is not implemented for entropy = 'CE' or 'PH'.")
    }
    d <- w0 * w_scale
    r <- if (!is.null(spec$r) && is.finite(spec$r)) spec$r else NA_real_
    intercept <- if (!is.na(r) && r != 0) rep(1 / r, n) else rep(0, n)
    if (!isTRUE(all.equal(w_scale, 1))) {
      warning("w.scale has no effect in method = 'DS' (it cancels analytically).")
    }
  } else {
    stop("Unknown method: ", method)
  }

  # Initial lambda: zeros, with a small nudge for GE on g(w0) column if present
  init <- rep(0, p)
  if (identical(method, "GE")) {
    g_name <- grep("^g\\(.*\\)$", colnames(X), value = TRUE)
    if (length(g_name) == 1) {
      init[which(colnames(X) == g_name)] <- 1
    }
  }

  fn <- function(lambda) {
    .residual_from_lambda(lambda, d = d, X = X, const = const, spec = spec,
                          intercept = intercept, w_scale = w_scale, G_scale = G_scale)
  }

  jac <- function(lambda) {
    .jacobian_from_lambda(lambda, d = d, X = X, spec = spec, intercept = intercept,
                          w_scale = w_scale, G_scale = G_scale, ridge = ridge)
  }

  # Stage 1: attempt solve
  nleqslv_res <- try(
    nleqslv::nleqslv(
      x = init,
      fn = fn,
      jac = jac,
      control = list(maxit = maxit, allowSingular = allowSingular, xtol = xtol),
      xscalm = "auto"
    ),
    silent = TRUE
  )

  diagnostics <- list(
    converged = FALSE,
    solver = "nleqslv",
    message = NA_character_,
    iterations = NA_integer_,
    termcd = NA_integer_,
    max_abs_constraint_gap = Inf,
    residual_norm2 = Inf
  )

  if (inherits(nleqslv_res, "try-error")) {
    diagnostics$message <- as.character(nleqslv_res)
    return(list(w = rep(NA_real_, n), lambda = rep(NA_real_, p), diagnostics = diagnostics, solver_raw = NULL))
  }

  diagnostics$message <- nleqslv_res$message
  diagnostics$termcd <- nleqslv_res$termcd
  diagnostics$iterations <- if (!is.null(nleqslv_res$iter)) nleqslv_res$iter else NA_integer_

  res <- fn(nleqslv_res$x)
  diagnostics$max_abs_constraint_gap <- max(abs(res))
  diagnostics$residual_norm2 <- sqrt(sum(res^2))

  w <- .weights_from_lambda(nleqslv_res$x, d = d, X = X, spec = spec,
                           intercept = intercept, w_scale = w_scale, G_scale = G_scale)

  good_w <- all(is.finite(w))
  close_enough <- is.finite(diagnostics$max_abs_constraint_gap) && diagnostics$max_abs_constraint_gap <= tol_gap

  diagnostics$converged <- isTRUE(close_enough) && isTRUE(good_w)

  if (!diagnostics$converged) {
    w <- rep(NA_real_, n)
  }

  list(w = w, lambda = nleqslv_res$x, diagnostics = diagnostics, solver_raw = nleqslv_res)
}
