# Internal Newton / line-search backend (no external dependencies)

#' @keywords internal
.requires_positive_weights <- function(spec) {
  if (!is.list(spec) || is.null(spec$family)) return(TRUE)

  if (identical(spec$family, "custom")) {
    # Conservative default
    return(TRUE)
  }

  if (identical(spec$family, "CE") || identical(spec$family, "PH")) {
    return(TRUE)
  }

  # Numeric Renyi and string-coded Renyi special cases use spec$r
  r <- spec$r
  if (is.null(r) || !is.finite(r)) return(TRUE)

  # Odd positive integers allow negative weights.
  if (r > 0 && abs(r - round(r)) < .Machine$double.eps && (round(r) %% 2 == 1)) {
    return(FALSE)
  }

  TRUE
}

#' @keywords internal
.v_ok_impl <- function(v, spec) {
  if (identical(spec$family, "custom")) {
    return(rep(TRUE, length(v)))
  }

  if (identical(spec$family, "PH")) {
    del <- spec$del
    return(is.finite(v) & abs(v) < (del - 1e-10))
  }

  if (identical(spec$family, "CE")) {
    return(is.finite(v) & (v < -1e-10))
  }

  # Renyi (numeric or mapped) cases
  r <- spec$r
  if (is.null(r) || !is.finite(r)) {
    return(is.finite(v))
  }

  if (r == 0) {
    return(is.finite(v))
  }

  if (r > 0 && abs(r - round(r)) < .Machine$double.eps && (round(r) %% 2 == 1)) {
    # odd positive integer: real powers defined for negative too
    return(is.finite(v))
  }

  # require r*v > 0
  is.finite(v) & (r * v > 1e-12)
}

#' @keywords internal
.w_ok_impl <- function(w, spec) {
  if (identical(spec$family, "custom")) {
    return(all(is.finite(w)))
  }

  if (identical(spec$family, "CE")) {
    return(all(is.finite(w) & w > 1))
  }

  if (identical(spec$family, "PH")) {
    return(all(is.finite(w) & w > 0))
  }

  r <- spec$r
  if (is.null(r) || !is.finite(r)) {
    return(all(is.finite(w)))
  }

  if (r == 0) {
    return(all(is.finite(w) & w > 0))
  }

  if (r > 0 && abs(r - round(r)) < .Machine$double.eps && (round(r) %% 2 == 1)) {
    return(all(is.finite(w)))
  }

  all(is.finite(w) & w > 0)
}

#' @keywords internal
.solve_newton <- function(X, const, w0, method, spec,
                          w_scale = 1, G_scale = 1,
                          bounds = NULL,
                          solver_control = list()) {

  n <- nrow(X)
  p <- ncol(X)

  # Parse bounds into length-n vectors (NULL means no bound on that side).
  # We reuse the same parser as the CVXR backend so the user-facing API is consistent.
  bnd <- .parse_bounds(bounds, n)
  lower <- bnd$lower
  upper <- bnd$upper
  has_bounds <- !is.null(lower) || !is.null(upper)

  maxit <- if (!is.null(solver_control$maxit)) solver_control$maxit else 50
  step_halving <- if (!is.null(solver_control$step_halving)) solver_control$step_halving else 12
  ridge <- if (!is.null(solver_control$ridge)) solver_control$ridge else 1e-8
  tol_gap <- if (!is.null(solver_control$constraint_tol)) solver_control$constraint_tol else 1e-8
  armijo_c <- if (!is.null(solver_control$armijo_c)) solver_control$armijo_c else 1e-4

  # Method-specific setup (match the nleqslv backend)
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

  if (any(!is.finite(intercept))) {
    diagnostics <- list(
      converged = FALSE,
      solver = "newton",
      iterations = 0L,
      message = "Non-finite intercept; check that w0 lies in the divergence domain.",
      max_abs_constraint_gap = Inf,
      residual_norm2 = Inf
    )
    return(list(w = rep(NA_real_, n), lambda = rep(NA_real_, p), diagnostics = diagnostics, solver_raw = NULL))
  }

  # Initial lambda
  lambda <- rep(0, p)
  if (identical(method, "GE")) {
    g_name <- grep("^g\\(.*\\)$", colnames(X), value = TRUE)
    if (length(g_name) == 1) {
      lambda[which(colnames(X) == g_name)] <- 1
    }
  }

  # ------------------------------------------------------------
  # Bounds via duality (Huberized / truncated conjugate)
  # ------------------------------------------------------------
  # With bounds M1 <= w_i <= M2, the dual remains an unconstrained minimization
  # in lambda and the primal weights are obtained by the clipped inverse-gradient
  # map: w*(lambda) = clip(g^{-1}(nu_i(lambda)), M1, M2).
  # See the attached derivations in tmp2.pdf.

  # Pre-compute bound thresholds in the "omega_tilde" scale, where
  #   omega_tilde_i = g^{-1}(intercept_i + u_i)
  # and
  #   w_i = (d_i / w_scale) * omega_tilde_i.
  # So bounds on w translate to bounds on omega_tilde via:
  #   omega_tilde_i ∈ [lower_i * w_scale / d_i, upper_i * w_scale / d_i].
  low_tilde <- up_tilde <- v_low <- v_up <- NULL
  if (has_bounds) {
    # If d is non-positive, the transformation is ill-defined.
    if (any(!is.finite(d)) || any(d <= 0)) {
      stop("Bounds require positive scaling (d > 0) in the Newton backend.")
    }

    # Warn (softly) if initial weights fall outside the requested bounds.
    # The duality derivation based on truncating the generator typically assumes
    # w0 is within [lower, upper]. We still proceed but surface a message.
    if (!is.null(lower) && any(is.finite(lower) & (w0 < lower))) {
      warning("Some entries of w0 are below the requested lower bound; bounded calibration may be ill-posed.", call. = FALSE)
    }
    if (!is.null(upper) && any(is.finite(upper) & (w0 > upper))) {
      warning("Some entries of w0 are above the requested upper bound; bounded calibration may be ill-posed.", call. = FALSE)
    }

    low_tilde <- if (!is.null(lower)) lower * w_scale / d else rep(-Inf, n)
    up_tilde  <- if (!is.null(upper)) upper * w_scale / d else rep( Inf, n)

    # Sanity check after transforming bounds.
    bad <- is.finite(low_tilde) & is.finite(up_tilde) & (low_tilde > up_tilde)
    if (any(bad)) {
      stop("Transformed bounds are inconsistent (some lower > upper after scaling).",
           " This can happen if bounds are not compatible with w0/method scaling.")
    }

    # Thresholds in nu-space: g(M1) and g(M2). Outside this range, the
    # clipped map saturates at the corresponding bound (no need to evaluate g^{-1}).
    v_low <- rep(-Inf, n)
    if (!is.null(lower)) {
      idx <- is.finite(low_tilde)
      if (any(idx)) {
        tmp <- .g_impl(low_tilde[idx], spec)
        if (any(!is.finite(tmp))) {
          stop("Lower bound violates the divergence domain for the chosen entropy.")
        }
        v_low[idx] <- tmp
      }
    }

    v_up <- rep(Inf, n)
    if (!is.null(upper)) {
      idx <- is.finite(up_tilde)
      if (any(idx)) {
        tmp <- .g_impl(up_tilde[idx], spec)
        if (any(!is.finite(tmp))) {
          stop("Upper bound violates the divergence domain for the chosen entropy.")
        }
        v_up[idx] <- tmp
      }
    }
  }

  # Compute weights and the derivative of omega_tilde with respect to u.
  # The derivative is 0 for saturated units and equals (g^{-1})'(nu) for interior units.
  weight_map <- function(lam) {
    u <- drop(X %*% lam) / (G_scale * w_scale)
    nu <- intercept + u

    omega_tilde <- rep(NA_real_, n)
    gp <- rep(0, n)
    interior <- rep(TRUE, n)

    if (has_bounds) {
      if (!is.null(lower)) {
        idx <- is.finite(low_tilde) & (nu <= v_low)
        if (any(idx)) {
          omega_tilde[idx] <- low_tilde[idx]
          gp[idx] <- 0
          interior[idx] <- FALSE
        }
      }

      if (!is.null(upper)) {
        idx <- is.finite(up_tilde) & (nu >= v_up)
        if (any(idx)) {
          omega_tilde[idx] <- up_tilde[idx]
          gp[idx] <- 0
          interior[idx] <- FALSE
        }
      }
    }

    if (any(interior)) {
      omega_int <- .g_inv_impl(u[interior], spec = spec, intercept = intercept[interior])
      gp_int <- .g_prime_inv_impl(u[interior], spec = spec, intercept = intercept[interior])
      if (any(!is.finite(omega_int)) || any(!is.finite(gp_int))) {
        return(list(ok = FALSE))
      }
      omega_tilde[interior] <- omega_int
      gp[interior] <- gp_int
    }

    if (any(!is.finite(omega_tilde))) {
      return(list(ok = FALSE))
    }

    w <- d / w_scale * omega_tilde

    # Numerical guard: enforce the provided (finite) bounds directly on w.
    if (has_bounds) {
      if (!is.null(lower)) {
        idx <- is.finite(lower)
        if (any(idx)) w[idx] <- pmax(w[idx], lower[idx])
      }
      if (!is.null(upper)) {
        idx <- is.finite(upper)
        if (any(idx)) w[idx] <- pmin(w[idx], upper[idx])
      }
    }

    if (!.w_ok_impl(w, spec)) {
      return(list(ok = FALSE))
    }

    list(ok = TRUE, w = w, gp = gp)
  }

  # Main iterations
  converged <- FALSE
  message <- NA_character_
  res <- rep(Inf, p)

  for (it in seq_len(maxit)) {
    mw <- weight_map(lambda)
    if (!isTRUE(mw$ok)) {
      message <- "Weights left the divergence domain (or bounds/domain were incompatible) during Newton iterations."
      break
    }

    w <- mw$w
    res <- colSums(X * w) - const
    max_gap <- max(abs(res))

    if (is.finite(max_gap) && max_gap <= tol_gap) {
      converged <- TRUE
      message <- "Converged."
      break
    }

    # Jacobian: t(X) %*% diag(d/(w_scale^2*G_scale) * gp) %*% X
    fac <- d / (w_scale^2 * G_scale) * mw$gp
    if (any(!is.finite(fac))) {
      message <- "Non-finite Jacobian encountered during Newton iterations."
      break
    }

    J <- t(X) %*% (X * fac)
    if (is.finite(ridge) && ridge > 0) {
      J <- J + diag(ridge, p)
    }

    step_dir <- tryCatch(solve(J, res), error = function(e) qr.solve(J, res))

    g0 <- sum(res^2)
    step <- 1
    accepted <- FALSE

    for (ls in seq_len(step_halving)) {
      lam_new <- lambda - step * step_dir
      mw_new <- weight_map(lam_new)
      if (isTRUE(mw_new$ok)) {
        res_new <- colSums(X * mw_new$w) - const
        g1 <- sum(res_new^2)
        if (is.finite(g1) && g1 <= (1 - armijo_c * step) * g0) {
          lambda <- lam_new
          accepted <- TRUE
          break
        }
      }
      step <- step * 0.5
    }

    if (!accepted) {
      # Conservative fallback, similar to the reference implementation
      lambda <- lambda * 0.5
    }

    if (it == maxit) {
      message <- "Maximum iterations reached without convergence."
    }
  }

  # Final evaluation
  mw_final <- weight_map(lambda)
  w_final <- if (isTRUE(mw_final$ok)) mw_final$w else rep(NA_real_, n)
  res_final <- if (all(is.finite(w_final))) colSums(X * w_final) - const else rep(Inf, p)
  max_gap_final <- max(abs(res_final))

  diagnostics <- list(
    converged = isTRUE(converged) && all(is.finite(w_final)) && is.finite(max_gap_final) && (max_gap_final <= tol_gap),
    solver = "newton",
    iterations = if (exists("it")) as.integer(it) else NA_integer_,
    message = message,
    max_abs_constraint_gap = max_gap_final,
    residual_norm2 = sqrt(sum(res_final^2))
  )

  if (has_bounds && isTRUE(diagnostics$converged)) {
    btol <- if (!is.null(solver_control$bound_tol)) solver_control$bound_tol else 1e-8
    diagnostics$n_at_lower <- if (!is.null(lower)) {
      sum(is.finite(lower) & abs(w_final - lower) <= btol)
    } else 0L
    diagnostics$n_at_upper <- if (!is.null(upper)) {
      sum(is.finite(upper) & abs(w_final - upper) <= btol)
    } else 0L
  }

  if (!diagnostics$converged) {
    w_final <- rep(NA_real_, n)
    lambda <- rep(NA_real_, p)
  }

  list(w = w_final, lambda = lambda, diagnostics = diagnostics, solver_raw = NULL)
}
