# Main user-facing API

#' @keywords internal
.has_finite_bounds <- function(bounds) {
  if (is.null(bounds)) return(FALSE)

  if (is.numeric(bounds) && length(bounds) == 2) {
    # Treat c(-Inf, Inf) as "no bounds"
    return(any(is.finite(bounds)))
  }

  if (is.list(bounds)) {
    vals <- c(bounds$lower, bounds$upper)
    return(any(is.finite(vals)))
  }

  if (is.matrix(bounds) && ncol(bounds) == 2) {
    return(any(is.finite(bounds)))
  }

  # Let downstream parsing throw a clear error for unsupported forms.
  TRUE
}


#' Calibrate survey weights
#'
#' `calibrate()` computes calibration weights under linear constraints, and (optionally)
#' computes estimated totals and standard errors if the formula has a response.
#'
#' @param formula A formula. The right-hand side defines calibration constraints.
#'   If a response is provided on the left-hand side, the calibrated-weight total
#'   (and a linearized standard error) is returned in one call.
#' @param w0 Design weights (vector) or a scalar. If `w0` is a name and `data` is
#'   provided, it is looked up in `data`.
#' @param data Optional data frame.
#' @param const Numeric vector of known population totals. For safety, a *named*
#'   vector is recommended; names are matched to `colnames(model.matrix(...))`.
#' @param method Calibration method: "BD" (Bregman divergence), "GE" (generalized entropy),
#'   or "DS" (Deville-Sarndal style divergence).
#' @param entropy Entropy family code ("SL", "EL", "ET", "HD", "CE", "PH"),
#'   a numeric Renyi order, or a list with `code`/`family` and (optionally) `del`
#'   (e.g., `list(code = "PH", del = 0.5)`).
#' @param divergence Optional custom divergence specification. If supplied, this overrides `entropy`.
#'   Required functions depend on the solver: for `solver = "cvxr"`, supply `G(x)`; for
#'   `solver = "nleqslv"`, supply `g_inv(u, intercept = 0)`; for `solver = "newton"`/`"auto"`,
#'   supply `G(x)` and `g_inv(u, intercept = 0)`. Optional functions `g`, `g_prime_inv`,
#'   and `fprime` improve accuracy and avoid numeric fallbacks.
#' @param w.scale Weight scaling factor (phi). See package documentation.
#' @param G.scale Entropy scaling factor (q). See package documentation.
#' @param bounds Optional weight bounds.
#'   Can be a length-2 numeric vector `(lower, upper)`, a list `list(lower=, upper=)`,
#'   or an `n x 2` matrix.
#'
#'   When `solver = "auto"` (the default) or `solver = "newton"`, bounds are handled
#'   via the dual formulation: the weight map becomes a *clipped* inverse-gradient
#'   (piecewise) update, which avoids introducing separate inequality multipliers.
#'   See the duality derivation in the attached notes. 
#' @param solver Solver backend: "auto", "nleqslv", or "cvxr".
#' @param solver_control List of backend controls. Common entries include
#'   `maxit`, `xtol`, `ridge`. For `solver = "cvxr"`, you may also provide
#'   `cvxr_solver` (e.g., `"MOSEK"`, `"ECOS_BB"`), `cvxr_opts` (a list of options
#'   passed to [CVXR::solve()]), `bounds_tol` (tolerance for bound violations),
#'   and `domain_eps` (positivity floor).
#'
#' @return A `calibration_fit` object. Use [weights()] to extract weights.
#' @export
calibrate <- function(formula,
                      w0 = 1,
                      data = NULL,
                      const,
                      method = c("BD", "GE", "DS"),
                      entropy = c("EL", "ET", "HD", "CE", "PH", "SL"),
                      del = NULL,
                      divergence = NULL,
                      w.scale = 1,
                      G.scale = 1,
                      bounds = NULL,
                      solver = c("auto", "newton", "nleqslv", "cvxr"),
                      solver_control = list(),
                      ...) {

  method <- match.arg(method)
  solver <- match.arg(solver)

  if (!inherits(formula, "formula")) stop("'formula' must be a formula.")

  # Evaluate w0
  w0_expr <- substitute(w0)
  w0_name <- if (is.name(w0_expr)) deparse(w0_expr) else "w0"
  if (is.null(data)) {
    w0_val <- eval(w0_expr, envir = parent.frame())
  } else {
    if (is.name(w0_expr) && w0_name %in% names(data)) {
      w0_val <- data[[w0_name]]
    } else {
      w0_val <- eval(w0_expr, envir = data, enclos = parent.frame())
    }
  }

  # Recycle scalar w0
  if (length(w0_val) == 1) {
    # We'll expand after we know n
    w0_is_scalar <- TRUE
  } else {
    w0_is_scalar <- FALSE
  }

  if (!is.null(del)) {
    warning(
      "`del` is deprecated; supply `entropy = list(code = \"PH\", del = ...)` instead.",
      call. = FALSE
    )
  }

  # Normalize entropy / divergence spec
  solver_for_validation <- if (identical(solver, "auto")) "newton" else solver
  spec <- if (!is.null(divergence)) {
    .validate_custom_divergence(divergence, solver = solver_for_validation)
  } else {
    entropy_spec(entropy, del = del)
  }

  # Build an evaluation env for model.frame/model.matrix that:
  # * injects a local g() closure bound to (spec)
  # * injects w0 into the formula environment so users can write g(w0) in the formula
  eval_env <- new.env(parent = environment(formula))

  # Hint for constant terms: if users write g(1) in a formula, model.frame() expects
  # a length-n vector. We'll recycle scalar outputs to nrow(data) when 'data' is provided.
  eval_env$.__n <- if (!is.null(data)) NROW(data) else NA_integer_

  eval_env$g <- function(x) {
    out <- .g_impl(x, spec)
    if (length(out) == 1L && is.finite(eval_env$.__n) && eval_env$.__n > 1L) {
      out <- rep(out, eval_env$.__n)
    }
    out
  }

  # Put w0 into the evaluation environment
  # (only if it's a valid name)
  if (is.name(w0_expr) && nzchar(w0_name) && make.names(w0_name) == w0_name) {
    assign(w0_name, w0_val, envir = eval_env)
  } else {
    assign("w0", w0_val, envir = eval_env)
    w0_name <- "w0"
  }

  environment(formula) <- eval_env

  # Build model frame and design matrix for calibration constraints
  if (is.null(data)) {
    mf <- stats::model.frame(formula, na.action = stats::na.fail)
  } else {
    mf <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
  }

  X <- stats::model.matrix(stats::terms(mf), mf)

  n <- nrow(X)
  if (w0_is_scalar) w0_val <- rep(as.numeric(w0_val), n)
  if (length(w0_val) != n) {
    stop("Length of w0 (", length(w0_val), ") does not match nrow(model.matrix)=", n)
  }

  # Align constraints by name (recommended) or position (fallback)
  const_aligned <- .match_const(const, X)

  # Detect whether the user supplied an outcome on the left-hand side.
  # (Used for the one-call workflow.)
  has_response <- length(formula) >= 3 && !is.null(formula[[2]])

  # Special case: no calibration constraints (p = 0), e.g., formula is `~ 0`.
  # In this case, the calibrated weights are simply w0.
  if (ncol(X) == 0) {
    if (length(const_aligned) != 0) {
      stop("When the calibration formula produces no constraint columns (e.g., '~ 0'), 'const' must be length 0.")
    }
    if (identical(method, "GE")) {
      stop("method='GE' requires a g(w0) term in the formula; formula '~ 0' has no constraint columns.")
    }

    w_out <- as.numeric(w0_val)

    # If bounds are supplied, the bounded solution with p=0 is the projection of w0
    # onto the box constraints.
    if (.has_finite_bounds(bounds)) {
      bnd0 <- .parse_bounds(bounds, n)
      if (!is.null(bnd0$lower)) {
        idx <- is.finite(bnd0$lower)
        if (any(idx)) w_out[idx] <- pmax(w_out[idx], bnd0$lower[idx])
      }
      if (!is.null(bnd0$upper)) {
        idx <- is.finite(bnd0$upper)
        if (any(idx)) w_out[idx] <- pmin(w_out[idx], bnd0$upper[idx])
      }
    }

    fit0 <- list(
      call = match.call(),
      formula = formula,
      data_provided = !is.null(data),
      method = method,
      entropy_spec = spec,
      w0_name = w0_name,
      w0 = w0_val,
      w = w_out,
      lambda = numeric(0),
      X = X,
      const = const_aligned,
      w_scale = w.scale,
      G_scale = G.scale,
      bounds = bounds,
      solver = "none",
      diagnostics = list(
        converged = TRUE,
        solver = "none",
        iterations = 0L,
        message = if (.has_finite_bounds(bounds)) {
          "No calibration constraints (p = 0): returning w0 after applying bounds."
        } else {
          "No calibration constraints (p = 0): returning w0 unchanged."
        },
        max_abs_constraint_gap = 0,
        residual_norm2 = 0,
        solver_status = NA_character_,
        status = "solved"
      ),
      solver_raw = NULL,
      constraints = .constraint_table(X, as.numeric(w0_val), const_aligned),
      estimate = NULL,
      cov = NULL
    )
    class(fit0) <- "calibration_fit"

    if (has_response) {
      est0 <- estimate(formula, data = data, calibration = fit0)
      fit0$estimate <- est0$estimate
      fit0$cov <- est0$cov
    }

    return(fit0)
  }

  # For method GE, enforce presence of g(w0) column and apply the scaling implied by q*phi*g(phi*d)
  if (identical(method, "GE")) {
    g_col_name <- paste0("g(", w0_name, ")")
    if (!(g_col_name %in% colnames(X))) {
      stop("method='GE' requires a g(w0) term in the formula; expected column name: ", g_col_name)
    }
    # Replace the computed column with q*phi*g(phi*w0)
    X[, g_col_name] <- G.scale * w.scale * .g_impl(w0_val * w.scale, spec)
  }

  # Solver selection
  if (identical(solver, "auto")) {
    # Prefer the internal Newton solver by default (fast, no external dependencies).
    # Bounds are handled via the dual clipped inverse-gradient map in the Newton backend.
    solver <- "newton"
  }

  # Stage 1 solve
  sol <- if (identical(solver, "newton")) {
    .solve_newton(X = X, const = const_aligned, w0 = w0_val, method = method, spec = spec,
                  w_scale = w.scale, G_scale = G.scale, bounds = bounds, solver_control = solver_control)
  } else if (identical(solver, "nleqslv")) {
    if (!requireNamespace("nleqslv", quietly = TRUE)) {
      stop("solver='nleqslv' requested, but package 'nleqslv' is not installed. Install it or use solver='newton'.")
    }
    .solve_nleqslv(X = X, const = const_aligned, w0 = w0_val, method = method, spec = spec,
                  w_scale = w.scale, G_scale = G.scale, bounds = bounds, solver_control = solver_control)
  } else {
    .solve_cvxr(X = X, const = const_aligned, w0 = w0_val, method = method, spec = spec,
                w_scale = w.scale, G_scale = G.scale, bounds = bounds, solver_control = solver_control)
  }

  diagnostics <- sol$diagnostics

  # Keep backend-specific status (e.g., CVXR's sol$status) separate from the
  # package-level status used for UX and feasibility reporting.
  if (is.null(diagnostics$solver_status)) diagnostics$solver_status <- NA_character_

  # Stage 2 diagnostics if needed
  if (!isTRUE(diagnostics$converged)) {
    feas <- .diagnose_feasibility(X, const_aligned, require_positive = .requires_positive_weights(spec))
    diagnostics$feasibility <- feas
    diagnostics$status <- if (isTRUE(feas$likely_infeasible)) "infeasible" else "solver_failed"
  } else {
    diagnostics$status <- "solved"
  }

  # Warn if NA weights are produced (requested user experience behavior)
  if (!isTRUE(diagnostics$converged)) {
    reasons <- NULL
    if (!is.null(diagnostics$feasibility$reasons) && length(diagnostics$feasibility$reasons) > 0) {
      reasons <- diagnostics$feasibility$reasons[1]
    }

    solver_status_txt <- if (!is.null(diagnostics$solver_status) && !is.na(diagnostics$solver_status)) {
      paste0(", solver_status='", diagnostics$solver_status, "'")
    } else {
      ""
    }

    msg_txt <- if (!is.null(diagnostics$message) && !is.na(diagnostics$message) && nzchar(diagnostics$message)) {
      if (!is.null(diagnostics$solver_status) && !is.na(diagnostics$solver_status) &&
          identical(as.character(diagnostics$message), as.character(diagnostics$solver_status))) {
        ""
      } else {
        paste0(" ", diagnostics$message)
      }
    } else {
      ""
    }

    warning(
      paste0(
        "Calibration did not converge (solver='", solver, "', status='", diagnostics$status, "'). ",
        "Returning NA weights.",
        solver_status_txt,
        msg_txt,
        " ",
        if (!is.null(reasons)) paste0("Hint: ", reasons) else ""
      ),
      call. = FALSE
    )
  }

  fit <- list(
    call = match.call(),
    formula = formula,
    data_provided = !is.null(data),
    method = method,
    entropy_spec = spec,
    w0_name = w0_name,
    w0 = w0_val,
    w = sol$w,
    lambda = sol$lambda,
    X = X,
    const = const_aligned,
    w_scale = w.scale,
    G_scale = G.scale,
    bounds = bounds,
    solver = solver,
    diagnostics = diagnostics,
    solver_raw = sol$solver_raw
  )

  if (isTRUE(diagnostics$converged)) {
    fit$constraints <- .constraint_table(X, fit$w, fit$const)
  } else {
    fit$constraints <- NULL
  }

  # One-call workflow: if response exists, compute totals & SEs
  if (has_response && isTRUE(diagnostics$converged)) {
    est <- estimate(formula, data = data, calibration = structure(fit, class = "calibration_fit"))
    fit$estimate <- est$estimate
    fit$cov <- est$cov
  }

  class(fit) <- "calibration_fit"
  fit
}

#' Calibrate survey weights (alias)
#'
#' This is a simple alias for [calibrate()] that avoids name clashes with
#' `survey::calibrate()` in interactive sessions.
#'
#' @inheritParams calibrate
#' @return A `calibration_fit` object.
#' @export
calibrate_weights <- function(...) {
  calibrate(...)
}

#' S3 method for `calibrate()` when the `survey` package masks the name
#'
#' The `survey` package defines `calibrate()` as an S3 generic that dispatches
#' on its first argument (usually a survey design object). If `survey` is loaded
#' *after* this package, `survey::calibrate()` can mask [calibrate()].
#'
#' Registering `calibrate.formula()` allows calls like `calibrate(~ x, ...)`
#' to still work, by dispatching to this package's implementation.
#'
#' @param x A formula.
#' @param ... Passed to [calibrate()].
#' @export
calibrate.formula <- function(x, ...) {
  calibrate(formula = x, ...)
}

# ---- S3 methods --------------------------------------------------------------

#' @export
print.calibration_fit <- function(x, ...) {
  cat("Calibration fit\n")
  cat("  method:", x$method, "\n")
  cat("  solver:", x$solver, "\n")
  cat("  converged:", isTRUE(x$diagnostics$converged), " (", x$diagnostics$status, ")\n", sep = "")
  cat("  n:", length(x$w0), "  p:", ncol(x$X), "\n")
  if (!is.null(x$constraints)) {
    cat("  max |constraint gap|:", format(max(abs(x$constraints$gap)), digits = 3), "\n")
  }
  invisible(x)
}

#' @export
summary.calibration_fit <- function(object, ...) {
  x <- object
  out <- list(
    call = x$call,
    method = x$method,
    solver = x$solver,
    diagnostics = x$diagnostics,
    constraints = x$constraints,
    estimate = x$estimate,
    cov = x$cov
  )
  class(out) <- "summary.calibration_fit"
  out
}

#' @export
weights.calibration_fit <- function(object, ...) {
  object$w
}

#' @export
coef.calibration_fit <- function(object, ...) {
  object$lambda
}

#' @export
vcov.calibration_fit <- function(object, ...) {
  object$cov
}

#' @export
print.summary.calibration_fit <- function(x, ...) {
  cat("Calibration summary\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
  }
  cat("\n")
  cat("method:", x$method, "\n")
  cat("solver:", x$solver, "\n")
  if (!is.null(x$diagnostics) && !is.null(x$diagnostics$solver_status) && !is.na(x$diagnostics$solver_status)) {
    cat("solver_status:", x$diagnostics$solver_status, "\n")
  }
  if (!is.null(x$diagnostics)) {
    cat("converged:", isTRUE(x$diagnostics$converged), " (", x$diagnostics$status, ")\n", sep = "")
    if (!is.null(x$diagnostics$max_abs_constraint_gap)) {
      cat("max |constraint gap|:", format(x$diagnostics$max_abs_constraint_gap, digits = 4), "\n")
    }
    if (!is.null(x$diagnostics$message) && !is.na(x$diagnostics$message)) {
      cat("message:", x$diagnostics$message, "\n")
    }
  }
  if (!is.null(x$constraints)) {
    cat("\nConstraint gaps (largest 10 by |gap|):\n")
    tab <- x$constraints
    ord <- order(abs(tab$gap), decreasing = TRUE)
    tab <- tab[ord, , drop = FALSE]
    print(utils::head(tab, 10), row.names = FALSE)
  }
  if (!is.null(x$estimate)) {
    cat("\nEstimates:\n")
    print(x$estimate)
  }
  invisible(x)
}
