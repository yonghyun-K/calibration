# Optional CVXR backend (robust, supports bounds)

#' @keywords internal
.parse_bounds <- function(bounds, n) {
  if (is.null(bounds)) return(list(lower = NULL, upper = NULL))

  # Treat c(-Inf, Inf) (or all-infinite bounds) as "no bounds"
  if (is.numeric(bounds) && length(bounds) == 2 && all(!is.finite(bounds))) {
    return(list(lower = NULL, upper = NULL))
  }

  lower <- NULL
  upper <- NULL

  if (is.numeric(bounds) && length(bounds) == 2) {
    lower <- rep(bounds[1], n)
    upper <- rep(bounds[2], n)
  } else if (is.list(bounds)) {
    if (!is.null(bounds$lower)) lower <- bounds$lower
    if (!is.null(bounds$upper)) upper <- bounds$upper
  } else if (is.matrix(bounds) && ncol(bounds) == 2) {
    lower <- bounds[, 1]
    upper <- bounds[, 2]
  } else {
    stop("'bounds' must be NULL, a length-2 numeric vector, a list(lower=, upper=), or an n x 2 matrix.")
  }

  if (!is.null(lower)) {
    if (length(lower) == 1) lower <- rep(lower, n)
    if (length(lower) != n) stop("bounds$lower must be length 1 or n.")
  }
  if (!is.null(upper)) {
    if (length(upper) == 1) upper <- rep(upper, n)
    if (length(upper) != n) stop("bounds$upper must be length 1 or n.")
  }

  if (!is.null(lower) && !is.null(upper)) {
    if (any(lower > upper)) stop("bounds$lower must be <= bounds$upper elementwise.")
  }

  list(lower = lower, upper = upper)
}

#' @keywords internal
.solve_cvxr <- function(X, const, w0, method, spec,
                        w_scale = 1, G_scale = 1,
                        bounds = NULL,
                        solver_control = list()) {

  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("CVXR backend requested, but package 'CVXR' is not installed. Install it or use solver='nleqslv'.")
  }

  n <- nrow(X)
  p <- ncol(X)

  if (!(method %in% c("BD", "DS"))) {
    stop("CVXR backend currently supports method = 'BD' or 'DS' only.")
  }

  # --- Determine which built-in entropy family we can represent with CVXR ---
  # We currently support the common GECal families that admit standard DCP-safe
  # expressions: SL, ET, EL, and HD.
  fam <- spec$family
  r <- spec$r
  code <- NA_character_

  if (identical(fam, "custom")) {
    code <- "custom"
  } else if (identical(fam, "renyi")) {
    if (isTRUE(all.equal(r, 1))) code <- "SL"
    else if (isTRUE(all.equal(r, 0))) code <- "ET"
    else if (isTRUE(all.equal(r, -1))) code <- "EL"
    else if (isTRUE(all.equal(r, -1/2))) code <- "HD"
  } else if (fam %in% c("SL", "ET", "EL", "HD")) {
    code <- fam
  }

  if (!isTRUE(nzchar(code))) {
    stop(
      "CVXR backend currently supports entropy in { 'SL', 'ET', 'EL', 'HD' } only.\n",
      "(You requested: ", fam, if (!is.null(r) && is.finite(r)) paste0(", r=", r) else "", ")"
    )
  }

  bnd <- .parse_bounds(bounds, n)

  w <- CVXR::Variable(n)

constraints <- list(t(X) %*% w == const)

if (!is.null(bnd$lower)) {
  idx <- is.finite(bnd$lower)
  if (any(idx)) constraints <- c(constraints, list(w[idx] >= bnd$lower[idx]))
}
if (!is.null(bnd$upper)) {
  idx <- is.finite(bnd$upper)
  if (any(idx)) constraints <- c(constraints, list(w[idx] <= bnd$upper[idx]))
}

  # Domain constraints for positive-domain divergences
  domain_eps <- if (!is.null(solver_control$domain_eps)) solver_control$domain_eps else 1e-8
  if (code %in% c("ET", "EL", "HD")) {
    if (!is.finite(w_scale) || w_scale <= 0) {
      stop("For entropy in {ET, EL, HD}, 'w.scale' must be a positive number when using solver='cvxr'.")
    }
    # Enforce strict positivity (numerically) unless bounds already do it.
    constraints <- c(constraints, list(w >= domain_eps))
  }
  if (identical(code, "custom")) {
    if (!is.function(spec$G)) {
      stop("Custom divergence for solver='cvxr' requires a 'G' function.")
    }
    if (!identical(method, "BD")) {
      stop("Custom divergence with solver='cvxr' currently supports method = 'BD' only.")
    }
    g_neg <- try(spec$G(-1), silent = TRUE)
    g_pos <- try(spec$G(1), silent = TRUE)
    if (!inherits(g_neg, "try-error") && !inherits(g_pos, "try-error")) {
      if (!is.finite(g_neg) && is.finite(g_pos)) {
        constraints <- c(constraints, list(w >= domain_eps))
      } else if (is.finite(g_neg) && !is.finite(g_pos)) {
        constraints <- c(constraints, list(w <= -domain_eps))
      }
    }
  }

  # Choose a CVXR solver if the user didn't specify one.
  # * For pure QP (SL) we prefer OSQP when available.
  # * For exponential-cone problems (ET/EL) we prefer ECOS, then SCS.
  cvxr_solver <- if (!is.null(solver_control$cvxr_solver)) solver_control$cvxr_solver else NULL
  if (is.null(cvxr_solver)) {
    if (identical(code, "SL") && requireNamespace("osqp", quietly = TRUE)) {
      cvxr_solver <- "OSQP"
    } else if (requireNamespace("ECOS", quietly = TRUE)) {
      cvxr_solver <- "ECOS"
    } else if (requireNamespace("scs", quietly = TRUE)) {
      cvxr_solver <- "SCS"
    } else {
      cvxr_solver <- NULL
    }
  }

  # Build the convex objective.
  objective_expr <- NULL

  if (identical(method, "DS")) {
    if (any(!is.finite(w0)) || any(w0 <= 0)) {
      stop("For method='DS' with solver='cvxr', w0 must be positive and finite.")
    }
    if (identical(code, "SL")) {
      # Chi-square distance on the ratio w/w0
      objective_expr <- CVXR::square(w - w0) / (2 * w0)
    } else if (identical(code, "ET")) {
      # Relative entropy: sum kl_div(w, w0) = w log(w/w0) - w + w0
      objective_expr <- CVXR::kl_div(w, w0)
    } else if (identical(code, "EL")) {
      # Dual KL form: kl_div(w0, w) = w0 log(w0/w) - w0 + w
      objective_expr <- CVXR::kl_div(w0, w)
    } else if (identical(code, "HD")) {
      # r = -1/2 implies a Hellinger-type distance on the ratio.
      # Up to constants, objective is: 2*w - 4*sqrt(w*w0)  (convex for w>=0, w0>0)
      objective_expr <- 2 * w - 4 * sqrt(w * w0)
    }
  } else {
    # method == "BD"
    z <- w * w_scale
    z0 <- w0 * w_scale
    g_z0 <- .g_impl(z0, spec)
    if (any(!is.finite(g_z0))) {
      stop("Non-finite g(w0) encountered. Check that w0 lies in the divergence domain for the chosen entropy.")
    }

    if (identical(code, "custom")) {
      objective_expr <- tryCatch(spec$G(z), error = function(e) NULL)
      if (is.null(objective_expr)) {
        g_vals <- try(spec$G(c(0, 1, 2)), silent = TRUE)
        g_neg <- try(spec$G(-1), silent = TRUE)
        if (!inherits(g_vals, "try-error") &&
            isTRUE(all.equal(g_vals, c(0, 1, 4))) &&
            !inherits(g_neg, "try-error") && !is.finite(g_neg)) {
          objective_expr <- CVXR::square(z)
        } else {
          stop("Custom divergence G must return a CVXR-compatible expression.")
        }
      }
      objective_expr <- objective_expr - z * g_z0
    } else if (identical(code, "SL")) {
      # G(z)=0.5 z^2, g(z0)=z0
      objective_expr <- 0.5 * CVXR::square(z) - z * g_z0
    } else if (identical(code, "ET")) {
      # G(z)=z log z - z, g(z0)=log(z0)
      objective_expr <- (-CVXR::entr(z)) - z - z * g_z0
    } else if (identical(code, "EL")) {
      # G(z)=-log z, g(z0)=-1/z0
      objective_expr <- (-log(z)) - z * g_z0
    } else if (identical(code, "HD")) {
      # r=-1/2: G(z)=-4*sqrt(z), g(z0)=-2/sqrt(z0)
      objective_expr <- (-4 * sqrt(z)) - z * g_z0
    }
  }

  objective <- CVXR::Minimize(G_scale * CVXR::sum_entries(objective_expr))

  problem <- CVXR::Problem(objective, constraints)

  cvxr_opts <- if (!is.null(solver_control$cvxr_opts)) solver_control$cvxr_opts else list()

  sol <- try(
    do.call(CVXR::solve, c(list(problem, solver = cvxr_solver), cvxr_opts)),
    silent = TRUE
  )

  diagnostics <- list(
    converged = FALSE,
    solver = "cvxr",
    solver_status = NA_character_,
    message = NA_character_,
    status = NA_character_,
    max_abs_constraint_gap = Inf
  )

  if (inherits(sol, "try-error")) {
    diagnostics$message <- as.character(sol)
    return(list(w = rep(NA_real_, n), lambda = rep(NA_real_, p), diagnostics = diagnostics, solver_raw = NULL))
  }

  diagnostics$solver_status <- sol$status
  diagnostics$message <- sol$status

  w_hat <- as.numeric(sol$getValue(w))
  achieved <- colSums(X * w_hat)
  gap <- achieved - const
  diagnostics$max_abs_constraint_gap <- max(abs(gap))

  bounds_tol <- if (!is.null(solver_control$bounds_tol)) solver_control$bounds_tol else 1e-6
  bounds_ok <- TRUE
  bounds_msg <- NULL
  if (!is.null(bnd$lower)) {
    idx <- is.finite(bnd$lower)
    if (any(idx)) {
      lower_slack <- bnd$lower[idx] - w_hat[idx]
      max_lower_slack <- max(lower_slack, na.rm = TRUE)
      if (is.finite(max_lower_slack) && max_lower_slack > bounds_tol) {
        bounds_ok <- FALSE
        bounds_msg <- paste0(
          "weights violate lower bounds (max slack: ",
          format(max_lower_slack, digits = 6), ")"
        )
      }
    }
  }
  if (!is.null(bnd$upper)) {
    idx <- is.finite(bnd$upper)
    if (any(idx)) {
      upper_slack <- w_hat[idx] - bnd$upper[idx]
      max_upper_slack <- max(upper_slack, na.rm = TRUE)
      if (is.finite(max_upper_slack) && max_upper_slack > bounds_tol) {
        bounds_ok <- FALSE
        upper_msg <- paste0(
          "weights violate upper bounds (max slack: ",
          format(max_upper_slack, digits = 6), ")"
        )
        bounds_msg <- if (is.null(bounds_msg)) upper_msg else paste(bounds_msg, upper_msg, sep = "; ")
      }
    }
  }
  if (!bounds_ok) diagnostics$message <- bounds_msg

  # Decide convergence primarily based on CVXR's own termination status.
  # CVXR backends vary widely in how tightly they enforce equality constraints
  # (e.g., OSQP/SCS defaults can be much looser than 1e-8). If the solver reports
  # an optimal solution, we return weights and surface the achieved constraint gaps
  # in diagnostics/summary().
  status_txt <- as.character(sol$status)
  status_ok <- isTRUE(nzchar(status_txt)) && grepl("^optimal", tolower(status_txt))

  diagnostics$converged <- status_ok &&
    bounds_ok &&
    all(is.finite(w_hat)) &&
    is.finite(diagnostics$max_abs_constraint_gap)

list(w = if (diagnostics$converged) w_hat else rep(NA_real_, n),
     lambda = rep(NA_real_, p),
     diagnostics = diagnostics,
     solver_raw = sol)
}
