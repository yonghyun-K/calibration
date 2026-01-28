# Estimation utilities (totals + variance) after calibration

.coerce_outcome_matrix <- function(y, response_name = "y") {
  # Convert various response types into a numeric matrix.
  #
  # Supported:
  # * numeric/integer -> 1-column matrix
  # * logical -> 1-column 0/1 matrix
  # * factor/character -> indicator matrix (one column per level)
  # * numeric matrix -> returned as-is
  #
  # This is used by estimate() so that factor outcomes can be handled
  # (counts/weighted totals by level).

  if (is.matrix(y)) {
    if (is.numeric(y)) return(y)
    if (is.logical(y)) return(matrix(as.numeric(y), nrow = nrow(y), dimnames = dimnames(y)))
    stop("Outcome matrix is not numeric. Use numeric outcomes or `cbind()` numeric outcomes.")
  }

  if (is.data.frame(y)) {
    # Common case: model.response can return a data.frame in some multivariate setups.
    # We only allow fully numeric/logical columns here.
    ok <- vapply(y, function(col) is.numeric(col) || is.integer(col) || is.logical(col), logical(1))
    if (!all(ok)) {
      stop("Outcome data.frame contains non-numeric columns. ",
           "Use `cbind()` with numeric outcomes, or convert factors to numeric/indicators explicitly.")
    }
    ymat <- as.matrix(y)
    storage.mode(ymat) <- "numeric"
    return(ymat)
  }

  if (is.factor(y) || is.character(y) || is.ordered(y)) {
    f <- if (is.factor(y) || is.ordered(y)) y else factor(y)
    mm <- stats::model.matrix(~ f - 1)
    # Nice, stable column names: <response>:<level>
    lev <- levels(f)
    if (length(lev) == ncol(mm)) {
      colnames(mm) <- paste0(response_name, ":", lev)
    }
    storage.mode(mm) <- "numeric"
    return(mm)
  }

  if (is.logical(y)) {
    mm <- matrix(as.numeric(y), ncol = 1)
    colnames(mm) <- response_name
    return(mm)
  }

  if (is.numeric(y) || is.integer(y)) {
    mm <- matrix(as.numeric(y), ncol = 1)
    colnames(mm) <- response_name
    return(mm)
  }

  stop("Unsupported outcome type: ", paste(class(y), collapse = "/"), ".")
}

#' Estimate totals and standard errors after calibration
#'
#' @param formula A formula with outcomes on the left-hand side.
#'   Use `y ~ 1` for a single outcome, or `cbind(y1, y2) ~ 1` for multiple.
#' @param data Optional data frame.
#' @param calibration A `calibration_fit` object returned by [calibrate()].
#' @param pimat Optional joint inclusion probability structure used for variance.
#'   This can be:
#'   * `NULL` (use a default diagonal approximation),
#'   * a numeric vector of length `n` (interpreted as the diagonal of `pimat`),
#'   * a matrix (dense or sparse) with dimension `n x n`.
#'
#' @return A list with elements `estimate` (table) and `cov` (covariance matrix).
#' @export
estimate <- function(formula, data = NULL, calibration, pimat = NULL) {
  if (!inherits(calibration, "calibration_fit")) {
    stop("'calibration' must be a calibration_fit object returned by calibrate().")
  }

  # Extract outcome(s)
  if (is.null(data)) {
    mf <- stats::model.frame(formula, na.action = stats::na.fail)
  } else {
    mf <- stats::model.frame(formula, data = data, na.action = stats::na.fail)
  }
  ys_raw <- stats::model.response(mf)
  if (is.null(ys_raw)) stop("No response found on the left-hand side of 'formula'.")

  # Try to infer a stable response name for labeling.
  resp_idx <- attr(stats::terms(mf), "response")
  resp_name <- if (!is.null(resp_idx) && length(resp_idx) == 1 && resp_idx > 0) {
    names(mf)[resp_idx]
  } else {
    "y"
  }

  ys <- .coerce_outcome_matrix(ys_raw, response_name = resp_name)

  X <- calibration$X
  w <- calibration$w
  w0 <- calibration$w0
  method <- calibration$method
  spec <- calibration$entropy_spec
  w_scale <- calibration$w_scale
  G_scale <- calibration$G_scale

  if (length(w) != nrow(ys)) {
    stop("Length of calibrated weights does not match number of outcome rows.")
  }

  # Default diagonal approximation (avoids dense n x n matrices)
  #
  # NOTE:
  # * A Poisson-style diagonal uses (1 - pi_i) / pi_i^2.
  # * If weights are on a *total* scale (sum(w) ~= N), a common proxy is w*(w-1).
  # * If weights are on a *mean* scale (sum(w) ~= 1), w*(w-1) can be negative.
  #   In that case, fall back to design-weight-based scaling when possible.
  if (is.null(pimat)) {
    sum_w <- sum(w)
    intercept_target <- NA_real_
    if (!is.null(calibration$const) && "(Intercept)" %in% names(calibration$const)) {
      intercept_target <- as.numeric(calibration$const["(Intercept)"])
    }

    # Heuristic: treat as "mean-scale" if the intercept constraint (and achieved sum of weights)
    # is approximately 1.
    is_mean_scale <- is.finite(sum_w) && abs(sum_w - 1) <= 1e-6
    if (is.finite(intercept_target)) {
      is_mean_scale <- is_mean_scale && abs(intercept_target - 1) <= 1e-8
    }

    if (is_mean_scale) {
      # For mean-scale estimators, approximate Var(mean) = Var(total) / N^2,
      # with N estimated by sum(w0) when N is not explicitly supplied.
      N_hat <- sum(w0)
      if (!is.finite(N_hat) || N_hat <= 0) N_hat <- nrow(ys)

      # If w0 varies (typical design weights), use it to form the diagonal proxy.
      # Otherwise, fall back to scaling calibrated weights back to a "total" scale.
      w0_var <- stats::var(w0)
      w0_nonconst <- is.finite(w0_var) && w0_var > 0
      base_w <- if (w0_nonconst || any(w0 > 1 + 1e-8, na.rm = TRUE)) {
        w0
      } else {
        w * N_hat
      }

      pimat <- base_w * (base_w - 1) / (N_hat^2)
    } else {
      # Total-scale default: use calibrated weights as a proxy.
      pimat <- w * (w - 1)
    }

    if (any(!is.finite(pimat))) {
      warning(
        "Default variance approximation produced non-finite values; consider providing 'pimat'.",
        call. = FALSE
      )
    }

    if (any(pimat < 0, na.rm = TRUE)) {
      warning(
        "Default variance diagonal had negative entries (common when weights are on a mean scale). ",
        "Negative values were truncated to zero. Consider providing 'pimat' for design-based variance.",
        call. = FALSE
      )
      pimat <- pmax(pimat, 0)
    }
  }

  # Regression-type adjustment (GREG-style)
  if (ncol(X) == 0) {
    yhat <- matrix(0, nrow = nrow(ys), ncol = ncol(ys))
  } else {
    if (identical(method, "DS")) {
      v <- w0 / G_scale
      A <- crossprod(X, X * v)
      B <- crossprod(X, ys * v)
      gammahat <- solve(A, B)
    } else {
      fp <- .fprime_impl(w0 * w_scale, spec = spec) / (G_scale * w_scale^2)
      A <- crossprod(X, X * fp)
      B <- crossprod(X, ys * fp)
      gammahat <- solve(A, B)
    }
    yhat <- X %*% gammahat
  }

  resid <- ys - yhat

  # Variance: t(resid) %*% pimat %*% resid
  if (is.numeric(pimat) && length(pimat) == nrow(resid)) {
    # diagonal represented as a vector
    Varhat <- crossprod(resid, resid * pimat)
  } else {
    Varhat <- crossprod(resid, pimat %*% resid)
  }

  est <- colSums(ys * w)
  # Guard against small negative diagonals due to numerical error
  dV <- diag(Varhat)
  if (any(dV < -sqrt(.Machine$double.eps), na.rm = TRUE)) {
    warning(
      "Negative variance estimates encountered; truncating negatives to zero. ",
      "Provide 'pimat' for a more appropriate design-based variance if needed.",
      call. = FALSE
    )
  }
  se <- sqrt(pmax(dV, 0))

  out <- list(
    cov = Varhat,
    estimate = cbind(Estimate = est, `Std. Error` = se)
  )
  class(out) <- "calibration_estimate"
  out
}
