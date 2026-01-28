# Internal helpers for constraint handling and diagnostics

#' @keywords internal
.match_const <- function(const, X) {
  if (is.null(const)) stop("'const' must be provided.")
  if (!is.numeric(const)) stop("'const' must be numeric.")
  if (anyNA(const)) stop("'const' must not contain NA values.")

  p <- ncol(X)
  x_names <- colnames(X)
  if (is.null(x_names) || any(x_names == "")) {
    # model.matrix should usually provide names; if it didn't, fall back to index names
    x_names <- paste0("V", seq_len(p))
    colnames(X) <- x_names
  }

  cn <- names(const)
  is_named <- !is.null(cn) && any(nzchar(cn))

  if (!is_named) {
    if (length(const) != p) {
      stop("'const' must have length ncol(model.matrix(...)) when not named.")
    }
    names(const) <- x_names
    return(const)
  }

  # Enforce exact matching by name for safety
  missing <- setdiff(x_names, cn)
  extra <- setdiff(cn, x_names)
  if (length(missing) > 0) {
    stop(
      "'const' is missing required names: ",
      paste(missing, collapse = ", ")
    )
  }
  if (length(extra) > 0) {
    stop(
      "'const' has extra names not in model matrix: ",
      paste(extra, collapse = ", ")
    )
  }

  const[x_names]
}

#' @keywords internal
.constraint_table <- function(X, w, const) {
  achieved <- colSums(X * w)
  data.frame(
    term = colnames(X),
    target = as.numeric(const),
    achieved = as.numeric(achieved),
    gap = as.numeric(achieved - const),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.rank_check <- function(X, tol = sqrt(.Machine$double.eps)) {
  qrX <- qr(X, tol = tol)
  list(rank = qrX$rank, p = ncol(X), full_rank = (qrX$rank == ncol(X)))
}

#' Very lightweight feasibility diagnostics
#'
#' These are *diagnostics*, not proofs.
#'
#' @keywords internal
.diagnose_feasibility <- function(X, const, require_positive = TRUE) {
  out <- list(
    likely_infeasible = FALSE,
    reasons = character(0),
    details = list()
  )

  rk <- .rank_check(X)
  out$details$rank <- rk
  if (!rk$full_rank) {
    out$likely_infeasible <- TRUE
    out$reasons <- c(out$reasons, sprintf("Design matrix is rank-deficient (rank=%s < p=%s).", rk$rank, rk$p))
  }

  # Simple marginal convex-hull check when an intercept is present and weights are constrained to be positive.
  if (require_positive && "(Intercept)" %in% colnames(X) && "(Intercept)" %in% names(const)) {
    N <- as.numeric(const["(Intercept)"])
    if (is.finite(N) && N > 0) {
      idx <- which(colnames(X) != "(Intercept)")
      if (length(idx) > 0) {
        mins <- apply(X[, idx, drop = FALSE], 2, min)
        maxs <- apply(X[, idx, drop = FALSE], 2, max)
        targets_mean <- const[colnames(X)[idx]] / N
        # compare to sample range of each column
        outside <- which(targets_mean < mins | targets_mean > maxs)
        if (length(outside) > 0) {
          out$likely_infeasible <- TRUE
          out$reasons <- c(out$reasons,
                           "Some target means lie outside the sample min/max range (necessary condition for positive weights when (Intercept) constraint fixes sum(w)).")
          out$details$marginal_outside <- data.frame(
            term = colnames(X)[idx][outside],
            target_mean = as.numeric(targets_mean[outside]),
            min = as.numeric(mins[outside]),
            max = as.numeric(maxs[outside]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  out
}
