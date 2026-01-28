# Internal entropy / divergence helpers

#' Normalize an entropy specification
#'
#' @param entropy Either a numeric value (Renyi order), a string code,
#'   or a list with `code`/`family` and (optionally) `del`.
#'   Supported codes: "SL", "EL", "ET", "HD", "CE", "PH".
#' @param del Optional threshold used by "PH". Defaults to 1 when omitted.
#'
#' @keywords internal
entropy_spec <- function(entropy, del = NULL) {
  if (is.list(entropy)) {
    code <- entropy$family
    if (is.null(code)) code <- entropy$code
    if (is.null(code)) code <- entropy$entropy
    if (is.null(code)) {
      stop("When entropy is a list, it must include 'family' or 'code'.")
    }
    if (!is.null(entropy$del)) del <- entropy$del
    entropy <- code
  }
  if (is.numeric(entropy)) {
    return(list(family = "renyi", r = as.numeric(entropy), del = del))
  }

  code <- toupper(as.character(entropy))
  if (!(code %in% c("SL", "EL", "ET", "HD", "CE", "PH"))) {
    stop("Invalid entropy value. Supported: SL, EL, ET, HD, CE, PH, or a numeric Renyi order.")
  }

  r <- switch(
    code,
    SL = 1,
    EL = -1,
    ET = 0,
    HD = -1 / 2,
    CE = NA_real_,
    PH = NA_real_
  )

  if (code == "PH") {
    if (is.null(del)) del <- 1
    if (length(del) != 1 || !is.finite(del) || del <= 0) {
      stop("When entropy = 'PH', 'del' must be a single positive number.")
    }
  }

  list(family = code, r = r, del = del)
}


#' Generator G(x) for built-in families
#'
#' This is used internally to evaluate dual objectives (via convex conjugates)
#' when bounds are handled through the truncated-conjugate / Huberization idea.
#'
#' The generator is only defined up to an additive linear term; Bregman
#' divergence is invariant to such shifts, so we pick convenient closed forms.
#'
#' @keywords internal
.G_impl <- function(x, spec) {
  if (!is.list(spec) || is.null(spec$family)) stop("Internal error: invalid entropy spec")

  if (identical(spec$family, "custom")) {
    # Optional: user-supplied generator (not required for unbounded problems)
    if (is.function(spec$G)) return(spec$G(x))
    return(rep(NA_real_, length(x)))
  }

  if (identical(spec$family, "CE")) {
    if (any(x <= 1, na.rm = TRUE)) return(rep(Inf, length(x)))
    return((x - 1) * log(x - 1) - x * log(x))
  }

  if (identical(spec$family, "PH")) {
    del <- spec$del
    if (is.null(del) || !is.finite(del) || del <= 0) return(rep(Inf, length(x)))
    return(del^2 * (sqrt(1 + (x / del)^2) - 1))
  }

  r <- spec$r
  if (is.null(r) || !is.finite(r)) {
    return(rep(NA_real_, length(x)))
  }

  # r == 0: exponential tilting / KL
  if (isTRUE(all.equal(r, 0))) {
    if (any(x <= 0, na.rm = TRUE)) return(rep(Inf, length(x)))
    return(x * log(x) - x)
  }

  # r == -1: empirical likelihood
  if (isTRUE(all.equal(r, -1))) {
    if (any(x <= 0, na.rm = TRUE)) return(rep(Inf, length(x)))
    return(-log(x))
  }

  # General Renyi-style family
  if (any(x < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
    return(rep(Inf, length(x)))
  }

  denom <- r * (r + 1)
  if (!is.finite(denom) || abs(denom) < .Machine$double.eps) {
    return(rep(Inf, length(x)))
  }
  x^(r + 1) / denom
}


#' Validate a user-supplied divergence specification
#'
#' A custom divergence must provide the following functions:
#' * g(x)
#' * g_inv(u, intercept = 0)
#' * g_prime_inv(u, intercept = 0)
#' * fprime(x)
#'
#' @keywords internal
.validate_custom_divergence <- function(divergence) {
  if (!is.list(divergence)) stop("'divergence' must be a list.")
  req <- c("g", "g_inv", "g_prime_inv", "fprime")
  missing <- req[!vapply(req, function(n) is.function(divergence[[n]]), logical(1))]
  if (length(missing) > 0) {
    stop("Custom divergence is missing required functions: ", paste(missing, collapse = ", "))
  }
  divergence$family <- "custom"
  divergence
}

# Helpers for numeric r: allow negative x only for odd positive integers.
.is_odd_positive_integer <- function(r) {
  if (!is.finite(r)) return(FALSE)
  rr <- round(r)
  isTRUE(all.equal(r, rr)) && rr > 0 && (rr %% 2 == 1)
}

#' Entropy derivative g(x) for built-in families
#'
#' @keywords internal
.g_impl <- function(x, spec) {
  if (!is.list(spec) || is.null(spec$family)) stop("Internal error: invalid entropy spec")

  if (identical(spec$family, "custom")) {
    return(spec$g(x))
  }

  if (identical(spec$family, "renyi")) {
    r <- spec$r
    if (r == 0) {
      if (any(x <= 0, na.rm = TRUE)) return(rep(Inf, length(x)))
      return(log(x))
    }
    if (any(x < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
      return(rep(Inf, length(x)))
    }
    return(x^r / r)
  }

  if (identical(spec$family, "CE")) {
    if (any(x <= 1, na.rm = TRUE)) return(rep(Inf, length(x)))
    return(log(1 - 1 / x))
  }

  if (identical(spec$family, "PH")) {
    del <- spec$del
    return(x / sqrt(1 + (x / del)^2))
  }

  # String-coded Renyi special cases are mapped to numeric r
  # (SL, EL, ET, HD)
  r <- spec$r
  if (r == 0) {
    if (any(x <= 0, na.rm = TRUE)) return(rep(Inf, length(x)))
    return(log(x))
  }
  if (any(x < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
    return(rep(Inf, length(x)))
  }
  x^r / r
}

#' Inverse derivative: g^{-1}(u) where g(w) = u
#'
#' This implementation supports an additive `intercept`, so that
#' g(w) = intercept + u.
#'
#' @keywords internal
.g_inv_impl <- function(u, spec, intercept = 0) {
  if (length(intercept) == 1) intercept <- rep(intercept, length(u))

  if (identical(spec$family, "custom")) {
    return(spec$g_inv(u, intercept = intercept))
  }

  if (identical(spec$family, "renyi")) {
    r <- spec$r
    if (r == 0) {
      return(exp(intercept + u))
    }
    base <- r * (intercept + u)
    if (any(base < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
      return(rep(Inf, length(u)))
    }
    return(base^(1 / r))
  }

  if (identical(spec$family, "CE")) {
    t <- intercept + u
    # Need exp(t) < 1 to keep weights positive and > 1
    if (any(t >= 0, na.rm = TRUE)) return(rep(Inf, length(u)))
    w <- 1 / (1 - exp(t))
    if (any(w <= 1, na.rm = TRUE)) return(rep(Inf, length(u)))
    return(w)
  }

  if (identical(spec$family, "PH")) {
    del <- spec$del
    t <- intercept + u
    if (any(abs(t) >= del, na.rm = TRUE)) return(rep(Inf, length(u)))
    # Equivalent to: x = t / sqrt(1 - (t/del)^2)
    return(1 / sqrt(1 / t^2 - 1 / del^2))
  }

  # String-coded Renyi special cases
  r <- spec$r
  if (r == 0) {
    return(exp(intercept + u))
  }
  base <- r * (intercept + u)
  if (any(base < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
    return(rep(Inf, length(u)))
  }
  base^(1 / r)
}

#' Derivative of g^{-1} with respect to its input u
#'
#' @keywords internal
.g_prime_inv_impl <- function(u, spec, intercept = 0) {
  if (length(intercept) == 1) intercept <- rep(intercept, length(u))

  if (identical(spec$family, "custom")) {
    return(spec$g_prime_inv(u, intercept = intercept))
  }

  if (identical(spec$family, "renyi")) {
    r <- spec$r
    if (r == 0) {
      return(exp(intercept + u))
    }
    base <- r * (intercept + u)
    if (any(base < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
      return(rep(Inf, length(u)))
    }
    return(base^(1 / r - 1))
  }

  if (identical(spec$family, "CE")) {
    t <- intercept + u
    if (any(t >= 0, na.rm = TRUE)) return(rep(Inf, length(u)))
    return(exp(t) / (1 - exp(t))^2)
  }

  if (identical(spec$family, "PH")) {
    del <- spec$del
    t <- intercept + u
    if (any(abs(t) >= del, na.rm = TRUE)) return(rep(Inf, length(u)))
    return((1 - (t / del)^2)^(-1.5))
  }

  # String-coded Renyi special cases
  r <- spec$r
  if (r == 0) {
    return(exp(intercept + u))
  }
  base <- r * (intercept + u)
  if (any(base < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
    return(rep(Inf, length(u)))
  }
  base^(1 / r - 1)
}

#' 1 / g'(x) (used in variance estimation)
#'
#' @keywords internal
.fprime_impl <- function(x, spec) {
  if (identical(spec$family, "custom")) {
    return(spec$fprime(x))
  }

  if (identical(spec$family, "renyi")) {
    r <- spec$r
    if (r == 0) {
      if (any(x <= 0, na.rm = TRUE)) return(rep(Inf, length(x)))
      return(x)
    }
    if (any(x < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
      return(rep(Inf, length(x)))
    }
    return(x^(1 - r))
  }

  if (identical(spec$family, "CE")) {
    return(x * (x - 1))
  }

  if (identical(spec$family, "PH")) {
    del <- spec$del
    return((1 + (x / del)^2)^(1.5))
  }

  # String-coded Renyi
  r <- spec$r
  if (r == 0) {
    if (any(x <= 0, na.rm = TRUE)) return(rep(Inf, length(x)))
    return(x)
  }
  if (any(x < 0, na.rm = TRUE) && !.is_odd_positive_integer(r)) {
    return(rep(Inf, length(x)))
  }
  x^(1 - r)
}

#' Exported helper: entropy derivative g(x)
#'
#' This is a simple, explicit wrapper around the internal implementation.
#' It does **not** look for values in the caller environment.
#'
#' @param x Numeric vector.
#' @param entropy Entropy family ("SL", "EL", "ET", "HD", "CE", "PH")
#'   or a numeric Renyi order.
#' @param del Optional threshold used when entropy = "PH". Defaults to 1 when omitted.
#'
#' @return Numeric vector.
#' @export
calibration_g <- function(x, entropy = c("SL", "EL", "ET", "HD", "CE", "PH"), del = NULL) {
  spec <- entropy_spec(entropy, del = del)
  .g_impl(x, spec)
}
