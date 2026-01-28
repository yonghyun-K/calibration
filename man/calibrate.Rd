\name{calibrate}
\alias{calibrate}
\title{Calibrate survey weights}
\description{
\code{calibrate} computes calibration weights under linear constraints and can optionally
compute totals and standard errors when a response is supplied on the left-hand side
of the formula.
}
\usage{
calibrate(
  formula,
  w0 = 1,
  data = NULL,
  const,
  method = c("BD", "GE", "DS"),
  entropy = c("EL", "ET", "HD", "CE", "PH", "SL"),
  divergence = NULL,
  w.scale = 1,
  G.scale = 1,
  bounds = NULL,
  solver = c("auto", "newton", "nleqslv", "cvxr"),
  solver_control = list(),
  ...
)
}
\arguments{
\item{formula}{A formula. The right-hand side defines calibration constraints.
If a response is provided, totals and standard errors are computed in one call.
Multivariate outcomes are supported via \code{cbind(y1, y2)}.
Factor/character outcomes are expanded to indicator columns (one per level), so the
returned totals correspond to calibrated counts by level.}
\item{w0}{Design weights (vector) or a scalar. If \code{w0} is a name and \code{data} is provided, it is looked up in \code{data}.}
\item{data}{Optional data frame.}
\item{const}{Numeric vector of known population totals. A named vector is recommended; names are matched to \code{colnames(model.matrix(...))}. \code{const} must not contain \code{NA}.}
\item{method}{Calibration method: \code{"BD"} (Bregman divergence), \code{"GE"} (generalized entropy), or \code{"DS"} (Deville--Sarndal divergence).}
\item{entropy}{Entropy family code (\code{"SL"}, \code{"EL"}, \code{"ET"}, \code{"HD"}, \code{"CE"}, \code{"PH"}), a numeric Renyi order, or a list with \code{code}/\code{family} and (optionally) \code{del} (e.g., \code{list(code = "PH", del = 0.5)}).}
\item{divergence}{Optional custom divergence specification. If supplied, this overrides \code{entropy}. Requirements depend on the solver: for \code{solver = "cvxr"}, supply \code{G(x)}; for \code{solver = "nleqslv"}, supply \code{g_inv(u, intercept = 0)}; for \code{solver = "newton"} or \code{"auto"}, supply \code{G(x)} and \code{g_inv(u, intercept = 0)}. Optional functions \code{g}, \code{g_prime_inv}, and \code{fprime} improve accuracy and avoid numeric fallbacks.}
\item{w.scale}{Weight scaling factor (phi).}
\item{G.scale}{Entropy scaling factor (q).}
\item{bounds}{Optional weight bounds. Can be a length-2 numeric vector \code{c(lower, upper)}, a list \code{list(lower=, upper=)}, or an \code{n x 2} matrix.

When \code{solver = "auto"} (the default) or \code{solver = "newton"}, bounds are handled via the dual clipped inverse-gradient map (piecewise), avoiding explicit inequality multipliers.

Bounds are also supported in \code{solver = "cvxr"} when the optional \pkg{CVXR} dependency is installed.}
\item{solver}{Solver backend. \code{"auto"} chooses the internal \code{"newton"} solver (no external dependencies).}
\item{solver_control}{List of backend controls (e.g. \code{maxit}, \code{xtol}, \code{ridge}).}
\item{...}{Reserved for future backend-specific options.}
}
\value{
A \code{calibration_fit} object.
}
\seealso{\code{\link{estimate}}, \code{\link{weights}}}
