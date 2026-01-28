\name{estimate}
\alias{estimate}
\title{Estimate totals and standard errors after calibration}
\description{
Computes calibrated-weight totals and a linearized covariance estimator.
}
\usage{
estimate(formula, data = NULL, calibration, pimat = NULL)
}
\arguments{
\item{formula}{A formula with outcomes on the left-hand side. The right-hand side is ignored.
Multivariate outcomes can be supplied via \code{cbind(y1, y2)}.
If the outcome is a factor/character vector, it is expanded to indicator columns (one per level)
and totals/standard errors are returned for each level.}
\item{data}{Optional data frame.}
\item{calibration}{A \code{calibration_fit} object returned by \code{calibrate}.}
\item{pimat}{Optional joint inclusion probability structure for variance estimation. 
This can be \code{NULL} (a default diagonal approximation), a numeric vector of length \code{n} (interpreted as the diagonal of \code{pimat}), or an \code{n x n} matrix.}
}
\value{
A list with elements \code{estimate} (point estimates and standard errors) and \code{cov} (covariance matrix).
}
\seealso{\code{\link{calibrate}}}
