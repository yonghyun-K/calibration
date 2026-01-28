\name{calibration_g}
\alias{calibration_g}
\title{Entropy derivative g(x)}
\description{
Explicit helper for computing the derivative of the built-in entropy generator.
Unlike older designs, this function does not depend on variables in the caller environment.
}
\usage{
calibration_g(x, entropy = c("SL", "EL", "ET", "HD", "CE", "PH"), del = NULL)
}
\arguments{
\item{x}{Numeric vector.}
\item{entropy}{Entropy family code or numeric Renyi order.}
\item{del}{Optional threshold used when entropy = "PH".}
}
\value{
A numeric vector.
}
