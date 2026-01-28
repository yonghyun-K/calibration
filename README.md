# calibration

`calibration` provides an integrated workflow for survey-weight calibration with a modern API,
structured solver diagnostics, and extensible solver backends.

It is a refactoring / redesign inspired by the existing **GECal** package.

## Installation (from source)

```r
install.packages("nleqslv")
# Optional (extra robustness)
# install.packages("CVXR")

install.packages("calibration_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick start

### 1) Calibrate weights only

```r
library(calibration)

pop_totals <- c(`(Intercept)` = 6194, stypeH = 755, stypeM = 1018)
fit <- calibrate(~ stype, w0 = 1, data = apiclus1,
                 const = pop_totals,
                 method = "BD",
                 entropy = "SL")

w <- weights(fit)
summary(fit)
```

### 2) One-call workflow: weights + estimated totals and SEs

```r
fit2 <- calibrate(api00 ~ stype, w0 = 1, data = apiclus1,
                  const = pop_totals,
                  method = "BD",
                  entropy = "SL")

fit2$estimate
fit2$cov
```

## Notes

* `const` **must not** contain any `NA` values.
* For safety, `const` should be a **named** numeric vector; names are matched to
  `colnames(model.matrix(...))`.
* Solver diagnostics are stored in `fit$diagnostics`.
* Bounds (`bounds = c(lower, upper)`) are supported in the default `solver = "auto"`
  via a dual clipped inverse-gradient map (no CVXR required).
