.onAttach <- function(libname, pkgname) {
  # Helpful interactive note: avoid confusion with survey::calibrate()
  packageStartupMessage(
    "Loaded 'calibration' (Bregman/entropy weight calibration) v",
    utils::packageVersion(pkgname),
    ".\n",
    "Note: 'survey' also defines calibrate(). If you see method-dispatch errors, use ",
    "calibration::calibrate() or calibrate_weights()."
  )
}

# Register an S3 method for the survey::calibrate() generic when that package is attached.
# This lets users call calibrate(~x, ...) even if survey::calibrate masks calibration::calibrate.
#
# We *avoid* requireNamespace("survey") here to prevent side-effects (loading a heavy namespace).
.onLoad <- function(libname, pkgname) {
  register_for_survey <- function(...) {
    if (!requireNamespace("survey", quietly = TRUE)) {
      return(invisible(FALSE))
    }

    ok <- FALSE

    # 1) Standard S3 registration (works on most R versions)
    ok <- ok || tryCatch({
      utils::registerS3method("calibrate", "formula", calibrate.formula)
      TRUE
    }, error = function(e) FALSE)

    # 2) Some setups are picky about where the generic lives.
    ok <- ok || tryCatch({
      utils::registerS3method("calibrate", "formula", calibrate.formula, envir = asNamespace("survey"))
      TRUE
    }, error = function(e) FALSE)

    # 3) Last-resort: place the method directly into the survey namespace so that
    #    UseMethod() can find it even without registration.
    ok <- ok || tryCatch({
      utils::assignInNamespace("calibrate.formula", calibrate.formula, ns = "survey")
      TRUE
    }, error = function(e) FALSE)

    invisible(ok)
  }

  # If survey is already loaded, register immediately.
  if ("survey" %in% loadedNamespaces()) {
    try(register_for_survey(), silent = TRUE)
  }

  # Register when survey is loaded or attached later.
  try(setHook(packageEvent("survey", "onLoad"), register_for_survey, action = "append"), silent = TRUE)
  try(setHook(packageEvent("survey", "onAttach"), register_for_survey, action = "append"), silent = TRUE)
}
