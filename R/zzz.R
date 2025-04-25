.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "NOTE: The package '", pkgname, "' is DEPRECATED and will be archived on July 1st, 2025.\n",
    "Please migrate to the new package: 'joker'."
  )
}
