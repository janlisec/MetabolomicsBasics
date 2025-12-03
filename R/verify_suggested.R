#' @title verify_suggested.
#' @description Check if packages are available and stop function otherwise.
#' @param pkg Package names to be checked.
#' @param default A default return value might be passed to the function for the case of missing packages.
#' @return NULL.
#' @keywords internal
#' @noRd
verify_suggested <- function(pkg, default = NULL) {
  # verify that suggested packages are available
  check_pkg <- sapply(pkg, requireNamespace, quietly = TRUE)
  if (!all(check_pkg)) {
    warning(paste0(
      "The use of this function requires package", ifelse(sum(!check_pkg)>1, "s ", " "),
      paste(names(check_pkg)[!check_pkg], collapse=", "),
      ". Please install."
    ))
    invisible(default)
  }
  invisible(NULL)
}
