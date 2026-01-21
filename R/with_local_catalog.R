#' @title with_local_catalog.
#'
#' @details Will modify system parameter value "XML_CATALOG_FILES" (for xml2)
#'     to include a local catalog file for function `recalib_mzML()`.
#'
#' @param expr expr.
#' @param catalog catalog.
#'
#' @keywords internal
#' @noRd
with_local_catalog <- function(expr, catalog) {
  old <- Sys.getenv("XML_CATALOG_FILES", unset = "")
  cat_val <- if (nzchar(old)) paste(c(catalog, old), collapse = .Platform$path.sep) else catalog
  withr::local_envvar(c(XML_CATALOG_FILES = cat_val))
  force(expr)
}
