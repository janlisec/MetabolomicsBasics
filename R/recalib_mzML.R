#' @title recalib_mzML
#'
#' @description If calibration fails, mass drift may occur within a measurement
#'     batch, leading to m/z values in the raw data being offset with respect to
#'     the true value. This often can be fixed using software tools provided by
#'     the vendor. However, it may become necessary to re-calibrate mzML data
#'     after conversion from (incorrect) raw data. This function allows to shift
#'     m/z values in MS1 and/or MS2 scans of mzML files by a fixed value.
#'     Except for base_peak and precursor values assigned to spectra, which are
#'     updated, all other nodes of the mzML file will be preserved. The binary
#'     encoding (32-bit or 64-bit of the measurement data is respected by the
#'     function.)
#'
#' @param infile mzML file path.
#' @param outfile If no target file is specified, input file name will be used (extended by '_calib').
#' @param shift_da The desired m/z shift in Dalton.
#' @param apply_to Use keyword 'all' to modify MS1 and MS2 or specify selectively.
#'
#' @returns The modified XML invisibly (for potential checks).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   fld <- "C:/Users/jlisec/Documents/Projects/PFAS_nontargeted/TSo_GMW/mzML/original/"
#'   fls <- dir(fld, pattern = "mzML$", full.names = TRUE)
#'   fl <- fls[1]
#'   for (fl in fls) {
#'     print(fl)
#'     recalib_mzML(infile = fl, shift_da = -0.012)
#'   }
#' }
recalib_mzML <- function(
    infile,
    outfile = sub("\\.[mM][zZ][mM][lL]$", "_calib.mzML", infile, ignore.case = TRUE),
    shift_da = -0.019,
    apply_to = c("all", "MS1", "MS2")
) {

  apply_to <- match.arg(apply_to)

  # read validated mzML File
  # set schema catalog for xml2
  #Sys.setenv(XML_CATALOG_FILES = normalizePath(system.file("extdata", "mzml_xsd", "catalog.xml", package = "MetabolomicsBasics")))

  # helper function in validation/xml loading
  to_file_uri <- function(x) {
    if (is.na(x) || !nzchar(x)) return(x)

    # 1) Backslashes -> forward slashes
    y <- gsub("\\\\", "/", x, fixed = FALSE)

    # 2) Falls es mit "file://" beginnt, sicherstellen, dass es "file:///" ist
    #    (leere Authority für lokalen Pfad). Für UNC behandeln wir separat.
    if (grepl("^file://", y, ignore.case = TRUE)) {
      # UNC-Host? (file://server/...)
      # Wenn nach "file://" direkt ein Buchstabe und ":" kommt, ist es ein Laufwerk → fehlendes "/"
      y <- sub("^file://([A-Za-z]:)", "file:///$1", y)
    } else {
      # Reiner Windows-Pfad wie "C:/..." -> als lokale file-URI
      if (grepl("^[A-Za-z]:/", y)) y <- paste0("file:///", y)
      # UNC-Pfad //server/share/... -> zu file://server/share/...
      if (grepl("^//[^/]+/", y))   y <- sub("^//", "file://", y)
    }

    # 3) Percent-encode für Leerzeichen und weitere problematische Zeichen
    #    utils::URLencode(..., reserved = TRUE) encodiert u.a. Space -> %20
    y <- utils::URLencode(y, reserved = TRUE)  # vorsicht: nicht doppelt encoden!
    y
  }

  catalog <- normalizePath(system.file("extdata", "mzml_xsd", "catalog.xml", package = "MetabolomicsBasics"))

  # validate mzML
  with_local_catalog({
    # read mzML
    tmp <- xml2::read_xml(x = infile)
    ns <- c(mzml = "http://psi.hupo.org/ms/mzml")
    schema <- xml2::read_xml(system.file("extdata", "mzml_xsd", "mzML1.1.0_idx.xsd", package = "MetabolomicsBasics"))
    # find sourceFile/@location nodes and transform
    loc_nodes <- xml2::xml_find_all(tmp, "//mzml:sourceFile/@location", ns = ns)
    old_vals <- xml2::xml_text(loc_nodes)
    new_vals <- unname(vapply(old_vals, to_file_uri, character(1)))
    xml2::xml_set_attr(xml2::xml_parent(loc_nodes), "location", new_vals)
    stopifnot(isTRUE(xml2::xml_validate(tmp, schema)))
  }, catalog = catalog)

  # extract spectra
  # helper function to get precision and compression info per scan
  bda_info <- function(bda_node, ns) {
    acc <- function(a) xml2::xml_find_first(
      bda_node, sprintf(".//mzml:cvParam[@accession='%s']", a), ns = ns
    )
    is_mz   <- !is.na(acc("MS:1000514"))  # m/z array
    is_int  <- !is.na(acc("MS:1000515"))  # intensity array

    # precision
    is_64   <- !is.na(acc("MS:1000523"))  # 64-bit float
    is_32   <- !is.na(acc("MS:1000521"))  # 32-bit float

    # compression
    no_comp <- !is.na(acc("MS:1000576"))  # no compression
    zlib    <- !is.na(acc("MS:1000574"))  # zlib compression (falls gesetzt)

    # (optional) recognize MS-Numpress
    numpress <- !is.na(xml2::xml_find_first(
      bda_node,
      ".//mzml:cvParam[starts-with(@name,'MS-Numpress') or contains(@name,'numpress')]",
      ns = ns
    ))

    list(
      is_mz = is_mz, is_int = is_int,
      precision = if (is_64) "64" else if (is_32) "32" else NA_character_,
      compression = if (zlib) "zlib" else if (no_comp) "none" else "unknown",
      numpress = numpress
    )
  }

  # helper function to encode and decode binary
  decode_binary <- function(b_node, n, precision = c("32","64"), compression = c("none","zlib")) {
    precision   <- match.arg(precision)
    compression <- match.arg(compression)

    if (is.null(b_node) || length(b_node) == 0L || n <= 0L) return(numeric(0))

    btxt <- xml2::xml_text(b_node)
    if (is.na(btxt) || !nzchar(btxt)) return(numeric(0))

    raw0 <- base64enc::base64decode(btxt)
    if (length(raw0) == 0L) return(numeric(0))

    if (compression == "zlib") {
      raw0 <- memDecompress(raw0, type = "gzip")
      if (length(raw0) == 0L) return(numeric(0))
    }

    con <- rawConnection(raw0, "rb")
    on.exit(close(con), add = TRUE)
    size <- if (precision == "32") 4L else 8L
    readBin(con, what = "numeric", n = n, size = size, endian = "little")
  }
  encode_binary <- function(x, precision = c("32","64"), compression = c("none","zlib")) {
    precision  <- match.arg(precision)
    compression <- match.arg(compression)

    con <- rawConnection(raw(), "wb")
    on.exit(close(con))
    size <- if (precision == "32") 4L else 8L
    writeBin(as.numeric(x), con, size = size, endian = "little")
    raw0 <- rawConnectionValue(con)

    if (compression == "zlib") {
      raw0 <- memCompress(raw0, type = "gzip")
    }
    base64enc::base64encode(raw0)
  }

  # helper function to extract ms spectra
  safe_array_length <- function(sp, ns) {
    n <- suppressWarnings(as.integer(xml2::xml_attr(sp, "defaultArrayLength")))
    if (!is.na(n) && n > 0L) return(n)

    # Fallback: arrayLength an den binaryDataArray-Knoten
    bdas <- xml2::xml_find_all(sp, ".//mzml:binaryDataArrayList/mzml:binaryDataArray", ns = ns)
    if (length(bdas)) {
      n_bda <- suppressWarnings(as.integer(xml2::xml_attr(bdas, "arrayLength")))
      n_bda <- n_bda[!is.na(n_bda) & n_bda > 0L]
      if (length(n_bda)) return(n_bda[1L])
    }
    0L
  }
  extract_ms <- function(spec_nodes, ns) {
    out <- vector("list", length(spec_nodes))

    for (i in seq_along(spec_nodes)) {
      sp <- spec_nodes[[i]]
      n  <- safe_array_length(sp, ns)

      bdas <- xml2::xml_find_all(sp, ".//mzml:binaryDataArrayList/mzml:binaryDataArray", ns = ns)
      infos <- lapply(bdas, bda_info, ns = ns)

      # nodes & info for m/z and Int
      idx_mz  <- which(vapply(infos, `[[`, logical(1), "is_mz"))
      idx_int <- which(vapply(infos, `[[`, logical(1), "is_int"))

      # Default return
      mz  <- numeric(0)
      int <- numeric(0)

      if (n > 0L && length(idx_mz) && length(idx_int)) {
        if (any(vapply(infos, `[[`, logical(1), "numpress"))) {
          warning("Found MS-Numpress encoded arrays. This function will NOT decode Numpress.")
        }

        bnode_mz  <- xml2::xml_find_first(bdas[[idx_mz]],  ".//mzml:binary", ns = ns)
        bnode_int <- xml2::xml_find_first(bdas[[idx_int]], ".//mzml:binary", ns = ns)

        # Dank defensivem decode_binary() sind diese Aufrufe jetzt safe
        mz  <- decode_binary(bnode_mz,  n, infos[[idx_mz]]$precision,  infos[[idx_mz]]$compression)
        int <- decode_binary(bnode_int, n, infos[[idx_int]]$precision, infos[[idx_int]]$compression)
      }

      out[[i]] <- list(
        node = sp,
        n    = n,
        mz   = mz,
        int  = int,
        meta = list(
          mz_info  = if (length(idx_mz))  infos[[idx_mz]] else NULL,
          int_info = if (length(idx_int)) infos[[idx_int]] else NULL
        ),
        empty = (n == 0L || length(mz) == 0L || length(int) == 0L)
      )
    }

    out
  }

  # helper function to modify spectra by applying mz shift
  apply_mz_shift <- function(ms_list, shift_da = -0.019) {
    for (i in seq_along(ms_list)) {
      ms_list[[i]]$mz <- ms_list[[i]]$mz + shift_da
    }
    ms_list
  }

  # helper function to write modified spectra back to mzML
  write_back_ms <- function(ms_list, ns) {
    for (i in seq_along(ms_list)) {
      sp   <- ms_list[[i]]$node
      n    <- ms_list[[i]]$n
      bdas <- xml2::xml_find_all(sp, ".//mzml:binaryDataArrayList/mzml:binaryDataArray", ns = ns)
      infos<- lapply(bdas, bda_info, ns = ns)

      # m/z-Array ersetzen
      idx_mz <- which(vapply(infos, `[[`, logical(1), "is_mz"))
      bnode_mz <- xml2::xml_find_first(bdas[[idx_mz]], ".//mzml:binary", ns = ns)
      new64 <- encode_binary(ms_list[[i]]$mz, precision  = infos[[idx_mz]]$precision, compression= infos[[idx_mz]]$compression)
      xml2::xml_set_text(bnode_mz, new64)
      # encodedLength (falls vorhanden) anpassen
      xml2::xml_set_attr(xml2::xml_parent(bnode_mz), "encodedLength", as.character(nchar(new64)))

      # if one day we modify intensities
      # idx_int <- which(vapply(infos, `[[`, logical(1), "is_int"))
      # bnode_int <- xml2::xml_find_first(bdas[[idx_int]], ".//mzml:binary", ns = ns)
      # newI <- encode_binary(ms_list[[i]]$int, precision=..., compression=...)
      # xml2::xml_set_text(bnode_int, newI)
      # xml2::xml_set_attr(xml2::xml_parent(bnode_int), "encodedLength", as.character(nchar(newI)))
    }
    invisible(ms_list)
  }

  # helper function to update spectra dependent information (base peak and precursor)
  update_base_peak_mz_cvparam <- function(spectrum_node, ns, shift_da = 0, mode = c("shift", "compute")) {
    mode <- match.arg(mode)

    bp_nodes <- xml2::xml_find_all(spectrum_node, ".//mzml:cvParam[@accession='MS:1000504']", ns = ns)
    if (length(bp_nodes) == 0L) return(invisible(FALSE))

    if (mode == "shift") {
      old_vals <- suppressWarnings(as.numeric(xml2::xml_attr(bp_nodes, "value")))
      # only update those we can parse as numeric
      ok <- !is.na(old_vals)
      if (any(ok)) {
        new_vals <- old_vals[ok] + shift_da
        xml2::xml_set_attr(bp_nodes[ok], "value", format(new_vals, scientific = FALSE, digits = 15))
        return(invisible(TRUE))
      } else {
        return(invisible(FALSE))
      }
    } else {
      # Recompute (if Intensities have been modified)
      bdas <- xml2::xml_find_all(
        spectrum_node, ".//mzml:binaryDataArrayList/mzml:binaryDataArray", ns = ns
      )
      if (!length(bdas)) return(invisible(FALSE))

      infos <- lapply(bdas, bda_info, ns = ns)

      idx_mz  <- which(vapply(infos, `[[`, logical(1), "is_mz"))
      idx_int <- which(vapply(infos, `[[`, logical(1), "is_int"))
      if (!length(idx_mz) || !length(idx_int)) return(invisible(FALSE))

      n <- as.integer(xml2::xml_attr(spectrum_node, "defaultArrayLength"))
      if (is.na(n) || n <= 0) return(invisible(FALSE))

      bnode_int <- xml2::xml_find_first(bdas[[idx_int]], ".//mzml:binary", ns = ns)
      bnode_mz  <- xml2::xml_find_first(bdas[[idx_mz ]], ".//mzml:binary", ns = ns)

      int <- decode_binary(bnode_int, n, infos[[idx_int]]$precision, infos[[idx_int]]$compression)
      mz  <- decode_binary(bnode_mz,  n, infos[[idx_mz ]]$precision,  infos[[idx_mz ]]$compression)

      if (!length(int) || !length(mz)) return(invisible(FALSE))
      k <- which.max(int)
      bp_new <- mz[k]

      xml2::xml_set_attr(bp_nodes, "value", format(bp_new, scientific = FALSE, digits = 15))
      return(invisible(TRUE))
    }
  }
  update_precursor_mz_cvparams <- function(spectrum_node, ns, shift_da = 0) {
    # selected ion m/z (MS:1000744) – im <selectedIon>
    sel_nodes <- xml2::xml_find_all(spectrum_node, ".//mzml:precursor//mzml:selectedIon/mzml:cvParam[@accession='MS:1000744']", ns = ns)
    if (length(sel_nodes)) {
      vals <- suppressWarnings(as.numeric(xml2::xml_attr(sel_nodes, "value")))
      ok <- !is.na(vals)
      if (any(ok)) {
        xml2::xml_set_attr(sel_nodes[ok], "value", format(vals[ok] + shift_da, scientific = FALSE, digits = 15))
      }
    }
    # isolation window target m/z (MS:1000827) – im <isolationWindow>
    iso_nodes <- xml2::xml_find_all(spectrum_node, ".//mzml:precursor//mzml:isolationWindow/mzml:cvParam[@accession='MS:1000827']", ns = ns)
    if (length(iso_nodes)) {
      vals <- suppressWarnings(as.numeric(xml2::xml_attr(iso_nodes, "value")))
      ok <- !is.na(vals)
      if (any(ok)) {
        xml2::xml_set_attr(iso_nodes[ok], "value", format(vals[ok] + shift_da, scientific = FALSE, digits = 15))
      }
    }
    invisible(TRUE)
  }

  if (apply_to %in% c("all", "MS1")) {
    ms_specs <- xml2::xml_find_all(tmp, "//mzml:spectrum[mzml:cvParam[@accession='MS:1000511' and @value='1']]", ns = ns)
    ms_list <- extract_ms(ms_specs, ns)
    #length(ms_list); lengths(ms_list[[1]][c("mz","int")])
    ms_list_mod <- apply_mz_shift(ms_list, shift_da = shift_da)
    ms_list_mod <- write_back_ms(ms_list_mod, ns = ns)
    for (sp in ms_specs) {
      #print(attr(xml2::as_list(sp)[[8]], "value"))
      update_base_peak_mz_cvparam(sp, ns = ns, shift_da = shift_da, mode = "shift")
    }
  }
  if (apply_to %in% c("all", "MS2")) {
    ms_specs <- xml2::xml_find_all(tmp, "//mzml:spectrum[mzml:cvParam[@accession='MS:1000511' and @value='2']]", ns = ns)
    ms_list <- extract_ms(ms_specs, ns)
    #length(ms_list); lengths(ms_list[[1]][c("mz","int")])
    ms_list_mod <- apply_mz_shift(ms_list, shift_da = shift_da)
    ms_list_mod <- write_back_ms(ms_list_mod, ns = ns)
    for (sp in ms_specs) {
      #print(attr(xml2::as_list(sp)[[8]], "value"))
      update_precursor_mz_cvparams(sp, ns = ns, shift_da = shift_da)
      update_base_peak_mz_cvparam(sp, ns = ns, shift_da = shift_da, mode = "shift")
    }
  }

  # save data
  xml2::write_xml(tmp, outfile)

  invisible(tmp)
}
