#' Each SWIMR object has a number of slots which store information. Key slots to access are listed below.
#'
#' @slot sc_eset_list A list of filtered scRNA-seq datasets, together with
#'   associated metadata, stored as SingleCellExperiment objects.
#' @slot spatial_count The filtered spatial transcriptomics count data.
#' @slot spatial_location A data frame specifying spatial coordinates for each location.
#' @slot Weight_est A matrix of estimated weights, where rows correspond to spatial
#'   locations and columns to scRNA-seq references.
#' @slot Proportion_est A matrix of estimated cell-type proportions from SWIMR,
#'   with rows representing spatial locations and columns representing cell types.
#' @slot Proportion_single A list of estimated cell-type proportions from SWIMR from each reference.
#' @slot info_parameters A list of parameters used during model fitting.
#'
setClass("SWIMR",
         slots = list(
           sc_eset_list = "list",
           spatial_count = "ANY",
           spatial_location = "data.frame",
           Weight_est = "matrix",
           Proportion_est = "matrix",
           Proportion_single = "list",
           info_parameters = "list")
)


#' Create a SWIMR object
#'
#' @param scRNAseq A **named or unnamed list** of references; each element must have
#'        `countData` (genes x cells) and `metaData` (cells x *).
#' @param spatial_count Gene-by-spot matrix (dense or sparse).
#' @param spatial_location data.frame/matrix with rownames matching colnames(spatial_count).
#' @param ct.varname_list Character vector/list giving the cell-type column in each ref meta.
#' @param ct.select_list Optional list of cell-type subsets for each ref.
#' @param sample.varname_list Optional list of sample column names for each ref.
#' @param minCountGene Minimum UMIs per spot (after row filter). Default 100.
#' @param minCountSpot Minimum nonzero spots per gene. Default 5.
#' @param maxCountGene Upper bound on UMIs per spot (to drop outliers). Default 1e6.
#' @param bypass_filters If TRUE, skips internal spatial filtering (assume pre-filtered).
#' @param verbose Print progress messages. Default TRUE.
#'
#' @return An S4 SWIMR object
#' @export
createSWIMRObject <- function(
    scRNAseq,
    spatial_count,
    spatial_location,
    ct.varname_list,
    ct.select_list = NULL,
    sample.varname_list = NULL,
    minCountGene = 100,
    minCountSpot = 5,
    maxCountGene = 1e6,
    bypass_filters = FALSE,
    verbose = TRUE
) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required.", call. = FALSE)
  if (!requireNamespace("CARD", quietly = TRUE))
    stop("Package 'CARD' is required.", call. = FALSE)

  msg <- function(...) if (isTRUE(verbose)) cat(sprintf(...))

  # --- helpers ---
  .as_dgc <- function(x) {
    if (inherits(x, "dgCMatrix")) return(x)
    if (inherits(x, "sparseMatrix")) return(as(x, "dgCMatrix"))
    if (is.matrix(x)) return(as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix"))
    if (is.vector(x)) return(as(Matrix::Matrix(t(as.matrix(x)), sparse = TRUE), "dgCMatrix"))
    stop("Counts must be vector, matrix, or sparseMatrix.", call. = FALSE)
  }

  # --- scRNA-seq references QC via CARD::sc_QC ---
  sc_eset_list <- list()
  if (is.null(names(scRNAseq))) names(scRNAseq) <- paste0("ref", seq_along(scRNAseq))

  for (i in seq_along(scRNAseq)) {
    ref_name <- names(scRNAseq)[i]
    ref <- scRNAseq[[i]]
    if (!all(c("countData", "metaData") %in% names(ref))) {
      stop(sprintf("Reference '%s' must contain 'countData' and 'metaData'.", ref_name))
    }
    sc_count <- ref[["countData"]]
    sc_meta  <- ref[["metaData"]]

    sample.varname <- if (length(sample.varname_list) >= i) sample.varname_list[[i]] else NULL
    ct.varname     <- if (length(ct.varname_list)     >= i) ct.varname_list[[i]]     else NULL
    ct.select      <- if (!is.null(ct.select_list) && length(ct.select_list) >= i) ct.select_list[[i]] else NULL

    sc_countMat <- .as_dgc(sc_count)

    if (is.null(rownames(sc_countMat)) || any(rownames(sc_countMat) == ""))
      stop(sprintf("Empty or missing gene names in '%s' countData.", ref_name), call. = FALSE)
    if (is.null(colnames(sc_countMat)))
      stop(sprintf("Missing cell names (colnames) in '%s' countData.", ref_name), call. = FALSE)
    if (is.null(rownames(sc_meta)))
      stop(sprintf("Missing rownames (cell IDs) in '%s' metaData.", ref_name), call. = FALSE)

    if (!identical(sort(rownames(sc_meta)), sort(colnames(sc_countMat)))) {
      stop(sprintf("Cell names mismatch between '%s' countData (cols) and metaData (rows).", ref_name))
    }
    # Reorder meta to countData columns
    sc_meta <- sc_meta[colnames(sc_countMat), , drop = FALSE]

    if (is.null(sample.varname)) {
      sample.varname <- "sampleID"
      sc_meta <- as.data.frame(sc_meta)
      sc_meta[[sample.varname]] <- "Sample"
    }
    if (is.null(ct.varname) || !ct.varname %in% colnames(sc_meta)) {
      stop(sprintf("Provide a valid 'ct.varname' present in '%s' metaData.", ref_name))
    }
    if (is.null(ct.select)) {
      ct.sel <- unique(as.character(sc_meta[[ct.varname]]))
    } else {
      ct.sel <- as.character(ct.select[!is.na(ct.select)])
    }

    sc_eset_list[[ref_name]] <- CARD::sc_QC(
      sc_countMat, sc_meta,
      ct.varname     = ct.varname,
      ct.select      = ct.sel,
      sample.varname = sample.varname
    )
  }

  # --- spatial dataset QC ---

  spatial_countMat <- .as_dgc(spatial_count)

  if (is.null(rownames(spatial_countMat)) || any(rownames(spatial_countMat) == ""))
    stop("Empty or missing gene rownames in spatial_count.", call. = FALSE)

  if (is.null(spatial_location))
    stop("Please provide the matching spatial_location data frame.", call. = FALSE)

  if (!is.data.frame(spatial_location) && !is.matrix(spatial_location))
    stop("spatial_location must be a data.frame or matrix with rownames.", call. = FALSE)

  spatial_location <- as.data.frame(spatial_location, stringsAsFactors = FALSE)

  if (is.null(rownames(spatial_location)))
    stop("spatial_location must have rownames matching spot/barcode IDs.", call. = FALSE)

  if (ncol(spatial_countMat) != nrow(spatial_location)) {
    stop("Mismatch in spatial locations: spatial_count is genes x n; spatial_location is n x p.", call. = FALSE)
  }
  if (!identical(colnames(spatial_countMat), rownames(spatial_location))) {
    stop("Column names of spatial_count must match row names of spatial_location (same order).", call. = FALSE)
  }

  if (!isTRUE(bypass_filters)) {
    # robust row filter: nonzeros per gene
    g_ok <- Matrix::rowSums(spatial_countMat != 0) >= minCountSpot
    if (!any(g_ok)) {
      stop(sprintf("No genes pass minCountSpot=%d. Consider lowering the threshold.", minCountSpot), call. = FALSE)
    }
    spatial_countMat <- spatial_countMat[g_ok, , drop = FALSE]

    # robust col filter: total UMIs per spot
    cs <- Matrix::colSums(spatial_countMat)
    s_ok <- (cs >= minCountGene) & (cs <= maxCountGene)
    if (!any(s_ok)) {
      stop(sprintf("No spots pass minCountGene=%d and <= %g. Consider adjusting thresholds.",
                   minCountGene, maxCountGene), call. = FALSE)
    }
    spatial_countMat  <- spatial_countMat[, s_ok, drop = FALSE]
    spatial_location  <- spatial_location[colnames(spatial_countMat), , drop = FALSE]
  } else {
    msg("Bypassing internal spatial filters (assuming pre-filtered inputs).\n")
  }

  # --- construct SWIMR object ---
  n_spots <- ncol(spatial_countMat)
  Weight_est <- matrix(NA_real_, nrow = n_spots, ncol = length(sc_eset_list))
  Proportion_est <- matrix(NA_real_, nrow = n_spots, ncol = 0)

  object <- new(
    Class = "SWIMR",
    sc_eset_list      = sc_eset_list,
    spatial_count     = spatial_countMat,
    spatial_location  = spatial_location,
    Weight_est        = Weight_est,
    Proportion_est    = Proportion_est,
    Proportion_single = list(),
    info_parameters   = list(
      ct.varname_list     = ct.varname_list,
      ct.select_list      = ct.select_list,
      sample.varname_list = sample.varname_list,
      minCountGene        = minCountGene,
      minCountSpot        = minCountSpot,
      maxCountGene        = maxCountGene,
      bypass_filters      = bypass_filters
    )
  )

  return(object)
}
