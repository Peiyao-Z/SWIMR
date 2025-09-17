#' Select Informative Genes used in the deconvolution
#'
#' This function selects a set of informative genes for deconvolution
#' based on a combination of differential expression and gene dispersion
#' within cell types.
#'
#' @param Basis A reference basis matrix where columns are cell types
#'   and rows are genes. Values represent gene expression.
#' @param sc_eset A SingleCellExperiment S4 object containing scRNA-seq
#'   data and metadata.
#' @param commonGene A character vector of genes common to both
#'   the single-cell and spatial transcriptomics datasets.
#' @param ct.select A character vector of cell type names to be
#'   used for deconvolution. If NULL, all cell types are used.
#' @param ct.varname A character string specifying the column name
#'   in `colData(sc_eset)` that contains cell type annotations.
#' @importFrom SummarizedExperiment assays colData
#' @importFrom stats quantile
#' @return A character vector of selected informative genes.
#' @export
#'
selectInfo <- function(Basis, sc_eset, commonGene, ct.select = NULL, ct.varname) {
  # --- 1. Identify differentially expressed genes (log-fold change) ---
  # If no cell types are selected, use all available cell types from the basis matrix.
  if (is.null(ct.select)) {
    ct.select <- colnames(Basis)
  }

  # Helper function to calculate row means, handling single-column matrices
  # without errors.
  .rowMeans_safe <- function(mat) {
    if (ncol(mat) == 1) {
      return(mat[, 1])
    } else {
      return(rowMeans(mat))
    }
  }

  # Find genes with log mean fold-change > 0.25 and positive expression
  # in each selected cell type.
  fc_genes <- lapply(ct.select, function(ct) {
    # Calculate the mean expression of all other cell types.
    rest_of_cells <- Basis[, colnames(Basis) != ct, drop = FALSE]
    mean_rest <- .rowMeans_safe(rest_of_cells)

    # Calculate log-fold change. Add a small constant for numerical stability.
    log_fc <- log(Basis[, ct] + 1e-6) - log(mean_rest + 1e-6)

    # Select genes with a positive log-fold change and positive expression
    # in the current cell type.
    rownames(Basis)[log_fc > 0.25 & Basis[, ct] > 0]
  })

  # Combine results and filter for genes present in both datasets.
  fc_genes <- unique(unlist(fc_genes))
  fc_genes <- intersect(fc_genes, commonGene)

  # --- 2. Filter genes by low intra-cell type dispersion ---

  # Get the raw counts for the differentially expressed genes.
  counts <- assays(sc_eset)$counts
  counts <- counts[rownames(counts) %in% fc_genes, , drop = FALSE]

  # Only consider cell types with at least 2 cells for variance calculation.
  cell_counts <- table(colData(sc_eset)[, ct.varname])
  valid_cts <- names(cell_counts[cell_counts > 1])

  # Calculate the variance-to-mean ratio (a measure of dispersion) within each valid cell type.
  dispersion_within_ct <- sapply(valid_cts, function(ct) {
    ct_counts <- counts[, colData(sc_eset)[, ct.varname] == ct, drop = FALSE]
    apply(ct_counts, 1, function(gene_counts) {
      var_counts <- var(gene_counts)
      mean_counts <- mean(gene_counts)
      # Avoid division by zero
      if (mean_counts > 0) {
        return(var_counts / mean_counts)
      } else {
        return(NA)
      }
    })
  })

  # Remove genes that are outliers in terms of dispersion.
  # This step filters out genes with high noise or high intrinsic variation.
  mean_dispersion <- apply(dispersion_within_ct, 1, mean, na.rm = TRUE)
  upper_quantile <- quantile(mean_dispersion, prob = 0.99, na.rm = TRUE)

  selected_genes <- names(mean_dispersion)[mean_dispersion < upper_quantile]

  return(selected_genes)
}


compute_kernel_subset <- function(spatial_location, isigma = 0.1) {
  # spatial_location rows must be in the desired spot order
  norm_c <- spatial_location[, c("x","y")]
  norm_c$x <- norm_c$x - min(norm_c$x)
  norm_c$y <- norm_c$y - min(norm_c$y)
  sf <- max(norm_c$x, norm_c$y)
  if (sf > 0) {
    norm_c$x <- norm_c$x / sf
    norm_c$y <- norm_c$y / sf
  }
  ED <- fields::rdist(as.matrix(norm_c))
  K  <- exp(-ED^2 / (2 * isigma^2))
  diag(K) <- 0
  dimnames(K) <- list(rownames(spatial_location), rownames(spatial_location))
  K
}


#' Run multi-reference deconvolution on a SWIMR object
#'
#' @param obj A \code{SWIMR} object containing:
#' \itemize{
#'   \item \code{sc_eset_list}: list of filtered scRNA-seq references
#'   \item \code{spatial_count}: filtered spatial count matrix
#'   \item \code{spatial_location}: matched coordinates data frame
#'   \item \code{info_parameters}: list of key parameters used for construction
#' }
#'
#' @return A \code{SWIMR} object with updated slots:
#' \itemize{
#'   \item \code{Weight_est}: reference weights per spatial location
#'   \item \code{Proportion_est}: fused cell type proportions per spatial location
#'   \item \code{algorithm_matrix}: intermediates from CARD
#' }
#' @importFrom SummarizedExperiment assays
#' @importFrom MCMCpack rdirichlet
#' @export
#'
#'
SWIMR_deconvolution <- function(obj) {

  phi     <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  epsilon <- 1e-04
  isigma  <- 0.1

  ## containers
  single_res <- vector("list", length(obj@sc_eset_list))
  names(single_res) <- names(obj@sc_eset_list)

  ## ---------- per-reference deconvolution ----------
  for (i in seq_along(obj@sc_eset_list)) {
    sc_eset <- obj@sc_eset_list[[i]]

    ct.select      <- obj@info_parameters$ct.select_list[[i]]
    ct.varname     <- obj@info_parameters$ct.varname_list[[i]]
    sample.varname <- obj@info_parameters$sample.varname_list[[i]]

    Basis_ref <- CARD::createscRef(sc_eset, ct.select, ct.varname, sample.varname)
    Basis <- Basis_ref$basis

    ## keep only requested cell types (preserve matrix-ness)
    keep <- colnames(Basis) %in% ct.select
    Basis <- Basis[, keep, drop = FALSE]
    idx <- match(ct.select, colnames(Basis))
    idx <- idx[!is.na(idx)]
    Basis <- Basis[, idx, drop = FALSE]

    spatial_count <- obj@spatial_count

    ## common genes (drop mitochondrial if named "mt-")
    commonGene <- intersect(rownames(spatial_count), rownames(Basis))
    commonGene <- commonGene[!(commonGene %in% commonGene[grep("mt-", commonGene)])]

    common <- selectInfo(Basis, sc_eset, commonGene, ct.select, ct.varname)

    Yinput <- spatial_count
    S      <- Basis

    ## order by gene names before intersecting
    Yinput <- Yinput[order(rownames(Yinput)), , drop = FALSE]
    S      <- S[order(rownames(S)), , drop = FALSE]

    S      <- S[rownames(S) %in% common, , drop = FALSE]
    Yinput <- Yinput[rownames(Yinput) %in% common, , drop = FALSE]

    ## guard: enough dimensions?
    if (nrow(Yinput) < 2 || ncol(Yinput) < 2 || nrow(S) < 1 || ncol(S) < 1) {
      stop(sprintf(
        "After filtering: Yinput is %d x %d; S is %d x %d.
Not enough genes/cell types survived. Check `selectInfo` thresholds or `ct.select`.",
        nrow(Yinput), ncol(Yinput), nrow(S), ncol(S)
      ))
    }

    ## remove all-zero rows/cols (keep 2-D)
    Yinput <- Yinput[Matrix::rowSums(Yinput) > 0, , drop = FALSE]
    Yinput <- Yinput[, Matrix::colSums(Yinput) > 0, drop = FALSE]

    colsumvec   <- Matrix::colSums(Yinput)
    Yinput_norm <- sweep(as.matrix(Yinput), 2, colsumvec, "/")

    ## align S rows to Yinput
    S <- S[rownames(S) %in% rownames(Yinput_norm), , drop = FALSE]
    S <- S[match(rownames(Yinput_norm), rownames(S)), , drop = FALSE]

    ## subset spatial_location and compute kernel
    spatial_location <- obj@spatial_location
    spatial_location <- spatial_location[rownames(spatial_location) %in% colnames(Yinput_norm), , drop = FALSE]
    spatial_location <- spatial_location[match(colnames(Yinput_norm), rownames(spatial_location)), , drop = FALSE]

    kernel_mat <- compute_kernel_subset(spatial_location, isigma)

    set.seed(20200107)
    P_int <- as.matrix(MCMCpack::rdirichlet(ncol(Yinput_norm), rep(10, ncol(S))))
    colnames(P_int) <- colnames(S)
    rownames(P_int) <- colnames(Yinput_norm)

    mean_Y <- mean(Yinput_norm)
    mean_S <- mean(S)
    Yinput_norm <- Yinput_norm * 1e-01 / mean_Y
    S <- S * 1e-01 / mean_S

    stopifnot(identical(rownames(S), rownames(Yinput_norm)))
    stopifnot(identical(colnames(S), colnames(P_int)))
    stopifnot(identical(colnames(Yinput_norm), rownames(kernel_mat)))

    ResList <- vector("list", length(phi))
    Obj     <- numeric(length(phi))

    for (iphi in seq_along(phi)) {
      res <- CARD::CARDref(
        XinputIn      = as.matrix(Yinput_norm),
        UIn           = as.matrix(S),
        WIn           = kernel_mat,
        phiIn         = phi[iphi],
        max_iterIn    = 1000,
        epsilonIn     = epsilon,
        initV         = P_int,
        initb         = rep(0, ncol(S)),
        initSigma_e2  = 0.1,
        initLambda    = rep(10, ncol(S))
      )
      rownames(res$V) <- colnames(Yinput_norm)
      colnames(res$V) <- colnames(S)
      ResList[[iphi]] <- res
      Obj[iphi] <- res$Obj
    }

    Optimal    <- max(which(Obj == max(Obj)))
    OptimalRes <- ResList[[Optimal]]

    # # hashes for reproducibility
    # h <- function(x) digest::digest(as.matrix(x), "xxhash64")
    # meta <- list(
    #   hash_X = h(Yinput_norm),
    #   hash_S = h(S),
    #   hash_W = h(kernel_mat),
    #   hash_V0 = h(P_int),
    #   phi = phi[Optimal]
    # )

    single_res[[i]] <- list(
      Proportion_single = sweep(OptimalRes$V, 1, rowSums(OptimalRes$V), "/"),
      scRef             = S * mean_S / 1e-01,
      spatial_location  = spatial_location
    )
  }

  ## ---------- multiple references ----------
  spatial_count    <- obj@spatial_count
  spatial_location <- obj@spatial_location

  tmp_name <- lapply(single_res, function(x) colnames(x$Proportion_single))
  cells    <- Reduce(intersect, tmp_name)

  ## reorder so shared cells come first
  for (i in seq_along(single_res)) {
    ps <- single_res[[i]]$Proportion_single
    ordered <- cbind(ps[, cells, drop = FALSE],
                     ps[, setdiff(colnames(ps), cells), drop = FALSE])
    colnames(ordered) <- c(cells, setdiff(colnames(ps), cells))
    single_res[[i]]$Proportion_single <- sweep(ordered, 1, rowSums(ordered), "/")
  }

  ## union of genes across references
  common <- Reduce(union, lapply(single_res, function(x) rownames(x$scRef)))

  Yinput      <- spatial_count
  colsumvec   <- Matrix::colSums(Yinput)
  Yinput_norm <- sweep(as.matrix(Yinput), 2, colsumvec, "/")
  Yinput_norm <- Yinput_norm[common, , drop = FALSE]

  spatial_location <- spatial_location[rownames(spatial_location) %in% colnames(Yinput_norm), , drop = FALSE]
  spatial_location <- spatial_location[match(colnames(Yinput_norm), rownames(spatial_location)), , drop = FALSE]

  Yinput_norm <- Yinput_norm[, match(rownames(single_res[[1]]$Proportion_single), colnames(Yinput_norm)), drop = FALSE]
  spatial_location <- spatial_location[match(rownames(single_res[[1]]$Proportion_single), rownames(spatial_location)), , drop = FALSE]

  ## recompute kernel for final spot order
  kernel_mat <- compute_kernel_subset(spatial_location, isigma)

  ## intersect genes present in all scRefs
  tmp_list <- lapply(single_res, function(x) intersect(common, rownames(x$scRef)))
  common2  <- Reduce(intersect, tmp_list)

  SListIn   <- lapply(single_res, function(x) x$scRef[common2, , drop = FALSE])
  initPList <- lapply(single_res, function(x) x$Proportion_single)
  Yinput_norm <- Yinput_norm[common2, , drop = FALSE]

  ## init Z (spots x references)
  Z <- matrix(1, nrow = ncol(Yinput_norm), ncol = length(single_res),
              dimnames = list(colnames(Yinput_norm), names(single_res)))

  initsMatList   <- lapply(single_res, function(x) colMeans(x$Proportion_single))
  initLambdaList <- lapply(SListIn, function(S) rep(10, ncol(S)))

  ## rescale
  mean_Y <- mean(as.matrix(Yinput_norm))
  Yinput_norm <- Yinput_norm * 1e-01 / mean_Y
  S_means <- vapply(SListIn, mean, numeric(1))
  for (i in seq_along(SListIn)) {
    SListIn[[i]] <- SListIn[[i]] * 1e-01 / S_means[i]
  }

  ## multi-ref step (call the R wrapper that binds to C++)
  result <- SWIMRMultiRef(
    YinputIn       = as.matrix(Yinput_norm),
    SListIn        = SListIn,
    KIn            = as.matrix(kernel_mat),
    phi1In         = 0.99,
    phi2In         = 0.99,
    max_iterIn     = 1000,
    epsilonIn      = epsilon,
    initPList      = initPList,      # NOTE: name matches wrapper formal
    initsMatList   = initsMatList,
    initSigma_e2   = 0.1,
    initLambdaList = initLambdaList,
    initZ          = Z
  )

  ## tidy names for PList
  for (i in seq_along(result$PList)) {
    rownames(result$PList[[i]]) <- colnames(Yinput_norm)
    colnames(result$PList[[i]]) <- colnames(SListIn[[i]])
  }

  ## assemble global proportions
  V_est <- result$Z
  rownames(V_est) <- colnames(Yinput_norm)
  colnames(V_est) <- names(single_res)

  total_name <- Reduce(union, tmp_name)
  final_V <- matrix(0, nrow = nrow(single_res[[1]]$Proportion_single),
                    ncol = length(total_name),
                    dimnames = list(rownames(single_res[[1]]$Proportion_single),
                                    total_name))

  for (nm in total_name) {
    for (i in seq_along(single_res)) {
      ps <- single_res[[i]]$Proportion_single
      if (tolower(nm) %in% tolower(colnames(ps))) {
        idx <- match(tolower(nm), tolower(colnames(ps)))
        final_V[, nm] <- final_V[, nm] + (t(ps) %*% diag(as.vector(V_est[, i])))[idx, ]
      }
    }
  }

  ## row-normalize safely
  rs <- rowSums(final_V)
  rs[rs == 0] <- 1
  final_V <- final_V / rs

  rs2 <- rowSums(V_est)
  rs2[rs2 == 0] <- 1
  V_est <- V_est / rs2

  ## write back
  obj@Weight_est        <- V_est
  obj@Proportion_est    <- final_V
  obj@Proportion_single <- result$PList

  return(obj)
}


