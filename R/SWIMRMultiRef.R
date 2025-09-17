#' Multi-reference spatial deconvolution (SWIMR C++ core)
#'
#' Runs the SWIMR multi-reference Conditional Autoregressive (CAR) model to
#' jointly estimate cell-type compositions per reference and location weights.
#'
#' @param YinputIn A numeric matrix (genes × spots): normalized spatial counts.
#' @param SListIn A list of basis matrices `S_r` (genes × K_r), one per reference.
#' @param KIn A numeric matrix (spots × spots): spatial affinity/adjacency (e.g., Gaussian kernel).
#' @param phi1In,phi2In Numeric scalars: CAR strength parameters for P and Z, respectively.
#' @param max_iterIn Integer: maximum iterations.
#' @param epsilonIn Numeric: relative tolerance for convergence check.
#' @param initPList A list of initial proportion matrices `P_r` (spots × K_r).
#' @param initsMatList A list of initial intercept vectors `s_rk` (length K_r) per reference.
#' @param initSigma_e2 Numeric scalar: initial residual variance.
#' @param initLambdaList A list of initial smoothing scales `lambda_rk` (length K_r) per reference.
#' @param initZ A numeric matrix (spots × R): initial location weights per reference.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{PList}: updated list of cell-type proportions per reference.
#'   \item \code{Z}: updated (spots × R) location weights.
#'   \item \code{sigma_e2}: estimated residual variance.
#'   \item \code{lambda}: updated list of \eqn{\lambda_{rk}} per reference.
#'   \item \code{sMatList}: updated list of intercept vectors per reference.
#'   \item \code{Obj}: final objective value.
#' }
#'
#' @export
SWIMRMultiRef <- function(YinputIn, SListIn, KIn, phi1In, phi2In, max_iterIn,
                          epsilonIn, initPList, initsMatList, initSigma_e2,
                          initLambdaList, initZ) {
  .Call(`_SWIMR_SWIMRMultiRef`, YinputIn, SListIn, KIn, phi1In, phi2In,
        max_iterIn, epsilonIn, initPList, initsMatList, initSigma_e2,
        initLambdaList, initZ)
}
