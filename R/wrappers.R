#' Matrix power function
#'
#' Computes the matrix power \eqn{A^\alpha} for a symmetric matrix.
#'
#' @param a A symmetric matrix.
#' @param alpha Power exponent.
#' @param active_dim Threshold to ignore small eigenvalues.
#'
#' @return A matrix raised to the power alpha.
#' @export
matpower_1 <- function(a, alpha, active_dim = 1e-6) {
  .Call(`_FAGLasso_matpower_1`, a, alpha, active_dim)
}

#' Functional PCA
#' 
#' @param coefs A list of coefficient matrices for explanatory functions.
#' @param GB A list of Gram matrices.
#'
#' @return A list containing scores, eigenvalues, and projections.
#' @export
fpca_cpp <- function(coefs, GB) {
  .Call(`_FAGLasso_fpca_cpp`, coefs, GB)
}

#' Soft-threshold function
#'
#' @param x A numeric value.
#' @return \code{x} if positive, 0 otherwise.
#' @export
soft_thresh <- function(x) {
  .Call(`_FAGLasso_soft_thresh`, x)
}

#' Functional adaptive group lasso estimator
#'
#' @param y Response vector.
#' @param design Design matrix.
#' @param cov_mat List of covariance matrices.
#' @param basis_num Number of basis functions per variable.
#' @param weights Weights for group lasso.
#' @param lambda_1 Tuning parameter for group lasso.
#' @param lambda_2 Vector of tuning parameters for roughness penalty.
#' @param Ps Block diagonal penalty matrix.
#' @param initN Initial value (optional).
#' @param centered Whether to center input.
#' @param active_dim Threshold for eigenvalue truncation.
#' @param max_iter Maximum number of iterations.
#' @param standardnorm Whether to use standard L2 norm.
#'
#' @return Estimated coefficients.
#' @export
f_agL <- function(y, design, cov_mat, basis_num, weights, lambda_1, lambda_2,
                  Ps, initN, centered = FALSE, active_dim = 1e-6, max_iter = 1000L, standardnorm = FALSE) {
  .Call(`_FAGLasso_f_agL`, y, design, cov_mat, basis_num, weights, lambda_1, lambda_2,
        Ps, initN, centered, active_dim, max_iter, standardnorm)
}
