#' Power function for matrix
#' 
#' @param a matrix to be powered.
#' @param alpha a number by which the given matrix is powered.
#' @param active_dim a num
#' @return a matrix powered by alpha
#' @export
matpower_1 <- function(a, alpha, active_dim = 1e-6) {
  .Call(`_FAGLasso_matpower_1`, a, alpha, active_dim)
}

#' Funtional principal component analysis
#' 
#' @param coefs list object with \code{coefs[[j]] = [x^j_{1:n}]} : coordinates of j-th explanatory functions.
#'        Hence, \code{dim( coefs[[j]] ) = c(m_j, n)}, where m_j is the number of basis and n is the number of observations.
#' @param GB list object with \code{GB[[j]]} is the j-th Gram matrix such as \code{fourierpen(basis, 0) or bsplinepen(basis, 0)}.
#' 
#' @return A list with the following components:
#'  \itemize{
#'    \item \code{pred}: ist object with \code{pred[[j]] = pca scores of j-th explanatory functions, dim(pred[[j]]) = (n, m_j)}.
#'    \item \code{eval}: list object with\code{eval[[j]]} = eigenvalues of j-th covariance operator.
#'    \item \code{GB}: list object with \code{GB[[j]]} = Gram Matrix for j-th explanatory functions.
#'    \item \code{mat}: list object with \code{mat[[j]] = [I]_B^C, operator changing basis system from the original(bspline or fourier)} to the eigenfunctions.
#'    \item \code{Sigma}:  list object with \code{Sigma[[j]]} is the coordinate representation of j-th covariance operator, whose basis system is eigenfunctions.
#'    \item \code{non_centered}: non_centered design matrix based on eigenfunctions.
#'  }
#'  @export
fpca_cpp <- function(coefs, GB) {
  .Call(`_FAGLasso_fpca_cpp`, coefs, GB)
}

#' Soft thresh function
#' @param x a numeric number.
#' @return \code{x * I(x > 0)}.
#' @export
soft_thresh <- function(x) {
  .Call(`_FAGLasso_soft_thresh`, x)
}

#' Optimization function
#' @param y a vector of response variable.
#' @param design a matrix of coefficients for explanatory functions.
#' @param cov_mat a list of coordinate representations for covariance operators.
#' @param basis_num  a vector of the numbers of basis functions to represent explanatory functions.
#' @param weights a vector of weights for penalties. If all elements of it is set to be one, it is the case for the non-adaptive group lasso.
#' @param lambda_1 a hyperparameter for the group lasso penalty.
#' @param lambda_2 a vector of hyperparameters for the roughness penalties.
#' @param Ps a diagonal block matrix of second level gram matrices. 
#' @param initN a vector of initial value. It defaults to be NULL.
#' @param centered a Boolean value. If true, y and design are centered. Otherwise, both y and design are not centered.
#' @param active_dim a threshold for inverse of matrix.
#' @param max_iter the maximum number of iterations allowed.
#' @param standardnorm If set to 'TRUE', the standard L2 norm is employed instead of our gamma norm. Note that it is just for the simulation study.
#' @export
f_agL <- function(y, design, cov_mat, basis_num, weights, lambda_1, lambda_2, Ps, initN, centered = FALSE, active_dim = 1e-6, max_iter = 1000L, standardnorm = FALSE) {
  .Call(`_FAGLasso_f_agL`, y, design, cov_mat, basis_num, weights, lambda_1, lambda_2, Ps, initN, centered, active_dim, max_iter, standardnorm)
}

