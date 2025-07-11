#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;

//' Power function for matrix
//' 
//' @param a matrix to be powered.
//' @param alpha a number by which the given matrix is powered.
//' @param active_dim a num
//' 
//' @return a matrix powered by alpha
// [[Rcpp::export]]
arma::mat matpower_1(const arma::mat& a, double alpha, double active_dim = 1e-6) {
  // Symmetrize the matrix
  arma::mat symA = (a + a.t()) / 2;
  
  // Compute eigenvalues and eigenvectors
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, symA);
  
  
  // Select the relevant eigenvalues and eigenvectors
  arma::uvec idx = arma::find(eigval > active_dim);
  arma::vec filtered_eigval = arma::pow(eigval.elem(idx), alpha);
  arma::mat filtered_eigvec = eigvec.cols(idx);
  
  // Reconstruct the matrix
  arma::mat result = filtered_eigvec * arma::diagmat(filtered_eigval) * filtered_eigvec.t();
  
  return result;
}



//' Funtional principal component analysis
//' 
//' @param coefs list object with \code{coefs[[j]] = [x^j_{1:n}]} : coordinates of j-th explanatory functions.
//'        Hence, \code{dim( coefs[[j]] ) = c(m_j, n)}, where m_j is the number of basis and n is the number of observations.
//' @param GB list object with \code{GB[[j]]} is the j-th Gram matrix such as \code{fourierpen(basis, 0) or bsplinepen(basis, 0)}.
//' 
//' @return A list with the following components:
//'  \itemize{
//'    \item \code{pred}: ist object with \code{pred[[j]] = pca scores of j-th explanatory functions, dim(pred[[j]]) = (n, m_j)}.
//'    \item \code{eval}: list object with\code{eval[[j]]} = eigenvalues of j-th covariance operator.
//'    \item \code{GB}: list object with \code{GB[[j]]} = Gram Matrix for j-th explanatory functions.
//'    \item \code{mat}: list object with \code{mat[[j]] = [I]_B^C, operator changing basis system from the original(bspline or fourier)} to the eigenfunctions.
//'    \item \code{Sigma}:  list object with \code{Sigma[[j]]} is the coordinate representation of j-th covariance operator, whose basis system is eigenfunctions.
//'    \item \code{non_centered}: non_centered design matrix based on eigenfunctions.
//'  }
// [[Rcpp::export]]
List fpca_cpp(const List& coefs, const List& GB) {
  int p = coefs.size();
  int n = as<arma::mat>(coefs[0]).n_cols;
  
  arma::mat one = arma::ones<arma::mat>(n, 1);
  arma::mat Q = arma::eye(n, n) - one * one.t() / n;
  
  List Sigma(p), pred(p), eval(p), non_centered(p), mat(p);
  
  for (int j = 0; j < p; j++) {
    arma::mat coef = as<arma::mat>(coefs[j]);
    arma::mat GBj = as<arma::mat>(GB[j]);
    
    arma::mat B_half = matpower_1(GBj, 0.5);
    arma::mat Sigma_j = B_half * coef * Q * coef.t() * B_half / n;
    
    
    arma::vec egn_val;
    arma::mat egn_vec;
    eig_sym(egn_val, egn_vec, Sigma_j);
    egn_vec = fliplr(egn_vec);
    arma::mat B_inv_half = matpower_1(GBj, -0.5);
    arma::mat mat_j = B_inv_half * egn_vec;
    
    mat[j] = mat_j;
    non_centered[j] = coef.t() * GBj * mat_j;
    pred[j] = Q * coef.t() * GBj * mat_j;
    eval[j] = egn_val;
    Sigma[j] = mat_j.t() * GBj * Sigma_j * mat_j;
  }
  
  return List::create(
    Named("pred") = pred,
    Named("eval") = eval,
    Named("GB") = GB,
    Named("mat") = mat,
    Named("Sigma") = Sigma,
    Named("non_centered") = non_centered
  );
}


//' Soft thresh function
//' @param x a numeric number.
//' @return \code{x * I(x > 0)}.
// [[Rcpp::export]]
double soft_thresh(double x) {
  if (x <= 0) {
    return 0;
  } else {
    return x;
  }
}


//' Optimization function
//' @param y a vector of response variable.
//' @param design a matrix of coefficients for explanatory functions.
//' @param cov_mat a list of coordinate representations for covariance operators.
//' @param basis_num  a vector of the numbers of basis functions to represent explanatory functions.
//' @param weights a vector of weights for penalties. If all elements of it is set to be one, it is the case for the non-adaptive group lasso.
//' @param lambda_1 a hyperparameter for the group lasso penalty.
//' @param lambda_2 a vector of hyperparameters for the roughness penalties.
//' @param Ps a diagonal block matrix of second level gram matrices. 
//' @param initN a vector of initial value. It defaults to be NULL.
//' @param centered a Boolean value. If true, y and design are centered. Otherwise, both y and design are not centered.
//' @param active_dim a threshold for inverse of matrix.
//' @param max_iter the maximum number of iterations allowed.
//' @param standardnorm If set to 'TRUE', the standard L2 norm is employed instead of our gamma norm. Note that it is just for the simulation study.
// [[Rcpp::export]]
arma::vec f_agL(const arma::vec& y, 
                const arma::mat& design, 
                const Rcpp::List& cov_mat, 
                const arma::vec& basis_num, 
                const arma::vec& weights, 
                const double& lambda_1, 
                const arma::vec& lambda_2, 
                const arma::mat& Ps, 
                Rcpp::Nullable<Rcpp::NumericVector> initN,
                bool centered = false,
                double active_dim = 1e-6, 
                int max_iter = 1000, 
                bool standardnorm = false) {
  
  arma::vec y_centered = y;
  arma::mat x = design;
  if (!centered) {
    y_centered -= arma::mean(y);
    x.each_row() -= arma::mean(x, 0);
  }
  
  int n = design.n_rows;
  int p = basis_num.n_elem;
  int p_total = arma::sum(basis_num);
  
  arma::mat A_root_inv = arma::zeros<arma::mat>(p_total, p_total);
  arma::vec deltas = arma::zeros<arma::vec>(p);
  arma::vec start_vec = arma::cumsum(arma::join_vert(arma::vec({1}), basis_num.head(p - 1)));
  arma::vec end_vec = arma::cumsum(basis_num);
  double e_star = 1e-6;
  
  
  Rcpp::List trans_D_y(p);
  for (int j = 0; j < p; ++j) {
    int start = start_vec[j] - 1;
    int end = end_vec[j] - 1;
    arma::mat x_j = x.cols(start, end);
    arma::mat Gamma_j = Rcpp::as<arma::mat>(cov_mat[j]);
    arma::mat A_j;
    
    if (!standardnorm) {
      A_j = Gamma_j + lambda_2[j] * Ps.submat(start, start, end, end);
    } else {
      A_j = arma::eye<arma::mat>(basis_num[j], basis_num[j]) + lambda_2[j] * Ps.submat(start, start, end, end);
    }
    
    arma::mat A_j_root_inv = matpower_1(A_j, -0.5, active_dim);
    A_root_inv.submat(start, start, end, end) = A_j_root_inv;
    arma::mat D_j = x_j * A_j_root_inv;
    
    x.cols(start, end) = D_j;
    
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, (1.0 / n) * D_j.t() * D_j);
    deltas[j] = (1 + e_star) * eigval[eigval.n_elem - 1]; 
    trans_D_y[j] = D_j.t() * y;
  }
  
  Rcpp::List trans_D_j_D(p);
  for(int j = 0; j < p; j++){
    int start = start_vec(j) - 1; 
    int end = end_vec(j) - 1;
    arma::mat x_j = x.cols(start, end);
    trans_D_j_D[j] = x_j.t() * x;
  }
  
  arma::vec init;
  if(initN.isNotNull()){
    init = Rcpp::as<arma::vec>(initN);
  } else {
    init = arma::zeros<arma::vec>(p_total);
  }
  
  
  arma::vec beta_old = init;
  arma::vec beta_new = init;
  
  for (int i = 0; i < max_iter; ++i) {
    for (int j = 0; j < p; ++j) {
      double delta = deltas[j];
      int start = start_vec[j] - 1;
      int end = end_vec[j] - 1;
      arma::mat D_j = x.cols(start, end);
      arma::vec trans_D_y_j = Rcpp::as<arma::vec>(trans_D_y[j]);
      arma::mat trans_D_j_D_j = Rcpp::as<arma::mat>(trans_D_j_D[j]);
      arma::vec temp = trans_D_y_j / n - (trans_D_j_D_j * beta_new) / n + delta * beta_new.subvec(start, end);
      arma::vec new_j = (1.0 / delta) * temp * soft_thresh(1 - lambda_1 * weights[j] / sqrt(arma::accu(arma::square(temp))));
      beta_new.subvec(start, end) = new_j;
    }
    double loss_old = 0.5 / n * arma::accu(arma::square(y - x * beta_old));
    double loss_new = 0.5 / n * arma::accu(arma::square(y - x * beta_new));
    for (int j = 0; j < p; ++j) {
      int start = start_vec[j] - 1;
      int end = end_vec[j] - 1;
      loss_old += lambda_1 * weights[j] * sqrt(arma::accu(arma::square(beta_old.subvec(start, end))));
      loss_new += lambda_1 * weights[j] * sqrt(arma::accu(arma::square(beta_new.subvec(start, end))));
    }
    
    if (std::abs(loss_old - loss_new) < 1e-4) {
      break;
    } else {
      beta_old = beta_new;
    }
  }
  
  beta_new = A_root_inv * beta_new;
  return beta_new;
}

