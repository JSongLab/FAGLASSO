#' @useDynLib FAGLasso, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL



#' numerical integration
#'
#' @param x vector of f(t_i), tt: vector of t_i
#' @param tt integration of f on the range of tt
#'
#' @return integration of f on the range of tt
#' @export
#'
#' @examples
numeric_integ <- function(x, tt){
  base <- (tt[2] - tt[1] + tt[4] - tt[3])/2
  nt <- length(tt)
  x[2:(nt-1)] <- 2*x[2:(nt-1)]
  return(sum(x)*base/2)
}




#' Funtional Adaptive Group Lasso function.
#'
#' @param y vector of response variable.
#' @param design  design matrix \code{[ [x^1_1:n]^T ... [x^p_1:n]^T ] in R^{ n x (m_1 + ... m_j) }}.
#'        Note that there is no gram matrix, so, we assume that orthonormal basis system is used such as fourier or eigenfunctions.
#' @param Ps A matris of roughness penalties, structured as diag(P_1, ... , P_p), where each P_i represents a specific penalty matrix.
#'           For instance, P_i could be 'bsplinepen(..., Lfdobj=2)'.
#' @param basis_num A vector of the number of basis functions used; (m_1, ...., m_p).
#' @param fpca_basis_num A vector indicating the number of eigenfunctions of FPCA for the basis system. It is allowed to differ from 'basis_num'.
#' @param basis_list A list where \code{'basis_list[[j]]'} contains the information for the basis used to represent the j-th explanatory function.
#'                   Therefore, \code{'basis_list[[j]]'} should be populated with basis creation functions such as 'create.bsplinebasis(...)' or 'create.fourierbasis(...)'.
#'                   If 'fpca' is set to 'TRUE', this argument is required.
#' @param weights weights of group lasso penalties. If you want to conduct non-adaptive group lasso, then set all weights to one.
#' @param fpca It has a Boolean value. If se to 'True', a basis system based on FPCA is used. Otherwise, it is assumed that the Fourier basis system will be employed.
#' @param lambda_1 It defines the search range for group lasso penalty, specified as a vector in the format `c(min, max)`. If set to `NULL`, it defaults to `c(1e-5, 1e-1)`.
#' @param lambda_2  It defines the search range for each roughness penalty. It is represented as a matrix with p rows and 2 columns, where each row corresponds to a specific roughness penalty.
#'                  Specifically, \code{lambda_2[j,]} contains the range with the form of c(min, max) for the j-th roughness penalty.
#' @param n_lamb_1 The number of grid for the range specified in lambda_1.
#' @param n_lamb_2 The common number of grid for the range specified in lambda_2.
#' @param cv It specifies the number of folds for K-fold cross-validation used to identify optival values for 'lambda_1' and 'lambda_2s'.
#'           Here, 'optimal' refers to minimizing the test MSE in each fold.
#' @param cv_average It has a Boolean value (True/False). If set to true, the final result is the average of the estimated functions across each fold.
#'                   Otherwise, the weights determined through cross-validation are applied to the entire trainning dataset to produce the final result.
#' @param active_dim This is identical argument of 'matpower'.
#' @param add_iter It specifies the number of additional iterations for optimizing 'lambda_2s'. By default, it is set to 'NULL', indicating no additional iterations.
#'                 Otherwise, it accepts a positive integer. If 'add_iter=k', then fixing the group lasso penalty value('lambda_1'), k additional optimization processes are carried out for 'lambda_2s'.
#' @param init Initial value of optimization, if null, optimization starts at (0, ..... , 0).
#' @param max_iter the number of max iteration.
#' @param standardnorm If set to 'TRUE', the standard L2 norm is employed instead of our gamma norm. Note that it is just for the simulation study.
#'
#' @return A list with the following components:
#'   \itemize{
#'     \item \code{beta_hat}: the coordinates of estimated functions. It is vector with its length being (m_1 + ... m_p).
#'     \item \code{lambda_1}: lambda_1 selected by cross validation.
#'     \item \code{lambda_2}: lambda_2s selected by cross validation.
#'     \item \code{cov_mat}: list object with \code{cov_mat[[j]]}=The coordinate representation of j-th covariance operators,
#'                           whose basis system is identical to the design specifed in input argument such as 'fourier', 'bspline'.
#'                           This output can be used to calculate \code{w_{1j}}s.
#'   }
#'
#' @export
#'
#' @examples
Fun_Reg_AGL <- function(y,
                        design,
                        Ps,
                        basis_num,
                        fpca_basis_num,
                        basis_list=NULL,
                        weights,
                        fpca=FALSE,
                        lambda_1=NULL,
                        lambda_2=NULL,
                        n_lamb_1=20,
                        n_lamb_2=10,
                        cv=5,
                        cv_average=FALSE,
                        active_dim=1e-6,
                        add_iter=NULL,
                        init=NULL,
                        max_iter=1000,
                        standardnorm=FALSE){
  
  n <- dim(design)[1]
  
  # If fpca==FALSE, then Fourier basis system is used.
  if(fpca == FALSE){
    ww <- weights
    p <- length(basis_num)
    
    # centering the design matrix by column.
    design_centered <- t(t(design) - apply(design, 2, mean))
    
    # list for the coordinate representation of covariance operator.
    cov_mat <- list()
    
    start_vec <- cumsum(c(1, basis_num))[1:p]
    end_vec <- cumsum(basis_num)
    for(j in 1:p){
      start <- start_vec[j]
      end <- end_vec[j]
      x_j <- design_centered[, start:end]
      # calculate the coordinate representation for the j-th sample covariance operator.
      cov_mat[[j]] <- t(x_j) %*% x_j / (n-1)
    }
    
    x_train_k <- list()
    x_test_k <- list()
    y_train_k <- list()
    y_test_k <- list()
    
    margin <- ceiling(n/cv)
    
    # Divide total data into 'cv' folds.
    for(i in 1:cv){
      start <- 1 + margin*(i-1)
      end <- min(margin*i, n)
      idx <- (1:n) %in% (start:end)
      
      x_train_k[[i]] <- design[!idx,]
      x_test_k[[i]] <- t(t(design[idx,]) - apply(design[!idx,], 2, mean))
      y_train_k[[i]] <- y[!idx]
      y_test_k[[i]] <- y[idx] - mean(y[!idx])
    }
    
    # Assign the searching range of the group lasso penalty.
    if(is.null(lambda_1)){
      lambdas_1 <- exp(seq(log(1e-5), log(1e-1), length=n_lamb_1))
    }else{
      lambdas_1 <- exp(seq(log(lambda_1[1]), log(lambda_1[2]), length=n_lamb_1))
    }
    
    
    
    # Assign searching ranges of roughness penalties.
    # If the standard L2 norm is used for the group lasso penalty, all weights of roughness penalty are the same.
    if(standardnorm == FALSE){
      if(is.null(lambda_2)){
        lambdas_2 <- matrix(exp(seq(log(1e-4), log(1e-0), length=n_lamb_2)), p, n_lamb_2, byrow=T)
      }else{
        lambdas_2 <- matrix(0, p, n_lamb_2)
        for(j in 1:p){
          lambdas_2[j,] <- exp(seq(log(lambda_2[j, 1]), log(lambda_2[j, 2]), length=n_lamb_2))
        }
      }
    }
    if(standardnorm == TRUE){
      if(is.null(lambda_2)){
        lambdas_2 <- exp(seq(log(1e-9), log(1e-5), length=n_lamb_2))
      }else{
        lambdas_2 <- exp(seq(log(lambda_2[1]), log(lambda_2[2]), length=n_lamb_2))
      }
    }
    
    
    # It is just for the progress bar.
    if(standardnorm == FALSE){
      total_iterations <- n_lamb_1 * p * n_lamb_2
    }
    if(standardnorm == TRUE){
      total_iterations <- n_lamb_1 * n_lamb_2
    }
    iteration_counter <- 0
    
    
    
    # Using our FC distance based norm, conduct the functional linear regression.
    if(standardnorm == FALSE){
      # Search the weight of the group Lasso penalty.
      for(i in 1:n_lamb_1){
        # Initialize weights of roughness penalties as the median value of `lambda_2` for each row.
        lambdas_2_vec <- apply(lambdas_2, 1, stats::median)
        # Search weights of roughness penalties, fixing lambda_{1i}.
        for(j in 1:p){
          if(j > 1){
            lambdas_2_vec <- best_lambda_2
          }
          # Search j-th weight of roughness penalties, fixing lambda_{1i}, lambda_{21}, ... , lambda_{2, (j-1)}.
          for(k in 1:n_lamb_2){
            iteration_counter <- iteration_counter + 1
            percent_complete <- (iteration_counter / total_iterations) * 100
            lambdas_2_vec[j] <- lambdas_2[j, k]
            
            
            if(i == 1 & j == 1 & k == 1){
              mse <- 0
              # Conduct K-fold cross-validation.
              for(l in 1:cv){
                y_train <- y_train_k[[l]]
                y_test <- y_test_k[[l]]
                x_train <- x_train_k[[l]]
                x_test <- x_test_k[[l]]
                
                fun_reg_obj <- f_agL(y=y_train, design=x_train, Ps=Ps, active_dim=active_dim, basis_num=basis_num, weights=ww,
                                     lambda_1=lambdas_1[i], lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=cov_mat, initN=NULL)
                
                # When cv_average is set to TRUE, the estimated function is averaged across all folds.
                if(cv_average==TRUE){
                  if(l == 1){
                    best_result <- fun_reg_obj/cv
                    cv_idx <- fun_reg_obj == 0
                  }else{
                    best_result <- best_result + fun_reg_obj/cv
                    cv_idx <- cv_idx | (fun_reg_obj == 0)
                    best_result[cv_idx] <- 0
                  }
                }
                
                # Calculate test MSE.
                y_hat <- x_test %*% fun_reg_obj
                mse <- mse + sum((y_test - y_hat)^2) /length(y_test)
              }
              
              
              best_mse <- mse
              # Assign hyperparameters that minimize test MSE.
              best_lambda_1 <- lambdas_1[i]
              best_lambda_2 <- lambdas_2_vec
            }else{
              
              mse <- 0
              # Conduct K-fold cross-validation.
              for(l in 1:cv){
                y_train <- y_train_k[[l]]
                y_test <- y_test_k[[l]]
                x_train <- x_train_k[[l]]
                x_test <- x_test_k[[l]]
                
                fun_reg_obj <- f_agL(y=y_train, design=x_train, Ps=Ps, active_dim=active_dim, basis_num=basis_num, weights=ww,
                                     lambda_1=lambdas_1[i], lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=cov_mat, initN=NULL)
                # When cv_average is set to TRUE, the estimated function is averaged across all folds.
                if(cv_average==TRUE){
                  if(l == 1){
                    temp_result <- fun_reg_obj/cv
                    cv_idx <- fun_reg_obj == 0
                  }else{
                    temp_result <- temp_result + fun_reg_obj/cv
                    cv_idx <- cv_idx | (fun_reg_obj == 0)
                    temp_result[cv_idx] <- 0
                  }
                }
                y_hat <- x_test %*% fun_reg_obj
                mse <- mse + sum((y_test - y_hat)^2)/length(y_test)
              }
              # Check if the current MSE is less than the smallest MSE recorded previously.
              if(mse < best_mse){
                # Assign hyperparameters that minimize test MSE.
                best_lambda_1 <- lambdas_1[i]
                best_lambda_2[j] <- lambdas_2[j, k]
                if(cv_average==TRUE){
                  best_result <- temp_result
                }
                
                best_mse <- mse
              }
            }
            # Display the progress bar.
            cat("percent:", paste0(round(percent_complete, 1), "%,"),
                paste0("lambda_1: ", sprintf("%.4e", lambdas_1[i]), ","),
                paste0("j-th function: ", j, ","),
                paste0("lambda_2: ", sprintf("%.4e", lambdas_2[j, k]), ","),
                "            \r")
            utils::flush.console()
          }
        }
      }
      # Conduct the additional iteration for weights of roughness penalties.
      if(is.null(add_iter) == FALSE){
        # Display the progress bar.
        cat("\nAdditional Iteration \n")
        total_iterations <- add_iter * p * n_lamb_2
        iteration_counter <- 0
        end_1 <- 0
        for(iter in 1:add_iter){
          # If there is no change to the test MSE, the additional iteration is terminated early.
          if(iter > 1 & end_1 == 0){
            cat("\nNo more additional iterations are needed.")
            break
          }
          if(iter > 1 & end_1 == 1) end_1 <- 0
          # Search weights of roughness penalties, fixing the weight of the group lasso penalty that minimizes the test MSE.
          for(j in 1:p){
            lambdas_2_vec <- best_lambda_2
            for(k in 1:n_lamb_2){
              iteration_counter <- iteration_counter + 1
              percent_complete <- (iteration_counter / total_iterations) * 100
              lambdas_2_vec[j] <- lambdas_2[j, k]
              mse <- 0
              # Conduct K-fold cross-validation.
              for(l in 1:cv){
                y_train <- y_train_k[[l]]
                y_test <- y_test_k[[l]]
                x_train <- x_train_k[[l]]
                x_test <- x_test_k[[l]]
                
                fun_reg_obj <- f_agL(y=y_train, design=x_train, Ps=Ps, active_dim=active_dim, basis_num=basis_num, weights=ww,
                                     lambda_1=best_lambda_1, lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=cov_mat, initN=NULL)
                
                # When cv_average is set to TRUE, the estimated function is averaged across all folds.
                if(cv_average==TRUE){
                  if(l == 1){
                    temp_result <- fun_reg_obj/cv
                    cv_idx <- fun_reg_obj == 0
                  }else{
                    temp_result <- temp_result + fun_reg_obj/cv
                    cv_idx <- cv_idx | (fun_reg_obj == 0)
                    temp_result[cv_idx] <- 0
                  }
                }
                
                # Calculate the test MSE.
                y_hat <- x_test %*% fun_reg_obj
                mse <- mse + sum((y_test - y_hat)^2)/length(y_test)
              }
              # Check if the current MSE is less than the smallest MSE recorded previously.
              if(mse < best_mse){
                # Assign hyperparameters that minimize test MSE.
                best_lambda_2[j] <- lambdas_2[j, k]
                
                if(cv_average==TRUE){
                  best_result <- temp_result
                }
                
                best_mse <- mse
                end_1 <- 1
              }
              # Disaaply the progress bar.
              cat("percent:", paste0(round(percent_complete, 1), "%,"),
                  paste0("Add.iteration: ", iter, ","),
                  paste0("j-th function: ", j, ","),
                  paste0("lambda_2: ", sprintf("%.4e", lambdas_2[j, k]), ","),
                  "            \r")
              utils::flush.console()
            }
          }
        }
      }
    }
    
    # Using the L2 standard norm, conduct the functional linear regression.
    if(standardnorm == TRUE){
      # Search the weight of the group Lasso penalty.
      for(i in 1:n_lamb_1){
        # Search weights of roughness penalties, fixing the weight of the group lasso penalty.
        for(j in 1:n_lamb_2){
          # It is just for the progress bar.
          iteration_counter <- iteration_counter + 1
          percent_complete <- (iteration_counter / total_iterations) * 100
          
          lambdas_2_vec <- rep(lambdas_2[j], p)
          if(i == 1 & j == 1){
            mse <- 0
            # Conduct K-fold cross-validation.
            for(l in 1:cv){
              y_train <- y_train_k[[l]]
              y_test <- y_test_k[[l]]
              x_train <- x_train_k[[l]]
              x_test <- x_test_k[[l]]
              fun_reg_obj <- f_agL(y=y_train, design=x_train, Ps=Ps, active_dim=active_dim, basis_num=basis_num, weights=ww,
                                   lambda_1=lambdas_1[i], lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=cov_mat, standardnorm=TRUE, initN=NULL)
              
              # When cv_average is set to TRUE, the estimated function is averaged across all folds.
              if(cv_average==TRUE){
                if(l == 1){
                  best_result <- fun_reg_obj/cv
                  cv_idx <- fun_reg_obj == 0
                }else{
                  best_result <- best_result + fun_reg_obj/cv
                  cv_idx <- cv_idx | (fun_reg_obj == 0)
                  best_result[cv_idx] <- 0
                }
              }
              # Calculate the test MSE.
              y_hat <- x_test %*% fun_reg_obj
              mse <- mse + sum((y_test - y_hat)^2) /length(y_test)
            }
            
            
            
            best_mse <- mse
            best_lambda_1 <- lambdas_1[i]
            best_lambda_2 <- lambdas_2[j]
          }else{
            mse <- 0
            # Conduct K-fold cross-validation.
            for(l in 1:cv){
              y_train <- y_train_k[[l]]
              y_test <- y_test_k[[l]]
              x_train <- x_train_k[[l]]
              x_test <- x_test_k[[l]]
              
              fun_reg_obj <- f_agL(y=y_train, design=x_train, Ps=Ps, active_dim=active_dim, basis_num=basis_num, weights=ww,
                                   lambda_1=lambdas_1[i], lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=cov_mat, standardnorm=TRUE, initN=NULL)
              
              # When cv_average is set to TRUE, the estimated function is averaged across all folds.
              if(cv_average==TRUE){
                if(l == 1){
                  temp_result <- fun_reg_obj/cv
                  cv_idx <- fun_reg_obj == 0
                }else{
                  temp_result <- temp_result + fun_reg_obj/cv
                  cv_idx <- cv_idx | (fun_reg_obj == 0)
                  temp_result[cv_idx] <- 0
                }
              }
              # Calculate the test MSE.
              y_hat <- x_test %*% fun_reg_obj
              mse <- mse + sum((y_test - y_hat)^2)/length(y_test)
            }
            
            # Check if the current MSE is less than the smallest MSE recorded previously.
            if(mse < best_mse){
              # Assign hyperparameters that minimize the test MSE.
              best_lambda_1 <- lambdas_1[i]
              best_lambda_2 <- lambdas_2[j]
              
              if(cv_average==TRUE){
                best_result <- temp_result
              }
              
              
              best_mse <- mse
            }
          }
          # Display the progress bar.
          cat("percent:", paste0(round(percent_complete, 1), "%,"),
              paste0("lambda_1: ", sprintf("%.4e", lambdas_1[i]), ","),
              paste0("lambda_2: ", sprintf("%.4e", lambdas_2[j]), ","),
              "            \r")
          utils::flush.console()
        }
      }
    }
    
    # If cv_average is set to FALSE, linear functional regression is conducted on the entire data using hyperparameters selected through K-fold cv.
    if(cv_average==FALSE){
      if(standardnorm == TRUE){
        best_result <- f_agL(y=y-mean(y), design=design_centered, centered=TRUE, Ps=Ps, active_dim=active_dim, basis_num=basis_num, weights=ww,
                             lambda_1=best_lambda_1, lambda_2=rep(best_lambda_2, p), max_iter=max_iter, cov_mat=cov_mat, standardnorm=TRUE, initN=NULL)
      }else{
        best_result <- f_agL(y=y-mean(y), design=design_centered, centered=TRUE, Ps=Ps, active_dim=active_dim, basis_num=basis_num, weights=ww,
                             lambda_1=best_lambda_1, lambda_2=best_lambda_2, max_iter=max_iter, cov_mat=cov_mat, standardnorm=FALSE, initN=NULL)
      }
    }
    
    out <- list(beta_hat=best_result, lambda_1=best_lambda_1, lambda_2=best_lambda_2, cov_mat=cov_mat, best_mse=best_mse)
    return(out)
  }
  # If fpca==TRUE, then we use eigenfunctions as basis functions.
  if(fpca == TRUE){
    ww <- weights
    p <- length(basis_num)
    start_vec_0 <- cumsum(c(1, basis_num))[1:p]
    end_vec_0 <- cumsum(basis_num)
    start_vec <- cumsum(c(1, fpca_basis_num))[1:p]
    end_vec <- cumsum(fpca_basis_num)
    coefs <- list()
    GB <- list()
    
    # For the FPCA, calculate 'coefs' and 'GB' which are arguments of fpca_cpp.
    for(j in 1:p){
      coefs[[j]] <- t(design[, (start_vec_0[j]):(end_vec_0[j])])
      basis <- basis_list[[j]]
      if(basis$type == 'fourier') GB[[j]] <- fda::fourierpen(basis, 0)
      if(basis$type == 'bspline') GB[[j]] <- fda::bsplinepen(basis, 0)
    }
    
    # Conduct the FPCA.
    fpca_obj <- fpca_cpp(coefs=coefs, GB=GB)
    
    
    # Assign coordinate representations of covariance operators based on their eigenfunctions.
    cov_mat <- fpca_obj$Sigma
    
    
    fpca_cov_mat <- list()
    Ps_fpca <- matrix(0, sum(fpca_basis_num), sum(fpca_basis_num))
    design_fpca <- matrix(0, n, sum(fpca_basis_num))
    design_fpca_noncentered <- matrix(0, n, sum(fpca_basis_num))
    B <- matrix(0, sum(basis_num), sum(basis_num))
    for(j in 1:p){
      start_0 <- start_vec_0[j]; end_0 <- end_vec_0[j]
      start <- start_vec[j]; end <- end_vec[j]
      # Coordinate representations of explanatory funtions based on their eigenfunctions.
      design_fpca[, start:end] <- fpca_obj$pred[[j]][, 1:fpca_basis_num[j]]
      # Centering
      design_fpca_noncentered[, start:end] <- fpca_obj$non_centered[[j]][, 1:fpca_basis_num[j]]
      # Coordinate representations of covariance functions based on fpca_basis_num[[j]] eigenfunctions.
      fpca_cov_mat[[j]] <- cov_mat[[j]][1:fpca_basis_num[j], 1:fpca_basis_num[j]]
      # Coordinate representations of second-order roughness penalty based on fpca_basis_num[[j]] eigenfunctions.
      Ps_fpca[start:end, start:end] <- t(fpca_obj$mat[[j]][,1:fpca_basis_num[j]]) %*% Ps[start_0:end_0, start_0:end_0] %*% fpca_obj$mat[[j]][,1:fpca_basis_num[j]]
      # Gram matrix of B-spline.
      B[start_0:end_0, start_0:end_0] <- fda::bsplinepen(basis_list[[j]], 0)
    }
    x_train_k <- list()
    x_test_k <- list()
    y_train_k <- list()
    y_test_k <- list()
    
    # Divide total data into 'cv' folds.
    margin <- ceiling(n/cv)
    for(i in 1:cv){
      start <- 1 + margin*(i-1)
      end <- min(margin*i, n)
      idx <- (1:n) %in% (start:end)
      
      
      x_train_k[[i]] <- t(t(design_fpca_noncentered[!idx,]) - apply(design_fpca_noncentered[!idx,], 2, mean))
      x_test_k[[i]] <- t(t(design_fpca_noncentered[idx,]) - apply(design_fpca_noncentered[!idx,], 2, mean))
      y_train_k[[i]] <- y[!idx] - mean(y[!idx])
      y_test_k[[i]] <- y[idx] - mean(y[!idx])
    }
    
    # Assign the searching range of the group lasso penalty.
    if(is.null(lambda_1)){
      lambdas_1 <- exp(seq(log(1e-5), log(1e-1), length=n_lamb_1))
    }else{
      lambdas_1 <- exp(seq(log(lambda_1[1]), log(lambda_1[2]), length=n_lamb_1))
    }
    
    
    # Assign searching ranges of roughness penalties.
    # If the standard L2 norm is used for the group lasso penalty, all weights of roughness penalties are the same.
    if(standardnorm == FALSE){
      if(is.null(lambda_2)){
        lambdas_2 <- matrix(exp(seq(log(1e-9), log(1e-5), length=n_lamb_2)), p, n_lamb_2, byrow=T)
      }else{
        lambdas_2 <- matrix(0, p, n_lamb_2)
        for(j in 1:p){
          lambdas_2[j,] <- exp(seq(log(lambda_2[j, 1]), log(lambda_2[j, 2]), length=n_lamb_2))
        }
      }
    }else{
      if(is.null(lambda_2)){
        lambdas_2 <- exp(seq(log(1e-9), log(1e-5), length=n_lamb_2))
      }else{
        lambdas_2 <- exp(seq(log(lambda_2[1]), log(lambda_2[2]), length=n_lamb_2))
      }
    }
    
    
    # It is just for the progreess bar.
    if(standardnorm == FALSE){
      total_iterations <- n_lamb_1 * p * n_lamb_2
    }else{
      total_iterations <- n_lamb_1 * n_lamb_2
    }
    
    iteration_counter <- 0
    
    
    
    
    # Use the our FC-distance as the group lasso penalty.
    if(standardnorm == FALSE){
      # Search the weight of the group Lasso penalty.
      for(i in 1:n_lamb_1){
        lambdas_2_vec <- apply(lambdas_2, 1, stats::median)
        # Search weights of roughness penalties, fixing the weight of the group lasso penalty.
        for(j in 1:p){
          if(j > 1){
            lambdas_2_vec <- best_lambda_2
          }
          # Search j-th weight of roughness penalties, fixing lambda_{1i}, lambda_{21}, ... , lambda_{2, (j-1)}.
          for(k in 1:n_lamb_2){
            iteration_counter <- iteration_counter + 1
            percent_complete <- (iteration_counter / total_iterations) * 100
            lambdas_2_vec[j] <- lambdas_2[j, k]
            if(i == 1 & j == 1 & k == 1){
              mse <- 0
              
              # Conduct K-fold cross-validation.
              for(l in 1:cv){
                y_train <- y_train_k[[l]]
                y_test <- y_test_k[[l]]
                x_train <- x_train_k[[l]]
                x_test <- x_test_k[[l]]
                fun_reg_obj <- f_agL(y=y_train, design=x_train, centered=TRUE, Ps=Ps_fpca, active_dim=active_dim, basis_num=fpca_basis_num, weights=ww,
                                     lambda_1=lambdas_1[i], lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=fpca_cov_mat, initN=NULL)
                
                
                
                
                # When cv_average is set to TRUE, the estimated function is averaged across all folds.
                if(cv_average==TRUE){
                  if(l == 1){
                    best_result <- fun_reg_obj/cv
                    cv_idx <- fun_reg_obj == 0
                  }else{
                    best_result <- best_result + fun_reg_obj/cv
                    cv_idx <- cv_idx | (fun_reg_obj == 0)
                    best_result[cv_idx] <- 0
                  }
                }
                
                
                # Calculate the test MSE.
                y_hat <- x_test %*% fun_reg_obj
                mse <- mse + sum((y_test - y_hat)^2) /length(y_test)
              }
              
              best_mse <- mse
              best_lambda_1 <- lambdas_1[i]
              best_lambda_2 <- lambdas_2_vec
            }else{
              mse <- 0
              # Conduct K-fold cross-validation.
              for(l in 1:cv){
                y_train <- y_train_k[[l]]
                y_test <- y_test_k[[l]]
                x_train <- x_train_k[[l]]
                x_test <- x_test_k[[l]]
                
                fun_reg_obj <- f_agL(y=y_train, design=x_train, centered=TRUE, Ps=Ps_fpca, active_dim=active_dim, basis_num=fpca_basis_num, weights=ww,
                                     lambda_1=lambdas_1[i], lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=fpca_cov_mat, initN=NULL)
                
                
                
                
                # When cv_average is set to TRUE, the estimated function is averaged across all folds.
                if(cv_average==TRUE){
                  if(l == 1){
                    
                    temp_result <- fun_reg_obj/cv
                    cv_idx <- fun_reg_obj == 0
                  }else{
                    
                    temp_result <- temp_result + fun_reg_obj/cv
                    cv_idx <- cv_idx | (fun_reg_obj == 0)
                    temp_result[cv_idx] <- 0
                  }
                }
                
                # Calculate the test MSE.
                y_hat <- x_test %*% fun_reg_obj
                mse <- mse + sum((y_test - y_hat)^2) /length(y_test)
              }
              
              # Check if the current MSE is less than the smallest MSE recorded previously.
              if(mse < best_mse){
                # Assign hyperparameters that minimize test MSE.
                best_lambda_1 <- lambdas_1[i]
                best_lambda_2[j] <- lambdas_2[j, k]
                if(cv_average==TRUE){
                  best_result <- temp_result
                }
                
                best_mse <- mse
              }
            }
            
            # Display the progress bar.
            cat("percent:", paste0(round(percent_complete, 1), "%,"),
                paste0("lambda_1: ", sprintf("%.4e", lambdas_1[i]), ","),
                paste0("j-th function: ", j, ","),
                paste0("lambda_2: ", sprintf("%.4e", lambdas_2[j, k]), ","),
                "            \r")
            utils::flush.console()
          }
        }
      }
      
      # Conduct the additional iteration for weights of roughness penalties.
      if(is.null(add_iter) == FALSE){
        # Display the progress bar.
        cat("\nAdditional Iteration \n")
        total_iterations <- add_iter * p * n_lamb_2
        iteration_counter <- 0
        end_1 <- 0
        
        # Conduct the additional iteration for weights of roughness penalties.
        for(iter in 1:add_iter){
          
          # If there is no change to the test MSE, the additional iteration is terminated early.
          if(iter > 1 & end_1 == 0){
            cat("\nNo more additional iterations are needed.")
            break
          }
          
          
          if(iter > 1 & end_1 == 1) end_1 <- 0
          
          # Search weights of roughness penalties, fixing the weight of the group lasso penalty that minimizes the test MSE.
          for(j in 1:p){
            lambdas_2_vec <- best_lambda_2
            for(k in 1:n_lamb_2){
              iteration_counter <- iteration_counter + 1
              percent_complete <- (iteration_counter / total_iterations) * 100
              lambdas_2_vec[j] <- lambdas_2[j, k]
              mse <- 0
              # Conduct the K-fold cross-validation.
              for(l in 1:cv){
                y_train <- y_train_k[[l]]
                y_test <- y_test_k[[l]]
                x_train <- x_train_k[[l]]
                x_test <- x_test_k[[l]]
                
                fun_reg_obj <- f_agL(y=y_train, design=x_train, centered=TRUE, Ps=Ps_fpca, active_dim=active_dim, basis_num=fpca_basis_num, weights=ww,
                                     lambda_1=best_lambda_1, lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=fpca_cov_mat, initN=NULL)
                
                
                # When cv_average is set to TRUE, the estimated function is averaged across all folds.
                if(cv_average==TRUE){
                  if(l == 1){
                    
                    temp_result <- fun_reg_obj/cv
                    cv_idx <- fun_reg_obj == 0
                  }else{
                    
                    temp_result <- temp_result + fun_reg_obj/cv
                    cv_idx <- cv_idx | (fun_reg_obj == 0)
                    temp_result[cv_idx] <- 0
                  }
                }
                
                
                # Calculate the test MSE.
                y_hat <- x_test %*% fun_reg_obj
                mse <- mse + sum((y_test - y_hat)^2) /length(y_test)
              }
              
              # Check if the current MSE is less than the smallest MSE recorded previously.
              if(mse < best_mse){
                # Assign hyperparameters that minimize test MSE.
                best_lambda_2[j] <- lambdas_2[j, k]
                if(cv_average==TRUE){
                  best_result <- temp_result
                }
                
                
                best_mse <- mse
                end_1 <- 1
              }
              # Display the progress bar.
              cat("percent:", paste0(round(percent_complete, 1), "%,"),
                  paste0("Add.iteration: ", iter, ","),
                  paste0("j-th function: ", j, ","),
                  paste0("lambda_2: ", sprintf("%.4e", lambdas_2[j, k]), ","),
                  "            \r")
              utils::flush.console()
            }
          }
        }
      }
    }
    
    
    # Using the L2 standard norm, conduct the functional linear regression.
    if(standardnorm == TRUE){
      # Search the weight of the group Lasso penalty.
      for(i in 1:n_lamb_1){
        # Search weights of roughness penalties, fixing the weight of the group lasso penalty.
        for(j in 1:n_lamb_2){
          iteration_counter <- iteration_counter + 1
          percent_complete <- (iteration_counter / total_iterations) * 100
          lambdas_2_vec <- rep(lambdas_2[j], p)
          if(i == 1 & j == 1){
            mse <- 0
            # K-fold cross-validation.
            for(l in 1:cv){
              y_train <- y_train_k[[l]]
              y_test <- y_test_k[[l]]
              x_train <- x_train_k[[l]]
              x_test <- x_test_k[[l]]
              fun_reg_obj <- f_agL(y=y_train, design=x_train, centered=TRUE, Ps=Ps_fpca, active_dim=active_dim, basis_num=fpca_basis_num, weights=ww,
                                   lambda_1=lambdas_1[i], lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=fpca_cov_mat, standardnorm=TRUE, initN=NULL)
              
              # When cv_average is set to TRUE, the estimated function is averaged across all folds.
              if(cv_average==TRUE){
                if(l == 1){
                  best_result <- fun_reg_obj/cv
                  cv_idx <- fun_reg_obj == 0
                }else{
                  best_result <- best_result + fun_reg_obj/cv
                  cv_idx <- cv_idx | (fun_reg_obj == 0)
                  best_result[cv_idx] <- 0
                }
              }
              # Test MSE.
              y_hat <- x_test %*% fun_reg_obj
              mse <- mse + sum((y_test - y_hat)^2) /length(y_test)
            }
            
            best_mse <- mse
            best_lambda_1 <- lambdas_1[i]
            best_lambda_2 <- lambdas_2[j]
          }else{
            mse <- 0
            
            # K-fold cross-validation.
            for(l in 1:cv){
              y_train <- y_train_k[[l]]
              y_test <- y_test_k[[l]]
              x_train <- x_train_k[[l]]
              x_test <- x_test_k[[l]]
              
              fun_reg_obj <- f_agL(y=y_train, design=x_train, centered=TRUE, Ps=Ps_fpca, active_dim=active_dim, basis_num=fpca_basis_num, weights=ww,
                                   lambda_1=lambdas_1[i], lambda_2=lambdas_2_vec, max_iter=max_iter, cov_mat=fpca_cov_mat, standardnorm=TRUE, initN=NULL)
              
              # When cv_average is set to TRUE, the estimated function is averaged across all folds.
              if(cv_average==TRUE){
                if(l == 1){
                  temp_result <- fun_reg_obj/cv
                  cv_idx <- fun_reg_obj == 0
                }else{
                  temp_result <- temp_result + fun_reg_obj/cv
                  cv_idx <- cv_idx | (fun_reg_obj == 0)
                  temp_result[cv_idx] <- 0
                }
              }
              # Test MSE.
              y_hat <- x_test %*% fun_reg_obj
              mse <- mse + sum((y_test - y_hat)^2)/length(y_test)
            }
            
            # Check if the current MSE is less than the smallest MSE recorded previously.
            if(mse < best_mse){
              # Assign hyperparameters that minimize test MSE.
              best_lambda_1 <- lambdas_1[i]
              best_lambda_2 <- lambdas_2[j]
              if(cv_average==TRUE){
                best_result <- temp_result
              }
              
              best_mse <- mse
            }
          }
          # Progress bar.
          cat("percent:", paste0(round(percent_complete, 1), "%,"),
              paste0("lambda_1: ", sprintf("%.4e", lambdas_1[i]), ","),
              paste0("lambda_2: ", sprintf("%.4e", lambdas_2[j]), ","),
              "            \r")
          utils::flush.console()
        }
      }
    }
    
    # If cv_average is set to FALSE, linear functional regression is conducted on the entire data using hyperparameters selected through K-fold cv.
    if(cv_average==FALSE){
      if(standardnorm == TRUE){
        best_result <- f_agL(y=y-mean(y), design=design_fpca, centered=TRUE, Ps=Ps_fpca, active_dim=active_dim, basis_num=fpca_basis_num, weights=ww,
                             lambda_1=best_lambda_1, lambda_2=rep(best_lambda_2, p), max_iter=max_iter, cov_mat=fpca_cov_mat, standardnorm=TRUE, initN=NULL)
      }else{
        best_result <- f_agL(y=y-mean(y), design=design_fpca, centered=TRUE, Ps=Ps_fpca, active_dim=active_dim, basis_num=fpca_basis_num, weights=ww,
                             lambda_1=best_lambda_1, lambda_2=best_lambda_2, max_iter=max_iter, cov_mat=fpca_cov_mat, standardnorm=FALSE, initN=NULL)
      }
    }
    
    # Transform the basis system from eigenfunctions to B-spline or Fourier basis functions.
    best_result_0 <- rep(0, sum(basis_num))
    for(j in 1:p){
      start_0 <- start_vec_0[j]
      end_0 <- end_vec_0[j]
      start <- start_vec[j]
      end <- end_vec[j]
      best_result_0[start_0:end_0] <- fpca_obj$mat[[j]][,1:fpca_basis_num[j]] %*% best_result[start:end]
      cov_mat[[j]] <- fpca_obj$mat[[j]] %*% cov_mat[[j]] %*% t(fpca_obj$mat[[j]]) %*% fpca_obj$GB[[j]]
    }
    out <- list(beta_hat=best_result_0, lambda_1=best_lambda_1, lambda_2=best_lambda_2, cov_mat=cov_mat, best_mse=best_mse)
    return(out)
  }
}
