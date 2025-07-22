# FAGLasso
This R-package is to conduct the scalar-on-function linear regeression with variable selection. 

The model is given by

$$
Y = \alpha + \sum_{j=1}^p \langle \beta_j, X_j \rangle_{H} + \epsilon,
$$

where $H$ is a separable Hilbert space of functions.

To perform variable selection, we penalize the variance of each projection $\langle \beta_j, X_j \rangle_{H}$ as follows:

$$
\hat \alpha, \hat \beta_1, \cdots, \hat\beta_p = \text{argmin}_{\alpha \in \mathbb R, \beta_j \in H} E_n \left\| Y - \alpha  - \sum _{j=1}^p \langle \beta_j, X_j \rangle _{H} \right\|^2 + \lambda \sum _{j=1}^p w _{1j} \left( \text{var} \left( \langle \beta_j, X_j \rangle _{H} \right) \right)^{1/2}.
$$


To prevent overfitting, we add a smoothness penalty as:

$$
\hat \alpha, \hat \beta_1, \cdots, \hat\beta_p = \text{argmin}_{\alpha \in \mathbb R, \beta_j \in H} E_n \left\| Y - \alpha  - \sum _{j=1}^p \langle \beta_j, X_j \rangle _{H} \right\|^2 + \lambda \sum _{j=1}^p w _{1j} \left( \text{var} \left( \langle \beta_j, X_j \rangle _{H} \right) + w _{2j} \| \beta_j ^{\prime\prime}\|^2 \right)^{1/2}.
$$



The package supports two modes: Fourier basis and B-spline basis. Example code for each mode is provided in the Examples section.




## Installiation

You can install the package using the following R code:

```r
library(devtools)
devtools::install_github("JSongLab/FAGLASSO", force=TRUE)
```


## Example: Fourier basis
The first example is the code when Fourier basis functions are used.

```r
############################################################################
##### This is the example code for functional adaptive group lasso. ########
############################################################################

library(fda)
library(FAGLasso)
library(fields)
n <- 100; p <- 5; s <- 2; m <- 15; basis_num <- rep(m, p)
init <- NULL
active_dim <- 1e-6
cv <- 5
max_iter <- 1000
#### Generate population data.
#### X_1, ... X_5 ~ i.i.d. zero-mean Materen process(sigma=1, rho=1, nu=0.5) whose domain is (0,1).
#### Y_i = <X_1i, sin(2 pi t)> + <X_2i, cos(2 pi t)> + e, where e ~ N(0, 0.1^2) 
#### Thus X_1, X_2 are active functions and X_3, X_4, X_5 are inactive functions.
nt0 <- 1001
X <- array(0, c(n, nt0, s))
times <- seq(0, 1, length=nt0)
cov_mat <- apply(abs(outer(times, times, "-")), c(1,2), FUN=Matern, range=1, nu=0.5)
cov_mat_sq <- matpower_1(cov_mat, 1/2)
for(i in 1:s){
  z <- matrix(rnorm(n*nt0, 0, 1), n, nt0)
  X[,,i] <- z %*% cov_mat_sq 
}

times <- seq(0, 1, length=nt0)
y <- rep(0, n)
beta_1 <- sin(2*pi*times); beta_2 <- cos(2*pi*times)
temp1 <- rep(0, n); temp2 <- rep(0, n)
for(i in 1:n){
  temp1[i] <- numeric_integ(beta_1*X[i,,1], times)
  temp2[i] <- numeric_integ(beta_2*X[i,,2], times)
  y[i] <- 1 + temp1[i] + temp2[i]
}

sigma <- 0.1
y <- y + rnorm(n, 0, sigma)


#### In real analysis, expnlanatory functions are not fully observed. 
#### Hence, we asuume that X_1 ~ X_5 are observed evenly at 21 points within the interval (0,1).
#### And then, 15 Fourier basis functions are used to represent the X_1 ~ X5 in their coordinate form.
t_idx <- seq(1, 1001, by=50)
times <- times[t_idx]
nt <- length(t_idx)
x <- array(0, c(n, nt, p))

for(i in 1:p){
  if(i <= s) x[,,i] <- X[,,i][,t_idx]
  else{
    z <- matrix(rnorm(n*nt, 0, 1), n, nt)
    x[,,i] <- z %*% cov_mat_sq[t_idx, t_idx]
  } 
}
par(mfrow=c(1,1))
idx <- sample(1:100, 1)
plot(times, x[idx,,4], type="l")

basis_list <- list()
for(j in 1:p){
  basis_list[[j]] <- create.fourier.basis(rangeval=c(0,1), nbasis=basis_num[j])
}
lambdas <- seq(1e-4, 1e-3, length=100)
gcvs <- rep(0, 100)
for(i in 1:n){
  temp <- smooth.basis(times, t(x[,,1]), fdPar(basis_list[[1]], 2, lambda=lambdas[i]))
  gcvs[i] <- mean(temp$gcv)
}
best_lambda <- lambdas[which.min(gcvs)]
plot(gcvs, type="l") ### If this plot don't show U-shape, change the range of lambdas


##### Now, contruct the arguments of Fun_Reg_AGL
coefs <- list()
Ps <- matrix(0, sum(basis_num), sum(basis_num))
design <- matrix(0, n, sum(basis_num))
start_vec <- cumsum(c(1, basis_num))[1:p]
end_vec <- cumsum(basis_num)
for(j in 1:p){
  coefs[[j]] <- smooth.basis(times, t(x[,,j]), fdPar(basis_list[[j]], 2, lambda=best_lambda))$fd$coefs 
  design[, (start_vec[j]):(end_vec[j])] <- t(coefs[[j]])
  Ps[start_vec[j]:end_vec[j], start_vec[j]:end_vec[j]] <- fourierpen(create.fourier.basis(c(0,1), nbasis=basis_num[j]), 2)
}

ww <- rep(1, p) #### It means non-adaptive group lasso.
lambda_2 <- matrix(rep(c(1e-4, 1e-1), p), ncol=2, byrow=T)

#### Conduct non-adapive group lasso.
fun_reg_obj <- Fun_Reg_AGL(y=y, design=design, Ps=Ps, basis_num=basis_num, weights=ww,
                           lambda_1=c(1e-3, 1e+1),
                           lambda_2=lambda_2, 
                           n_lamb_1=20,
                           n_lamb_2=10,
                           cv=5, 
                           active_dim=1e-6,
                           add_iter=3,
                           max_iter=2000)

#### Check the selected lambdas are near or exact at the bound of searcing range.
fun_reg_obj$lambda_1
fun_reg_obj$lambda_2

#### Construct weights_j = 1/sd(<X_j, beta_j>)
#### Note that the output 'cov_mat' is used calculating the weights.
for(j in 1:p){
  start <- m*(j-1)+1
  end <- m*j
  beta_j <- fun_reg_obj$beta_hat[start:end]
  ww[j] <- 1/sqrt(t(beta_j) %*% fun_reg_obj$cov_mat[[j]] %*% beta_j)
}

idx2 <- ww != Inf
idx <- rep(idx2, basis_num)
idx3 <- (1:p)[idx2]


#### Conduct the adaptive group lasso.
fun_reg_obj <- Fun_Reg_AGL(y=y, design=design[,idx], Ps=Ps[idx, idx], basis_num=basis_num[idx2], weights=ww[idx2],
                           lambda_1=c(1e-3, 1e+1),
                           lambda_2=lambda_2[idx2,], 
                           n_lamb_1=20,
                           n_lamb_2=10,
                           cv=5, 
                           active_dim=1e-6,
                           add_iter=3,
                           max_iter=1000)

#### Check the selected lambdas are near or exact at the bound of searcing range.
fun_reg_obj$lambda_1
fun_reg_obj$lambda_2



#### Performance part ####
selection <- rep(0, length(ww[idx2]))
for(i in 1:length(selection)){
  start <- start_vec[i]
  end <- end_vec[i]
  beta_j <- fun_reg_obj$beta_hat[start:end]
  selection[i] <- 1/(sqrt(sum(beta_j^2)))
}


idx3 <- idx3[selection != Inf]
sel_1 <- s - sum(idx3 %in% (1:s))
sel_2 <- sum(idx3 %in% ((s+1):p))

ttt <- seq(0, 1, by=0.01)
evals_1 <- eval.fd(ttt, fd(fun_reg_obj$beta_hat[start_vec[1]:end_vec[1]], basis_list[[1]]))
evals_2 <- eval.fd(ttt, fd(fun_reg_obj$beta_hat[start_vec[2]:end_vec[2]], basis_list[[2]]))


#### Estimation performance 
#### Black lines : estimated functions
#### Red lines: true functions
par(mfrow=c(1,2))
plot(ttt, evals_1, type="l", ylim=c(-1.1, 1.1))
lines(ttt, sin(2*pi*ttt), col='red')
plot(ttt, evals_2, type="l", ylim=c(-1.1, 1.1))
lines(ttt, cos(2*pi*ttt), col='red')

#### Selection performance 
print(paste("The numbdf of false exclusions is ", sel_1))
print(paste("The numbdf of false inclusions is ", sel_2))


```


## Example: B-spline basis

The seconde example is the code when B-spline basis functions are used.

```r
############################################################################
##### This is the example code for functional adaptive group lasso. ########
############################################################################

library(FAGLasso)
library(fda)
library(fields)
n <- 100; p <- 5; s <- 2; m <- 15; basis_num <- rep(m, p)


#### Generate population data.
#### X_1, ... X_5 ~ i.i.d. zero-mean Materen process(sigma=1, rho=1, nu=0.5) whose domain is (0,1).
#### Y_i = <X_1i, t> + <X_2i, 10(t-0.5)^2> + e, where e ~ N(0, 0.1^2) 
#### Thus X_1, X_2 are active functions and X_3, X_4, X_5 are inactive functions.
nt0 <- 1001
X <- array(0, c(n, nt0, s))
times <- seq(0, 1, length=nt0)
cov_mat <- apply(abs(outer(times, times, "-")), c(1,2), FUN=Matern, range=1, nu=0.5)
cov_mat_sq <- matpower_1(cov_mat, 1/2)

for(i in 1:s){
  z <- matrix(rnorm(n*nt0, 0, 1), n, nt0)
  X[,,i] <- z %*% cov_mat_sq 
}

times <- seq(0, 1, length=nt0)
y <- rep(0, n)
beta_1 <- times; beta_2 <- 10*(times - 0.5)^2
temp1 <- rep(0, n); temp2 <- rep(0, n)
for(i in 1:n){
  temp1[i] <- numeric_integ(beta_1*X[i,,1], times)
  temp2[i] <- numeric_integ(beta_2*X[i,,2], times)
  y[i] <- 1 + temp1[i] + temp2[i]
}
sigma <- 0.1
y <- y + rnorm(n, 0, sigma)

#### In real analysis, expnlanatory functions are not fully observed. 
#### Hence, we asuume that X_1 ~ X_5 are observed evenly at 21 points within the interval (0,1).
#### And then, 15 B-spline basis functions are used to represent the X_1 ~ X5 in their coordinate form.
t_idx <- seq(1, 1001, by=50)
times <- times[t_idx]

nt <- length(t_idx)
x <- array(0, c(n, nt, p))

for(i in 1:p){
  if(i <= s) x[,,i] <- X[,,i][,t_idx]
  else{
    z <- matrix(rnorm(n*nt, 0, 1), n, nt)
    x[,,i] <- z %*% cov_mat_sq[t_idx, t_idx]
  } 
}
par(mfrow=c(1,1))
idx <- sample(1:100, 1)
plot(times, x[idx,,4], type="l") 

basis_list <- list()
for(j in 1:p){
  basis_list[[j]] <- create.bspline.basis(rangeval=c(0,1), nbasis=basis_num[j])
}
lambdas <- seq(1e-5, 1e-4, length=100)
gcvs <- rep(0, 100)
for(i in 1:n){
  temp <- smooth.basis(times, t(x[,,1]), fdPar(basis_list[[1]], 2, lambda=lambdas[i]))
  gcvs[i] <- mean(temp$gcv)
}
best_lambda <- lambdas[which.min(gcvs)]
plot(gcvs, type="l") ### If this plot don't show U-shape, change the range of lambdas


##### Now, contruct the arguments of Fun_Reg_AGL
coefs <- list()
Ps <- matrix(0, sum(basis_num), sum(basis_num))
design <- matrix(0, n, sum(basis_num))
start_vec <- cumsum(c(1, basis_num))[1:p]
end_vec <- cumsum(basis_num)
for(j in 1:p){
  coefs[[j]] <- smooth.basis(times, t(x[,,j]), fdPar(basis_list[[j]], 2, lambda=best_lambda))$fd$coefs 
  design[, (start_vec[j]):(end_vec[j])] <- t(coefs[[j]])
  Ps[start_vec[j]:end_vec[j], start_vec[j]:end_vec[j]] <- bsplinepen(basis_list[[j]], 2)
}

#### Conduct non-adapive group lasso.
#### Note that the argument 'fpca' of Fun_Reg_AGL is TRUE.
#### Hence, basis system is eigenfunctions obtained fpca to estimate true functions.
#### Note that this fpca basis are only used inside of Fun_Reg_AGL.
#### Thus, coordinates of estimated functions are based on original system such as B-spline or Fourier.

ww <- rep(1, p)  #### It means non-adaptive group lasso.
lambda_2 <- matrix(rep(c(1e-4, 1e-0), p), ncol=2, byrow=T)

init <- NULL
active_dim <- 1e-10
cv <- 5
max_iter <- 2000
#### Conduct non-adapive group lasso.
fun_reg_obj <- Fun_Reg_AGL(y=y, design=design, Ps=Ps,
                           basis_num=basis_num, basis_list=basis_list,
                           weights=ww,
                           fpca=TRUE,
                           fpca_basis_num=rep(10, p),
                           lambda_1=c(1e-3, 1e+1),
                           lambda_2=lambda_2, 
                           n_lamb_1=20,
                           n_lamb_2=10,
                           cv=cv, 
                           cv_average=TRUE,
                           active_dim=active_dim,
                           add_iter=10,
                           max_iter=2*max_iter)

#### Check the selected lambdas are near or exact at the bound of searcing range.
fun_reg_obj$lambda_1
fun_reg_obj$lambda_2


#### Construct weights_j = 1/sd(<X_j, beta_j>)
#### Note that the output 'cov_mat' is used calculating the weights.
for(j in 1:p){
  start <- start_vec[j]
  end <- end_vec[j]
  beta_j <- fun_reg_obj$beta_hat[start:end]
  GB <- bsplinepen(basis_list[[j]], 0)
  ww[j] <- 1/sqrt(t(beta_j) %*% GB %*% fun_reg_obj$cov_mat[[j]] %*% beta_j)
}

idx2 <- ww != Inf
idx <- rep(idx2, basis_num)
idx3 <- (1:p)[idx2]

#### Conduct the adaptive group lasso.
fun_reg_obj <- Fun_Reg_AGL(y=y, design=design[,idx], Ps=Ps[idx, idx],
                           basis_num=basis_num[idx2], basis_list=basis_list,
                           weights=ww[idx2],
                           fpca=TRUE,
                           fpca_basis_num=rep(10, p)[idx2],
                           lambda_1=c(1e-3, 1e+1),
                           lambda_2=lambda_2[idx2,], 
                           n_lamb_1=20,
                           n_lamb_2=10,
                           cv=5, 
                           cv_average=TRUE,
                           #cv_average=FALSE,
                           active_dim=1e-6,
                           add_iter=3,
                           max_iter=2000)

#### Check the selected lambdas are near or exact at the bound of searcing range.
fun_reg_obj$lambda_1
fun_reg_obj$lambda_2


#### Performance part ####
selection <- rep(0, length(ww[idx2]))
for(i in 1:length(selection)){
  start <- m*i - (m-1)
  end <- m*i
  beta_j <- fun_reg_obj$beta_hat[start:end]
  selection[i] <- 1/(sqrt(sum(beta_j^2)))
}
idx3 <- idx3[selection != Inf]
sel_1 <- s - sum(idx3 %in% (1:s))
sel_2 <- sum(idx3 %in% ((s+1):p))

ttt <- seq(0, 1, by=0.01)
evals_1 <- eval.fd(ttt, fd(fun_reg_obj$beta_hat[start_vec[1]:end_vec[1]], basis_list[[1]]))
evals_2 <- eval.fd(ttt, fd(fun_reg_obj$beta_hat[start_vec[2]:end_vec[2]], basis_list[[2]]))

#### Estimation performance 
#### Black lines : estimated functions
#### Red lines: true functions
par(mfrow=c(1,2))
plot(ttt, evals_1, type="l", ylim=c(-0.1, 1.1))
lines(ttt, ttt, col='red')
plot(ttt, evals_2, type="l", ylim=c(-0.3, 2.8))
lines(ttt, 10*(ttt-0.5)^2, col='red')

#### Selection performance 
print(paste("The numbdf of false exclusions is ", sel_1))
print(paste("The numbdf of false inclusions is ", sel_2))


```




