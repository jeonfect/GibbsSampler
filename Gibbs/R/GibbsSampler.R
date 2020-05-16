# Contets of Gibbs_Sampler.R


#' Solve method for variance matrices.
#'
#' @param V Variance matrix
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
#' @return Matrix solving system of equations and optionally the log-determinant.
#' @details This function is faster and more stable than `base::solve()` when `V` is known to be positive-definite.
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V) # cholesky decomposition
  if(missing(x)) x <- diag(nrow(V))
  # solve is O(ncol(C) * ncol(x)) with triangular matrices
  # using backward subsitution
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}

#' Sampling method for observation means.
#'
#' @param x observation `n x p`
#' @param mu matrix of cluster means `k x p`
#' @param Sigma matrix of cluster variances `p x p x n`
#' @param V matrix of variance for each observation x `p x p x n`
#' @param z vector of each observation's cluster allocation, length `n`
#' @return matrix `n x p` of parameter values for observation x
#' @details `theta_i` has conditional distribution `N(G_i(y_i-mu_{z_i})+mu_{z_i},G_iV_i)` where `G_i=Sig_{z_i}(V_i+Sig_{z_i})^{-1}`. All sampling can be done through `rRxNorm()` from `mniw` package.
#' @references Martin Lysy and Bryan Yates (2019). mniw: The Matrix-Normal Inverse-Wishart Distribution. R package version 1.0. https://CRAN.R-project.org/package=mniw
sampling_theta <- function(x, mu, Sigma, V, z) {
  n <- nrow(x)
  p <- ncol(x)
  mu_z <- matrix(0,n,p)
  Sig_z <- array(0,dim=c(p,p,n))
  for(i in 1:n) {
    mu_z[i,] <- mu[z[i],]
    Sig_z[,,i] <- Sigma[,,z[i]]
  }
  return(rRxNorm(n,x,V,mu_z,Sig_z))
}

#' Sampling method for within-group sample mean.
#'
#' @param z vector length `n` of observation's cluster allocation
#' @param theta matrix `n x p` of observation means
#' @param mu matrix `k x p` of mean values for each cluster
#' @param Sigma matrix `p x p x k` of variances for each cluster
#' @return  matrix `k x p` of sampled mu
#' @details `mu_k` has conditional distribution `N(\hat\theta_k,\Sigma_k/N_k)` where `N_k=\sum_{i=1}^{N}I(z_i=k)`. Sampling is done using `rmNorm()` from `mniw` package.
#' @references Martin Lysy and Bryan Yates (2019). mniw: The Matrix-Normal Inverse-Wishart Distribution. R package version 1.0. https://CRAN.R-project.org/package=mniw
sampling_mu <- function(z, theta, mu, Sigma) {
  n <- nrow(theta)
  p <- ncol(mu)
  k <- nrow(mu)
  theta_k <- matrix(0,k,p)
  Sig_k <- array(0,dim=c(p,p,k))
  new_mu <- matrix(0,k,p)

  for(i in 1:k) {
    ind <- which(z==i)
    if(length(ind)==0){
      theta_k[i,] <- mu[i,]
      Sig_k <- Sigma[,,i]
    } else {
      if(length(ind)==1){ #only one observation has assigned group i
        theta_k[i,] <- theta[ind,]
      } else{
        theta_k[i,] <- colSums(theta[ind,])/length(ind)
      }
      Sig_k[,,i] <- Sigma[,,i]/length(ind)
    }
  }

  return(rmNorm(k,theta_k,Sig_k))
}

#' Sampling method for within-group variance-covariance matrix
#'
#' @param z vector length `n` of observation's cluster allocation
#' @param mu matrix `n x p` of mean values for each cluster
#' @param Omega matrix `p x p x k` of between row precisions
#' @param theta matrix `n x p` of observation means
#' @param v_k vector length `k` of prior contribution to each cluster
#' @return matrix `p x p x k` of sampled sigma
#' @details `Sigma_k` has conditional distribution `InvWish(\Omega_k + \sum_{i:I(z_i=k)}(\theta_i-\mu_k)(\theta_i-\mu_k)',N_k+v_k)`. Sampling is done using `riwish()` function from `mniw` package.
#' @references Martin Lysy and Bryan Yates (2019). mniw: The Matrix-Normal Inverse-Wishart Distribution. R package version 1.0. https://CRAN.R-project.org/package=mniw
sampling_sigma <- function(z, mu, Omega, theta, v_k) {
  p <- ncol(mu)
  k <- nrow(mu)
  psi <- array(0, dim=c(p,p,k))
  df <- rep(0,k)

  for(i in 1:k) {
    ind <- which(z==i)
    temp <- theta[ind,]-mu[i,]
    psi[,,i] <- Omega[,,i]
    for(j in 1:length(ind)){
      if(length(ind)==1){
        psi[,,i] <- psi[,,i] + crossprod(t(temp))
      } else{
        psi[,,i] <- psi[,,i] + crossprod(t(temp[j,]))
      }
    }
    df[i] <- length(ind) + v_k[i]
  }
  return(riwish(k, psi, df))
}

#' Sampling method for probability vector for group membership.
#'
#' @param z vector length `n` of observation's cluster allocation
#' @param k integer of number of clusters in the data
#' @return vector length `k` of prior distribution for clusters
#' @details `rho_k` has conditional distribution `Dirichlet(\alpha)`, where `\alpha=(N_k+1)`.
sampling_rho <- function(z, k) {
  tau <- rep(0,k)
  for(i in 1:k) {
    n_k <- length(which(z==i))
    tau[i] <- rgamma(1,n_k+1,1)
  }
  tau_sum <- sum(tau)

  return(tau/tau_sum) #equivalent to Dirichlet distribution
}

#' Sampling method for probability matrix for group membership and assign group membership based on these probabilities
#'
#' @param theta matrix `n x p` of parameter values for observation x
#' @param rho vector length `k` of prior distribution for clusters
#' @param mu matrix `k x p` of mean values for each clusters
#' @param Sigma matrix `p x p x k` of variance for each clusters
#' @param prior_z vector length `n` of clusters allocation, which will be returned instead of newly sampled cluster allocation when at least one N_k = 0
#' @param prior_lambda matrix `n x k` of posterior probabilities, which will be returned instead of newly sampled posteriro probabilities when at least one N_k = 0
#' @return matrix `n x k` of observations' posterior probabilities
#' @return vector length `n` of observations' pooled cluster allocation
#' @details Group membership `z` is chosen based on maximum a posteriori criteria. If there is at least 1 group that has not been assigned, it would not update `lambda` and `z`.
sampling_z <- function(theta, rho, mu, Sigma, prior_z, prior_lambda) {
  n <- nrow(theta)
  p <- ncol(theta)
  k <- length(rho)
  log_rho <- log(rho)
  inv_sig <- array(0,dim=c(p,p,k))
  log_det <- rep(0,k)
  lambda_mat <- matrix(0,n,k)
  Z <- rep(0,n)

  #Store inverse matrices and log determinants of Sigma
  for(i in 1:k){
    inv_sig[,,i] <- solveV(Sigma[,,i])
    log_det[i] <- log(det(Sigma[,,i]))
  }

  for(i in 1:n){
    for(j in 1:k){
      Siga <- crossprod(inv_sig[,,j],(theta[i,]-mu[j,]))
      aSiga <- crossprod((theta[i,]-mu[j,]),Siga)
      lambda_mat[i,j] <- log_rho[j] - 1/2*(aSiga + log_det[j]) #log density
    }
    lambda_mat[i,] <- exp(lambda_mat[i,])
    lambda_mat[i,] <- lambda_mat[i,]/sum(lambda_mat[i,])
    Z[i] <- which.max(lambda_mat[i,]) #Maximum a Posteriori
  }

  #If there is at least one group not chosen as posterior, do not update
  if(length(unique(Z))!=k) {
    return(list(lambda=prior_lambda, z=prior_z))
  }
  res <- list(lambda=lambda_mat, z=Z)
  return(res)
}


#' Gibbs Sampler for model-based clustering
#'
#' @param x matrix `n x p` of observations
#' @param k number of clusters
#' @param m_iter integer number of iterations used for Gibbs sampler
#' @param init list of initial values of parameters(`mu`,`Sigma`,`theta`,`rho`,`z`). If not specified, this function assigns initial values based on `x`
#' @param V martrix `p x p x N` of variances corresponding to each observation
#' @param burn_in positive integer number of iterations that are only used for throwing out garbage values. If not specified, 10% of `m_iter` is used as default
#' @param theta_out logical value whether to save every `theta` produced through all iterations.. Default is `FALSE`
#' @param z_out logical value whether to save every `z` updated through all iterations. Default is `FALSE`.
#' @return a list with elements:
#' \describe{
#'   \item{`mu`}{matrix `n x p x m_iter` of mean values for clutsers}
#'   \item{`sigma`}{matrix `p x p x k x m_iter` of within-group variance}
#'   \item{`rho`}{matrix `n x k x m_iter` of proportion of group membership}
#'   \item{`lambda`}{matrix `n x k` of average probability of group membership}
#'   \item{`z`}{vector length `n` of posterior group membership if `z_out=FALSE`. matrix `n x m_iter` is produced otherwise}
#'   \item{`theta`}{matrix `n x p` of observation means if `theta_out=FALSE`, matrix `n x p x m_iter` is produced otherwise}
#'   \item{`burn`}{integer number of burn-in periods}
#'   \item{`m_iter`}{integer number of total iterations}
#' }
#' @export
gibbs <- function(x, k, m_iter, init, V, burn_in, theta_out=FALSE, z_out=FALSE) {
  n <- nrow(x)
  p <- ncol(x)
  if(missing(burn_in)) {
    burn_in <- round(0.10*m_iter) #default is to assign 10% as burn_in
  }
  if(missing(init)) { # initialize parameters
    kmean <- kmeans(x,k)
    z <- kmean$cluster
    n_k <- rep(0,k)
    init_rho <- rep(0,k)
    init_var <- array(0, dim=c(p,p,k))
    init_mu <- matrix(0,k,p)
    for(i in 1:k){
      n_k[i] <- length(which(z==i))
      init_rho[i] <- n_k[i]/n
      if(n_k[i]==0){
        init_mu[i,] <- colSums(x)/n
      } else {
        init_mu[i,] <- colSums(x[z==i,])/n_k[i]
      }
      if(n_k[i] < p) {
        init_var[,,i] <- var(x)
      } else {
        init_var[,,i] <- var(x[z==i,])
      }
    }
    init_theta <- x
  }
  else{
    z <- init$z
    init_rho <- init$rho
    init_mu <- init$mu
    init_var <- init$var
    init_theta <- init$theta
  }

  iters <- m_iter - burn_in
  v_k <- rep(p+2,k)
  lambda <- matrix(1/k,n,k)
  Omega <- array(var(x),dim=c(p,p,k))

  result <- list(mu=array(0, dim=c(k,p,m_iter)),
                 sigma=array(0, dim=c(p,p,k,m_iter)),
                 rho=array(0, dim=c(k,m_iter)),
                 lambda=matrix(0,n,k))

  if(z_out){
    result$z <- matrix(0,n,m_iter)
    result$z[,1] <- z
  } else {
    result$z <- z
  }
  if(theta_out){
    result$theta <- array(0, dim=c(n,p,m_iter))
    result$theta[,,1] <- init_theta
  }
  result$mu[,,1] <- init_mu
  result$sigma[,,,1] <- init_var
  result$rho[,1] <- init_rho

  for(i in 2:m_iter) {
    #set priors
    temp_mu <- result$mu[,,i-1]
    temp_sigma <- result$sigma[,,,i-1]
    temp_rho <- result$rho[,i-1]
    if(theta_out) {
      temp_theta <- result$theta[,,i-1]
    } else {
      if(i == 2){
        temp_theta <- init_theta
      } else{
        temp_theta <- new_theta
      }
    }
    if(i == 2) {
      temp_lambda <- lambda
    } else{
      temp_lambda <- new_lambda
    }

    #Updating new parameters
    new_mu <- sampling_mu(z, temp_theta, temp_mu, temp_sigma)
    new_sigma <- sampling_sigma(z, new_mu, Omega, temp_theta, v_k)
    new_theta <- sampling_theta(x,new_mu,new_sigma,V,z)
    new_rho <- sampling_rho(z, k)
    sampled_z <- sampling_z(new_theta, new_rho, new_mu, new_sigma, z, temp_lambda)
    z <- sampled_z$z
    new_lambda <- sampled_z$lambda


    result$mu[,,i] <- new_mu
    result$sigma[,,,i] <- new_sigma
    result$rho[,i] <- new_rho
    if(theta_out) {
      result$theta[,,i] <- new_theta
    }
    if(i > burn_in) {
      result$lambda <- result$lambda + new_lambda #lambda is added after burn-in period
    }
    if(z_out){
      result$z[,i] <- z
    }
  }
  result$lambda <- result$lambda/iters #average of lambda produced after burn-in
  result$burn <- burn_in
  result$m_iter <- m_iter

  if(!z_out){
    result$z <- sapply(1:n,FUN=function(i){
      which.max(result$lambda[i,]) #Maximum a Posteriori estimate based on lambda estimates
    })
  }
  return(result)
}
