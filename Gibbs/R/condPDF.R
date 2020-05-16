# Contents for 'ConPDF.R'

#' Conditional log PDF of theta
#'
#' @param y matrix `n x p` of observation
#' @param theta matrix `n x p` of observation mean
#' @param z vector length `n` of group membership
#' @param mu matrix `k x p` of group mean
#' @param Sig matrix `p x p x k` of within-group variance
#' @param V matrix `p x p x n` of observation variance
#' @param log logical value whether to output log density or not. it defaults to `FALSE`
#' @param sum logical value whether to output summation of densities of observations. It defaults to `FALSE`
#' @return vector length `n` of (log) density of given `theta` values when `sum=FALSE`. It produces integer that sums up all (log) densities of `theta`
#' @export
dtheta <- function(y, theta, z, mu, Sig, V, log=FALSE, sum=FALSE) {
  n <- nrow(y)
  p <- ncol(y)
  k <- nrow(mu)

  mean_mat <- matrix(0,n,p)
  var_mat <- array(0,dim=c(p,p,n))
  for(i in 1:n) {
    z_i <- z[i]
    temp <- V[,,i] + Sig[,,z_i]
    inv_mat <- solveV(temp)
    G_mat <- crossprod(Sig[,,z_i],inv_mat)
    err <- y[i,] - mu[z_i,]
    mean_mat[i,] <- G_mat%*%err + mu[z_i,]
    var_mat[,,i] <- G_mat%*%V[,,i]
  }
  if(sum){
    return(sum(dmNorm(theta,mean_mat,var_mat,log=log)))
  } else{
  return(dmNorm(theta,mean_mat,var_mat,log=log))
  }
}

#' Conditional log PDF of mu
#'
#' @param theta matrix `n x p` of observation mean
#' @param z vector length `n` of group membership
#' @param mu matrix `k x p` of group mean
#' @param Sig matrix `p x p x k` of within-group variance
#' @param log logical value whether to output log density or not. it defaults to `FALSE`
#' @param sum logical value whether to output summation of densities of observations. It defaults to `FALSE`
#' @return vector length `n` of (log) density of given `mu` values when `sum=FALSE`. It produces integer that sums up all (log) densities of `mu`
#' @export
dmu <- function(theta, z, mu, Sig, log=FALSE, sum=FALSE) {
  n <- nrow(theta)
  p <- ncol(theta)
  k <- dim(Sig)[3]

  mean_mat <- matrix(0,k,p)
  var_mat <- array(0,dim=c(p,p,k))
  for(i in 1:k) {
    ind <- which(z==i)
    mean_mat[i,] <- colSums(theta[ind,])/length(ind)
    var_mat[,,i] <- Sig[,,i]/length(ind)
  }
  if(sum) {
    return(sum(dmNorm(mu,mean_mat,var_mat,log=log)))
  } else{
  return(dmNorm(mu,mean_mat,var_mat,log=log))
  }
}


#' Conditional log PDF of Sigma
#'
#' @param theta matrix `n x p` of observation mean
#' @param z vector length `n` of group membership
#' @param mu matrix `k x p` of group mean
#' @param Sigma matrix `p x p x k` of within-group variance
#' @param Omega matrix `p x p x k` of between row precision
#' @param v_k vector length `k` of prior contribution to each cluster
#' @param log logical value whether to output log density or not. it defaults to `FALSE`
#' @param sum logical value whether to output summation of densities of observations. It defaults to `FALSE`
#' @return vector length `n` of (log) density of given `Sigma` values when `sum=FALSE`. It produces integer that sums up all (log) densities of `Sigma`
#' @export
dSig <- function(theta, z, mu, Sigma, Omega, v_k, log=FALSE, sum=FALSE) {
  n <- nrow(theta)
  p <- ncol(theta)
  k <- nrow(mu)

  psi <- array(0, dim=c(p,p,k))
  df <- rep(0,k)
  for(i in 1:k) {
    ind <- which(z==i)
    temp <- theta[ind,]-mu[i,]
    psi[,,i] <- Omega[,,i]
    for(j in 1:dim(temp)[1]){
      if(length(temp)==p){
        psi[,,i] <- psi[,,i] + crossprod(t(temp))
      } else{
        psi[,,i] <- psi[,,i] + crossprod(t(temp[j,]))
      }
    }
    df[i] <- length(ind) + v_k[i]
  }
  if(sum) {
    return(sum(diwish(Sigma,psi,df,log=log)))
  } else{
  return(diwish(Sigma,psi,df,log=log))
  }
}

#' Conditional log PDF of rho
#'
#' @param z vector length `n` of group membership
#' @param rho vector length `k` of prior group membership probability
#' @param log logical value whether to output log density or not. it defaults to `FALSE`
#' @param sum logical value whether to output summation of densities of observations. It defaults to `FALSE`
#' @return vector length `n` of (log) density of given `rho` values when `sum=FALSE`. It produces integer that sums up all (log) densities of `rho`
#' @export
drho <- function(z, rho, log=FALSE, sum=FALSE) {
  k <- length(rho)

  alpha <- rep(0,k)
  for(i in 1:k) {
    ind <- which(z==i)
    alpha[i] <- length(ind)+1
  }
  if(sum){
    return(sum(ddirichlet(rho,alpha,log=log)))
  } else{
    return(ddirichlet(rho,alpha,log=log))
  }
}


#' Conditional log PDF of z
#'
#' @param z vector length `n` of group membership
#' @param rho vector length `n` of prior group membership probability
#' @param theta matrix `n x p` of observation mean
#' @param mu matrix `k x p` of group mean
#' @param Sigma matrix `p x p x k` of within-group variance
#' @param log logical value whether to output log density or not. it defaults to `FALSE`
#' @param sum logical value whether to output summation of densities of observations. It defaults to `FALSE`
#' @return vector length `n` of (log) density of given `z` values when `sum=FALSE`. It produces integer that sums up all (log) densities of `z`
#' @export
dz <- function(z, rho, theta, mu, Sigma, log=FALSE, sum=FALSE) {
  n <- nrow(theta)
  p <- ncol(theta)
  k <- length(rho)
  res <- rep(0,n)
  log_rho <- log(rho)
  inv_sig <- array(0,dim=c(p,p,k))
  log_det <- rep(0,k)
  lambda_mat <- matrix(0,n,k)
  memb <- matrix(0,n,k)
  for(i in 1:k){
    inv_sig[,,i] <- solveV(Sigma[,,i])
    log_det[i] <- log(det(Sigma[,,i]))
  }
  for(i in 1:n){
    z_i <- z[i]
    memb[i,z_i] <- 1
    for(j in 1:k){
      Siga <- crossprod(inv_sig[,,j],(theta[i,]-mu[j,]))
      aSiga <- crossprod((theta[i,]-mu[j,]),Siga)
      lambda_mat[i,j] <- log_rho[j] - 1/2*(aSiga + log_det[j])
    }
    lambda_mat[i,] <- exp(lambda_mat[i,])
    lambda_mat[i,] <- lambda_mat[i,]/sum(lambda_mat[i,])
    res[i] <- dmultinom(memb[i,],size=1,prob=lambda_mat[i,],log=log)
  }
  if(sum) {
    return(sum(res))
  } else{
    return(res)
  }
}

#' Unnormalized complete data log-posterior
#'
#' @param mu matrix `k x p` of group mean
#' @param Sigma matrix `p x p x k` of within-group variance
#' @param rho vector length `k` of prior group membership probability
#' @param y matrix `n x p` of observation
#' @param V matrix `p x p x n` of observation variance
#' @param theta matrix `n x p` of observation mean
#' @param z vector length `n` of group membership
#' @param Omega matrix `p x p x k` of between row precision
#' @param v_k vector length `k` of prior contribution to each cluster
#' @return integer value of log posterior for given parameters
#' @export
log_post <- function(mu, Sigma, rho, y, V, theta, z, Omega, v_k) {
  p <- ncol(mu)
  k <- nrow(mu)

  ll <- nnm_loglik(mu,Sigma,rho,y,V,theta,z)
  for(i in 1:k) {
    log_det <- log(det(Sigma[,,i]))
    inv_sig <- solveV(Sigma[,,i])
    trc <- sum(diag(crossprod(inv_sig,Omega[,,i])))
    ll <- ll - 1/2 * ((v_k[i]+p+1)*log_det + trc)
  }
  return(ll)
}
