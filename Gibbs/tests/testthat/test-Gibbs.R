# Contents of file `test-Gibbs.R`
context("GibbsSampler")

require(Gibbs)
test_that("conditional log PDF and complete data log-posterior differs by constant", {
  normalize <- function(x) {
    return(x/sum(x))
  }
  #generate test case
  n <- sample(50:100,1)
  p <- sample(2:5,1)  # x = (n x p) matrix
  k <- sample(2:5,1)  # number of clusters
  rho <- normalize(runif(k,1,3)) # prior membership distribution
  z <- sample(1:k, n, replace=T, prob=rho) # group membership
  mu <- matrix(rnorm(k*p, 0,10), k, p) # group mean

  sig <- array(0,dim=c(p,p,k)) # within-group variance
  for(i in 1:k){
    sig[,,i] <- crossprod(rMNorm(1,matrix(0,p,p)))
  }

  V <- array(0,dim=c(p,p,n)) # observation variance
  for(i in 1:n){
    V[,,i] <- crossprod(rMNorm(1,matrix(0,p,p)))/10
  }

  x <- matrix(0,n,p) # observation
  for(i in 1:n){
    cluster <- z[i]
    x[i,] <- rmNorm(1,mu[cluster,],sig[,,cluster])
  }
  m_iter <- sample(300:500,1) #iterations
  Omega <- array(var(x),dim=c(p,p,k))
  v_k <- rep(p+2,k)

  gbs <- gibbs(x,k,m_iter,V=V,theta_out=T,z_out=T)

  ind <- 1
  ii <- sample(1:gbs$m_iter,1) # which iteration to fix conditioned parameters
  jj <- sample(3:6,1) # which iteration to change parameter values
  results <- list('theta'=rep(NA,jj),'mu'=rep(NA,jj),
                  'sig'=rep(NA,jj),'rho'=rep(NA,jj),
                  'z'=rep(NA,jj)) # list of parameters to compare with complete log likelihood
  for(j in sample(1:gbs$m_iter,jj)) {
    X <- as.numeric(j)
    # compute complete log posterior - conditional theta
    lp_theta <- log_post(gbs$mu[,,ii],gbs$sigma[,,,ii],gbs$rho[,ii],
                         x,V,gbs$theta[,,X],gbs$z[,ii],Omega,v_k)
    theta_lpdf <- dtheta(x,gbs$theta[,,X],gbs$z[,ii],gbs$mu[,,ii],
                         gbs$sigma[,,,ii],V,log=TRUE,sum=TRUE)
    results$theta[ind] <- lp_theta - theta_lpdf

    # compute complete log posterior - conditional mu
    lp_mu <- log_post(gbs$mu[,,X],gbs$sigma[,,,ii],gbs$rho[,ii],
                      x,V,gbs$theta[,,ii],gbs$z[,ii],Omega,v_k)
    mu_lpdf <- dmu(gbs$theta[,,ii],gbs$z[,ii],gbs$mu[,,X],
                   gbs$sigma[,,,ii],log=TRUE,sum=TRUE)
    results$mu[ind] <- lp_mu - mu_lpdf

    # compute complete log posterior - conditional Sigma
    lp_sig <- log_post(gbs$mu[,,ii],gbs$sigma[,,,X],gbs$rho[,ii],
                      x,V,gbs$theta[,,ii],gbs$z[,ii],Omega,v_k)
    sig_lpdf <- dSig(gbs$theta[,,ii],gbs$z[,ii],gbs$mu[,,ii],
                     gbs$sigma[,,,X],Omega,v_k,log=TRUE,sum=TRUE)
    results$sig[ind] <- lp_sig - sig_lpdf

    # compute complete log posterior - conditional rho
    lp_rho <- log_post(gbs$mu[,,ii],gbs$sigma[,,,ii],gbs$rho[,X],
                      x,V,gbs$theta[,,ii],gbs$z[,ii],Omega,v_k)
    rho_lpdf <- drho(gbs$z[,ii],gbs$rho[,X],log=TRUE,sum=TRUE)
    results$rho[ind] <- lp_rho - sig_lpdf

    # compute complete log posterior - conditional z
    lp_z <- log_post(gbs$mu[,,ii],gbs$sigma[,,,ii],gbs$rho[,ii],
                      x,V,gbs$theta[,,ii],gbs$z[,j],Omega,v_k)
    z_lpdf <- dz(gbs$z[,X],gbs$rho[,ii],gbs$theta[,,ii],gbs$mu[,,ii],
                 gbs$sigma[,,,ii],log=TRUE,sum=TRUE)
    results$z[ind] <- lp_z - z_lpdf

    ind <- ind + 1
  }

  tol <- 1e-6
  expect_equal(length(diff(results$theta)<tol),jj-1)
  expect_equal(length(diff(results$mu)<tol),jj-1)
  expect_equal(length(diff(results$sig)<tol),jj-1)
  expect_equal(length(diff(results$rho)<tol),jj-1)
  expect_equal(length(diff(results$z)<tol),jj-1)

})
