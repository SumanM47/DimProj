#' @name Slice_to_3D
#' @title Draw posterior samples from the three dimensional process using two dimensional slices
#'
#' @description
#' Generates posterior samples for parameters of the three dimensional process using two dimensional slices
#'
#' @usage Slice_to_3D(data,K,
#' inits=NULL,
#' hyperparams=list("mu0" = rep(0,dim(data$Y)[2]), "rho0" = rep(1,K),
#' "A0" = matrix(0,dim(data$Y)[2],K),
#' "musig2"=0,"sigsig2"=10,
#' "sig2mu"=100,"sig2rho"=100,"sig2A"=100),
#' MCMCparams=list(niters=1.2e4,nthin=10,nburn=2e3))
#'
#' @param data X class object or a list containing the entries array Y (number of slices x number of species x number of points), matrix S (number of locations x 2), vector gridsize (2x1, size of the computation grid) and positive scalar bma (thickness of each slice)
#' @param K positive integer, the number of latent processes to be used. Recommended to be smaller than the number of species
#' @param inits named list containing initial values for the parameters of the model. The components are mu(numeric vector, same length as the number of different species: 3 dimensional constant mean for each species), rho(positive numeric vector, same length as K: 3 dimensional spatial range parameter for each latent process), A(numeric matrix of dimension number of species x K: mapping matrix from the latent processes to the different species), sig2(positive numeric vector, same length as the number of slices: slice specific variability) and tau2(positive numeric vector, same length as the number of difference species: species specific variability). If NULL (default), a quick initial value is selected within the function.
#' @param hyperparams named list containing the hyperparameters for the prior distributions of mu, A, rho, log(sig2). The components are mu0(numeric vector, same length as the number of different species: mean for the normal prior for the mu parameters), rho0(positive numeric vector, same length as K: exponentiated mean for the normal prior for the log(rho) parameters), A0(numeric matrix of dimension number of species x K: mean for the normal prior of the A parameters), musig2(scalar: mean parameter for the normal prior for log sig2), sigsig2(positive scalar: variance parameter for the normal prior for log sig2), sig2mu(positive scalar: variance parameter for the normal prior of the mu parameters), sig2A(positive scalar: variance parameter for the normal prior of the A parameters), sig2rho(positive scalar: variance parameter for the normal prior of the rho parameters).
#' @param MCMCparams named list containing details about the running of MCMC, the components are niters (integer, number of posterior samples to be drawn), nthin(integer, the thinning interval), nburn (integer, number of posterior samples to be discarded as burn-in samples). Total number of iterations = niters*nthin
#'
#' @useDynLib DimProj
#'
#' @import Rcpp
#' @import RcppArmadillo
#' @import stats
#' @import grDevices
#'
#'
#' @return list of posterior samples for the parameters
#' @export

Slice_to_3D <- function(data,K,
                      inits=NULL,
                      hyperparams=list("mu0" = rep(0,dim(data$Y)[2]), "rho0" = rep(1,K), "A0" = matrix(0,dim(data$Y)[2],K),
                                       "musig2"=0,"sigsig2"=10,
                                       "sig2mu"=100,"sig2rho"=100,"sig2A"=100),
                      MCMCparams=list(niters=1.2e4,nthin=10,nburn=2e3)){

  ## Unpacking the data
  Y <- data$Y
  datdim <- dim(Y)
  I <- datdim[1]
  J <- datdim[2]
  gridsize <- data$gridsize
  bma <- data$bma
  S <- as.matrix(data$S)

  ## Initial Values
  A <- is.null(inits)
  Bmu <- is.numeric(inits$mu)
  Bmus <- length(inits$mu) == J
  Brho <- is.numeric(inits$rho) & prod(inits$rho>0)
  Brhos <- length(inits$rho) == K
  BA <- is.matrix(inits$A)
  BAs <- as.logical(prod(dim(inits$A) == c(J,K)))
  Bsig2 <- is.numeric(inits$sig2)
  Bsig2x <- as.logical(prod(inits$sig2 > 0))
  Bsig2s1 <- nrow(inits$sig2) == I
  Bsig2s2 <- ncol(inits$sig2) == J
  B <- as.logical(Bmu*Bmus*Brho*Brhos*BA*BAs*Bsig2*Bsig2x*Bsig2s1*Bsig2s2)
  AB <- as.logical(A+!B)
  if(AB){
    mu_init <- apply(Y,2,"mean",na.rm=T)

    VY <- apply(Y,1:2,"var",na.rm=T)
    sig2_init <- 0.1*VY

    A_init <- matrix(sqrt(0.9*colMeans(VY,na.rm=T)/3),J,K)
    dmax <- sqrt(max(S[,1]-S[1,1])^2 + max(S[,2] - S[1,2])^2)
    rho_init <- rep(0.1*dmax,K)
  } else{
    mu_init <- inits$mu
    rho_init <- inits$rho
    A_init <- inits$A
    sig2_init <- inits$sig2
  }

  ## Need to sourceCpp
  ## Check how for a package

  ll <- Slicesto3D(Y,S,bma,as.integer(gridsize[1]),as.integer(gridsize[2]),
                   mu_init,rho_init,A_init,sig2_init,
                   MCMCparams$niters,MCMCparams$nburn,MCMCparams$nthin,
                   hyperparams$mu0,hyperparams$rho0,hyperparams$A0,hyperparams$musig2,hyperparams$sigsig2,
                   hyperparams$sig2mu,hyperparams$sig2rho,hyperparams$sig2A)

  return(ll)

}
