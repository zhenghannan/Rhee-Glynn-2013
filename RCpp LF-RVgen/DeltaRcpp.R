library(Rcpp)
library(RcppArmadillo)
sourceCpp('GlynnRheeLF.cpp')

gen_N = function(N,M){
  
  quantile = seq(from = 0, to = N, by = 1) 
  Nprob<-dpois(quantile,2)  
  Mat = matrix(0, nrow=M, ncol=2)
  
  for(m in 1:M){
    n <- sample(quantile,1,prob=Nprob,replace=TRUE)
    Mat[m,1] = n
    Mat[m,2] = Nprob[n+1]
  }
  
  Mat
  
}

gen_P = function(Delta, M, ell) {
  # tau is the length of the intervals, M is the # of MC simulation
  
  #P of the algorithm
  P_here = rpois(M, ell*Delta) # generate a poisson process of M length, l*Delta intensity
  
  res = matrix(Delta, nrow = M, ncol = max(P_here)+3) # construct a matrix with tau in it 
  res[,1] = P_here # replace the first column by the poisson realization
  res[,2] = rep(0,time=M)
  for (m in 1:M) {
    p = res[m,1] # number of poisson jumps for the montecarlo j
    if (p > 0) {
      # if there is at least a jump, fill in the row with the jump times between 0 and tau
      # these jump times are drawn through a uniform
      res[m,3:(p+2)] = sort(runif(p,0,Delta)) 
    }
  }
  
  res
}

gen_Marks = function(theta,M, Npois){
  
  Nmax = max(Npois[,1])
  muv = theta[14]
  res = matrix(0, nrow = M*2, ncol = Nmax )
  for(m in 1:M){
    res[2*m-1,] = rnorm(Nmax)
    res[2*m,] = runif(Nmax)
  }
  
  res
}

gen_Norm1 = function(N,M, Npois,Nrv, d){
  
  Nmax = max(Npois[,1])
  res = matrix(NA, nrow=M, ncol = (Nmax+1)*d*(2^N))
  
  for(m in 1:M){
    Nj = Npois[m,1]
    n = Nrv[m,1]
    res[m,1:(d*(Nj+1)*(2^n))] = rnorm(d*(Nj+1)*(2^n))
  }

  res
  
}

gen_Norm3 = function(N2,M,d){
  
  res = matrix(NA, nrow=M, ncol = N2*d)
  
  for(m in 1:M){
    res[m,] = rnorm(d*N2)
  }
  
  res
  
}

# Compute the likelihood function using the R loop
LF = function(M,Delta,N,theta,rho,ell,X,N2,Nrv,Npois,Marks,W1,W3){
  
  Lf = matrix(0, nrow = M, ncol = (ncol(X)-1))  
  
  for (m in 1:M){
    Nj = Npois[m,1]
    n=Nrv[m,1]
    Lf[m,] = LFVec(Delta,n,theta,rho,ell,X,N2,Npois[m,],Marks[(2*m-1):(2*m),],W1[m,1:(2^n*d*(Nj+1))],W3[m,])/Nrv[m,2]
  }

  lfvec = log(colSums(Lf)/M)
  
  sum(lfvec)
  
}
