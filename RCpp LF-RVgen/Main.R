# This code compute the likelihood function 
setwd("~/Desktop/RCpp LF-RVgen")
source("DeltaRcpp.R")
source("X_generation.R")
rho=c(0,0)
ell=40
Delta=1/12
T=4
theta=c(0.05,-0.5,0.1,0.2,-0.5,0.2,25,-0.01,0.02,0.02)
#theta=c(0.05,-0.5,0.1,0.2,-0.5,0.2,25,0,0,-0.01,0.02,0,1,0.02)

# First example with 3 points
# X = matrix(0,2,3)
# X[,1] = c(0,0.1)
# X[,2] = c(0,0.11)
# X[,3] = c(-0.3,0.17)
# Second example with more points
#X = matrix(0,2,61)
#X[2,] = seq(from=0.1, to=0.4, by = 0.005)
X=gen_X(theta,T,1/12,c(4.5,0.1))



N=5 #This is the max N
M = 20000 #MC simulations
N2 = 5 # discretization steps after last jump time
d = nrow(X)

Nrv = gen_N(N,M)
Npois = gen_P(Delta,M,ell)
Marks = gen_Marks(theta,M,Npois)
W1 = gen_Norm1(N,M,Npois,Nrv,d)
W3 = gen_Norm3(N2,M,d)

#Monte-Carlo part
#start=Sys.time()
# Compute the log-likelihood function using the R-C++ function
#LF = LF(M,Delta,N,theta,rho,ell,X,N2,Nrv,Npois,Marks,W1,W3)
# Compute the log-likelihood function using the C++ function
LF = LFcpp(ubounds,lbounds,M,Delta,N,theta,rho,ell,X,N2,Nrv,Npois,Marks,W1,W3)

end=Sys.time()
time1=end-start

lbounds = c(-0.2,-1,0  ,0  ,-1,  0,0,-0.05,0   ,0)
ubounds = c(0.2,  1,0.5,0.5,1 ,0.5,50,0.05,0.1,0.1 )



fr<-function(theta_hat){

    logllf = -1*LFcpp(ubounds,lbounds,M,Delta,N,theta_hat,rho,ell,X,N2,Nrv,Npois,Marks,W1,W3)

    if (is.na(logllf)){
      logllf=10000000
    }
    
    logllf
  
}

theta_false=theta*1.2
# mll=optim(par = theta_false, fn = fr, method="Nelder-Mead", control = list(fnscale=-1,maxit=1000,REPORT=1))
# 
# print(c(mll$value, mll$counts, mll$par))

mll=fminbnd(fun = fr,x0 = theta_false,xmin = lbounds,xmax = ubounds,options = list(MaxIter=100,MaxFunEvals=100000000))

