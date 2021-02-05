gen_X=function(theta,T,Delta,X0){
  
  # Drift of X
  mu_X = function(theta, X) {
    mu=c(theta[1]+theta[2]*X[2],theta[3]*(theta[4]-X[2]))
  }
  
  # Volatility of X
  sigma_X = function(theta,X) {
    
    rho = theta[5]
    sigmav = theta[6]
    
    mat = matrix(0,nrow=2,ncol=2)
    mat[1,] = c(1, 0)
    mat[2,] = c(rho*sigmav, sqrt(1-rho^2)*sigmav)
    
    sqrt(X[2])*mat
    
  }
  
  # Jump size of the process
  gamma_X = function(theta, x, z) {
    
    gamma0x1 = theta[8]
    gamma1x1 = theta[9]
    
    gamma=c(0,0)
    gamma[1] = gamma0x1 + gamma1x1*z[1]
    gamma[2] = z[2]
    
    gamma
    
  }

  
  lambda = theta[7]
  M = rpois(1, lambda*T)
  U = sort(runif(M,0,T))
  D = rnorm(M, 0, 1)
  muv = theta[10]
  
  Euler_times=seq(from=0, to=T, by=Delta)
  p = findInterval(U, Euler_times)
  
  Z=matrix(0,nrow=2,ncol=length(D))
  Z[1,]=D
  Z[2,]=rexp(M,1/muv)
  
  Norm = matrix(0,nrow=(T/Delta+M),ncol = 2)
  Norm[,1] = rnorm((T/Delta+M))
  Norm[,2] = rnorm((T/Delta+M))
  
  X=matrix(0,2,T/Delta+1)
  X[,1] = X0
  
  k=1
  n=1
  time<-c()
  for (i in 1:(T/Delta)){
    count=sum(p==i)
    if(count!=0){
      t=0
      Xj = X[,i]
      for (j in 1:count){
        t=U[k]-t-(i-1)/12
        time=c(time,t)
        Xj = Xj + mu_X(theta,Xj)*t + sqrt(t)*sigma_X(theta,Xj)%*%Norm[n,] + gamma_X(theta,Xj,Z[,k])
        k=k+1
        n=n+1
        Xj[2]=max(Xj[2],0)
      }
      t=i/12-U[k-1]
      time=c(time,t)
      X[,i+1] = Xj + mu_X(theta,Xj)*t + sqrt(t)*sigma_X(theta,Xj)%*%Norm[n,]
      X[2,i+1]=max(X[2,i+1],0)
       n=n+1
    }else{
      t=1/12
      time=c(time,t)
      X[,i+1] = X[,i] + mu_X(theta,X[,i])*t + sqrt(t)*sigma_X(theta,X[,i])%*%Norm[n,]
      X[2,i+1]=max(X[2,i+1],0)
      n=n+1
    }
    
  }
  X
}

# X0 = c(4.5,0.1)
# T=4
# Delta=1/12
# Euler_times=seq(from=0, to=T, by=Delta)
# theta=c(0.05,0,0.1,0.2,-0.5,0.2,3,0,0,-0.01,0.02,0,1,0.02)
# sample=gen_X(theta,T,Delta,X0)
# 
# plot(Euler_times,exp(sample[1,]),type='l')
