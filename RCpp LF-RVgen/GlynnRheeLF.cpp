// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// Drift of X
arma::vec  mu_Xcpp(const arma::vec& theta, const arma::vec& X) {
    
    double alpha(theta[0]), beta(theta[1]), kappa(theta[2]), vbar(theta[3]);
    arma::vec mu(2);
    
    mu[0] = alpha + beta*X[1];
    mu[1] = kappa*(vbar - X[1]);
    
    return mu;
  }
  
// Volatility of X 
arma::mat  sigma_Xcpp(const arma::vec& theta, const arma::vec& X) {
    
    double rho = theta[4];
    double sigmav = theta[5];
    
    arma::mat mat(2,2);
    mat(0,0)=1;
    mat(1,0) = rho*sigmav;
    mat(0,1) = 0;
    mat(1,1) = sigmav*pow((1-pow(rho,2)),0.5);
    
    return mat*pow(X[1],0.5);
    
  }
  
 
// Variance of X
arma::mat var_X(const arma::vec& theta, const arma::vec& X){
    
    double rho = theta[4];
    double sigmav = theta[5];
    
    arma::mat mat(2,2);
    mat(0,0)=1;
    mat(1,0) = rho*sigmav;
    mat(0,1) = rho*sigmav;
    mat(1,1) = pow(sigmav,2);
    
    return X[1]*mat;
    
  }

// psi function as defined in the paper 
arma::colvec  psi_Xcpp(const arma::vec& theta, const arma::vec& x, const arma::vec& rho) {
    
    arma::colvec omega(2);
    omega = mu_Xcpp(theta,x)-rho;
    
    return inv(sigma_Xcpp(theta,x))*omega;
    
  }
  
// Jump intensity
double lambda_Xcpp(const arma::vec& theta, const arma::vec& x) {
    
    //double lambda0(theta[6]), lambda1(theta[7]), lambda2(theta[8]);
    
    //return lambda0 + lambda1*x[0] + lambda2*x[1];
    return theta[6];
  }
  
// "Gamma" in the paper - Jump size in returns
arma::vec gamma_X(const arma::vec& theta, const arma::vec& z) {
    
    double gamma0x1(theta[7]), gamma1x1(theta[8]);
    //double gamma0x2(theta[11]), gamma1x2(theta[12]);
    arma::vec gammaX(2);
    
    gammaX[0] = gamma0x1 + gamma1x1*z[0];
    //gammaX[1] = gamma0x2 + gamma1x2*z[1];
    gammaX[1] = z[1];
    
    return gammaX;
    
  }
  
// It's the part between p_hat and first exp
double  Part2cpp(const arma::vec& theta, const arma::vec& x, double ell){
    
    return lambda_Xcpp(theta,x)/ell;
    
  }
  

// It's the second exp part
double  Part3cpp(double ell, const arma::vec& theta, const arma::vec& rho, const arma::vec& x, const arma::vec& Wt, int n, double time){
  
    double prod = as_scalar(trans(psi_Xcpp(theta,x,rho))*psi_Xcpp(theta,x,rho));
    return (ell-lambda_Xcpp(theta,x)-0.5*prod)*time/(pow(2,n)) + as_scalar(trans(psi_Xcpp(theta,x,rho))*Wt);

  }
  
// Compute the discretization between 2 jumps 
arma::mat JumpDcpp(const arma::vec& theta, double Delta, double T1, double T2, arma::mat Xstar, int n, const arma::mat& Wt, double ell, const arma::vec& rhod){
    
    //NumericVector wt(2);
    arma::colvec wt(2);
    double h(0);

    //NumericVector X(2);
    arma::vec X(2);

    if (n!=0){
      
      arma::mat Wt2(2,pow(2,n-1));
      for (int i = 0; i< Wt2.n_cols; i++) {
        Wt2(0,i) = Wt(0,2*i) + Wt(0,(2*i+1));
        Wt2(1,i) = Wt(1,2*i) + Wt(1,(2*i+1));
      }
      for (int i=n-1; i<(n+1); i++){
        h=(T2-T1)/pow(2,i);
        X[0] = Xstar(i-n+1,0);
        X[1] = Xstar(i-n+1,1);
        //std::cout<<X[0]<<" : "<<X[1]<<std::endl;
        for (int j=0; j<pow(2,i);j++){
          if(i==n-1){
            wt[0] = Wt2(0,j);
            wt[1] = Wt2(1,j);
            Xstar(0,3)= Xstar(0,3) + Part3cpp(ell,theta,rhod,X,wt,i,T2-T1);
            X = X + rhod*h + sigma_Xcpp(theta,X)*wt;
            if(X[1]<0){ X[1]=0.001;}
            Xstar(0,0) = X[0];
            Xstar(0,1) = X[1];
          }else{
            wt[0] = Wt(0,j);
            wt[1] = Wt(1,j);
            Xstar(1,3)= Xstar(1,3) + Part3cpp(ell,theta,rhod,X,wt,i,T2-T1);
            X = X + rhod*h + sigma_Xcpp(theta,X)*wt;
            if(X[1]<0){ X[1]=0.001;}
            Xstar(1,0) = X[0];
            Xstar(1,1) = X[1];
          }
        }
        if(std::abs(Delta-T2)>0.00001){
          Xstar(i-n+1,2) = Xstar(i-n+1,2)*Part2cpp(theta,X,ell);
        }
      }
    }else{
      h=T2-T1;
      X[0] = Xstar(1,0);
      X[1] = Xstar(1,1);
      wt[0] = Wt(0,0);
      wt[1] = Wt(1,0);
      Xstar(1,3) = Xstar(1,3) + Part3cpp(ell,theta,rhod,X,wt,0,T2-T1);
      X = X + rhod*h + sigma_Xcpp(theta,X)*wt;
      if(X[1]<0){ X[1]=0.001;}
      Xstar(1,0) = X[0];
      Xstar(1,1) = X[1];
      if(std::abs(Delta-T2)>0.00001){
          Xstar(1,2) = Xstar(1,2)*Part2cpp(theta,X,ell);
      }
    }
    
    return Xstar;
    
  }
  
double  Part1cpp(const arma::vec& theta, const arma::vec& XTN, const arma::vec& rhod, double Gap, const arma::mat& Wt, const arma::vec& w, int N2){

    double pi = 3.141592653589793238462643383280;
    arma::vec X_tilda(XTN);
    double h(Gap/N2);
    int dim(w.size());
    double pdf(0);
    arma::colvec wt(dim);

    // Simulate the path of X
    for (int i=0; i<(N2-1); i++) {
      wt[0] = Wt(0,i);
      wt[1] = Wt(1,i);
      X_tilda = X_tilda + rhod*h + sigma_Xcpp(theta,X_tilda)*wt;
    }

    // Compute the pdf with the multivariate normal distribution
    arma::mat V = var_X(theta,X_tilda)*h;
    double det = arma::det(V);

    arma::vec omega(dim);
    omega = w - X_tilda - rhod*h;
    pdf = 1/pow(pow(2*pi,dim)*det,0.5)*exp(-0.5*as_scalar(trans(omega)*(arma::inv(V))*omega));
  
    
    return pdf;
    
  }

//// [[Rcpp::export]]
double DELTAcpp(double Delta, int N, const arma::vec& theta, const arma::vec& rho, int ell, const arma::vec& v, const arma::vec& w, int N2, const arma::vec& Npois, const arma::mat& Z, const arma::vec& W1, const arma::vec& W3){
  
    // Nd discretization steps
  int M(pow(2,N));
  // Dimension of the process
  int d(v.n_elem);
  
  // Random jump related part
  int Nj(Npois[0]);
  NumericVector tau(Nj+2);
  for(int i=0;i<(Nj+2);i++){
   tau[i]=Npois[i+1];
   //std::cout<<tau[i]<<std::endl;
  }
  
  // X_star is used to store results of Y, V, Part2, 3, 4
  arma::mat X_star(2,d+2);
  arma::mat XTN(2,d);
  for(int i=0;i<d;i++){
    X_star(0,i) = v[i];
    XTN(0,i) = v[i];
    X_star(1,i) = v[i];
    XTN(1,i) = v[i];
  }
  X_star(0,d) = 1;
  X_star(1,d) = 1;
  X_star(0,d+1) = 0;
  X_star(1,d+1) = 0;
 
  arma::mat W(d,(Nj+1)*M);
  //arma::vec W1(arma::randn(d*(Nj+1)*M));
  for(int i=0;  i<d; i++){
    for(int j=0; j<(Nj+1)*M; j++){
      W(i,j) = W1[j+i*(Nj+1)*M];
      //std::cout<<W(i,j)<<std::endl;
    }
  }
  arma::mat Wt(d,M);
  
  //Here are the Euler discretizations between jumps
  if (Nj!=0){
    for (int i=0; i < (Nj+1); i++){
      double sd = pow((tau[i+1]-tau[i]),0.5)/M;
      //General Brownian Motions with sd=1
      for(int j  = 0 ; j < M ; j++){
        Wt(0,j) = sd*W(0,j+i*M);
        Wt(1,j) = sd*W(1,j+i*M);
      }
      
      // discretize between jump times
      X_star = JumpDcpp(theta,Delta,tau[i],tau[i+1],X_star,N,Wt,ell,rho);
      if (i<Nj){
        arma::vec Zj(d);
        Zj[0] = Z(0,i);
        Zj[1] = Z(1,i);
        arma::vec jump(gamma_X(theta,Zj));
        for(int k=0; k<d; k++){
          X_star(0,k) = X_star(0,k) + jump[k];
          XTN(0,k) = X_star(0,k);
          X_star(1,k) =  X_star(1,k) + jump[k];
          XTN(1,k) = X_star(1,k);
        }
      }
    }
  }else{
    double sd = pow(Delta,0.5)/M;
    //General Brownian Motions with sd=1
    for(int j  = 0 ; j < M ; j++){
        Wt(0,j) = sd*W(0,j);
        Wt(1,j) = sd*W(1,j);
      }
    // Discretize between 0 and T
    X_star = JumpDcpp(theta,Delta,0,Delta,X_star,N,Wt,ell,rho);
  }

  // Part2 values
  NumericVector Part2_V(2); 
  Part2_V[0] = X_star(0,d);
  Part2_V[1] = X_star(1,d);
  // Part3 values
  NumericVector Part3_V(2);
  Part3_V[0] = exp(X_star(0,d+1));
  Part3_V[1] = exp(X_star(1,d+1));

  // Time Gap between last jump and T
  double Gap;
  if (Nj!=0){
     Gap = Delta-tau[Nj];
  }else{ Gap=Delta;}
  
  arma::mat Wt3(d,N2);
  //arma::vec W3(arma::randn(d*N2));
  double std = pow(Gap/N2,0.5);
  for(int k=0;k<d;k++){
    for(int i=0; i<N2; i++){
      Wt3(k,i) = std*W3[i+k*N2];
    }
  }
  
  if (N!=0){
    arma::vec XNm1(d);
    arma::vec XN(d);
    for(int k=0;k<d;k++){
     XNm1[k] = XTN(0,k);
     XN[k] = XTN(1,k);
    }
    double delta1(Part1cpp(theta, XNm1, rho, Gap, Wt3, w,N2)*Part2_V[0]*Part3_V[0]);
    double delta2(Part1cpp(theta, XN, rho, Gap, Wt3, w,N2)*Part2_V[1]*Part3_V[1]);
    return (delta2-delta1);
  }else{
    arma::vec XN(d);
    for(int k=0;k<d;k++){
     XN[k] = XTN(1,k);
    }
    double delta2(Part1cpp(theta, XN, rho, Gap, Wt3, w,N2)*Part2_V[1]*Part3_V[1]);
    return delta2;
  }
  
}
   
arma::vec LFVecCpp(double Delta, int N, const arma::vec& theta, const arma::vec& rho, int ell, const arma::mat& X, int N2, const arma::vec& Npois, const arma::mat& Z, const arma::vec& W1, const arma::vec& W3){

  int Tf = X.n_cols;
  int d = X.n_rows;
  arma::vec LF(Tf-1); 
  for(int i=0; i< (Tf-1); i++){
    arma::vec Xt(d);
    arma::vec Xtp1(d);
    for(int j=0; j<d;j++){
      Xt[j]=X(j,i);
      Xtp1[j]=X(j,i+1);
    }

    LF[i] = DELTAcpp(Delta,N,theta,rho,ell,Xt,Xtp1,N2, Npois, Z,W1,W3);
  }
 
  return LF;

}

// [[Rcpp::export]] 
arma::vec LFVec(double Delta, int N, NumericVector theta, arma::vec rho, int ell, arma::mat X, int N2, arma::vec Npois, arma::mat Z, arma::vec W1, arma::vec W3){

  int Tf = X.n_cols;
  int d = X.n_rows;
  arma::vec LF(Tf-1); 
  for(int i=0; i< (Tf-1); i++){
    arma::vec Xt(d);
    arma::vec Xtp1(d);
    for(int j=0; j<d;j++){
      Xt[j]=X(j,i);
      Xtp1[j]=X(j,i+1);
    }

    LF[i] = DELTAcpp(Delta,N,theta,rho,ell,Xt,Xtp1,N2, Npois, Z,W1,W3);
  }
 
  return LF;

}

// [[Rcpp::export]]  
double LFcpp(arma::vec ubounds, arma::vec lbounds, int M, double Delta, int N, arma::vec theta, arma::vec rho, int ell, arma::mat X, int N2, arma::mat Nrv, arma::mat Npois, arma::mat Marks, arma::mat W1, arma::mat W3){

arma::vec hit_u(ubounds-theta);
arma::vec hit_l(theta-lbounds);
double hit_ub = hit_u.min();
double hit_lb = hit_l.min();
std::cout<<hit_ub<<" : "<<hit_lb<<std::endl;
double llf(0);

std::cout<<"Theta"<<std::endl;
std::cout<<theta[0]<<" : "<<theta[1]<<" : "<<theta[2]<<" : "<<theta[3]<<" : "<<theta[4]<<" : "<<theta[5]<<" : "<<theta[6]<<" : "<<theta[7]<<" : "<<theta[8]<<" : "<<theta[9]<<std::endl;



if(hit_ub>0 && hit_lb>0){

int Tf(X.n_cols);
int d(X.n_rows);
arma::vec Lf(arma::zeros(Tf-1));
arma::vec Np(Npois.n_cols);
arma::mat Z(d,Npois.n_cols-1);
double muv(theta[9]);

double rho_c = theta[4];
double sigmav = theta[5];
arma::mat mat(2,2);
mat(0,0)=1;
mat(1,0) = rho_c*sigmav;
mat(0,1) = 0;
mat(1,1) = sigmav*pow((1-pow(rho_c,2)),0.5);

double deter(arma::det(mat));
std::cout<<"determinant: "<<deter<<std::endl;



for (int m=0; m<M; m++){

  int n(Nrv(m,0));
  int Nj(Npois(m,0));

  for(int i=0;i<Npois.n_cols;i++){
    Np[i] = Npois(m,i);
    if(i<Nj){
      Z(0,i) = Marks(2*m,i);
      Z(1,i) = -log(Marks(2*m+1,i))*muv;
    }
  }

  int t1(pow(2,n)*d*(Nj+1));
  arma::vec BM1(t1);
  for(int i=0;i<t1;i++){
    BM1[i] = W1(m,i);
  }

  arma::vec BM3(d*N2);
  for(int i=0;i<(d*N2);i++){
    BM3[i] = W3(m,i);
  }
  
  Lf = Lf + LFVecCpp(Delta,n,theta,rho,ell,X,N2,Np,Z,BM1,BM3)/(M*Nrv(m,1));
}

//double llf(0);
for(int i=0;i<(Tf-1);i++){
  llf=llf+log(Lf[i]);
}

} else{
  llf = -10000000;
}
std::cout<<"llf: "<<llf<<std::endl;

return llf;

}


