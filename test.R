library('Rcpp')
library('RcppArmadillo')
setwd("C:\\Users\\Smart\\Desktop\\MoMA-prototype\\Source code")
sourceCpp('sfpca_extended.cpp')

norm_vec <- function(x) sqrt(sum(x^2))

SSD <- function(n){
  a <- 6*diag(n)
  for(i in 1:n){
    for(j in 1:n){
      if(abs(i-j) == 1) a[i,j] = -3;
      if(abs(i-j) == 2) a[i,j] = 1;
    }
  }      
  return(a);
}

uni <- function(n){
  u_1 <- as.vector(rnorm(n)) 
  return(u_1/norm_vec(u_1))
}

#-------------------
# Generate data
#-------------------
n <- 199
p <- 200
ind <- as.vector(seq(p))
u_1 <- uni(n)
u_2 <- uni(n)
u_3 <- uni(n)
eps <- matrix(rnorm(n*p),n,p)
eps <- eps/20

# Sinusoidal
v_1 <- sin((ind+15)*pi/17);v_1[floor(7/20*p):p]=0;
v_1 <- v_1/norm_vec(v_1);
# Gaussian-modulated sinusoidal
v_2 <- as.vector(exp(-(ind-100)^2/650)*sin((ind-100)*2*pi/21)); 
v_2[0:floor(7/20*p)]=0;v_2[floor(130/200*p):p] = 0;
v_2 <- v_2/norm_vec(v_2);
# Sinusoidal
v_3 <- sin((ind-40)*pi/30);v_3[0:floor(130/200*p)]=0;
v_3 <- v_3/norm_vec(v_3);

plot(v_1,type = 'l',ylim=c(-0.3,0.3));
lines(v_2,col='blue');
lines(v_3,col='red')

X <-  n/4*u_1 %*% t(v_1) +eps + n/5*u_2 %*% t(v_2) #+ n/6 * u_3 %*% t(u_3)
# Noise proportion
print(norm(X) /norm(eps))
O_u <- SSD(n)
O_v <- SSD(p)


cnt=1
  for(l in seq(1,10,0.5)){
    res1 <- sfpca("PCA",
                  X=X,
                  Y=matrix(runif(p*n),n,p),
                  Omega_u=O_u,Omega_v=O_v,
                  alpha_u=1,alpha_v=1,
                  lambda_u=l,lambda_v=l,
                  P_u="l1",P_v="l1",
                  non_neg=0,
                  scad_a=3.7,
                  EPS=1e-9,
                  MAX_ITER=1e+5, 
                  solver='ISTA',
                  SVD = 1)
    plot(res1$v,type='l',xlab=paste("scad//l1//Non-nega l1","lambda=",l),xlim=c(0,70))
    res1 <- sfpca("PCA",
                  X=X,
                  Y=matrix(runif(p*n),n,p),
                  Omega_u=O_u,Omega_v=O_v,
                  alpha_u=1,alpha_v=1,
                  lambda_u=l,lambda_v=l,
                  P_u="scad",P_v="scad",
                  non_neg=0,
                  scad_a=3.7,
                  EPS=1e-9,
                  MAX_ITER=1e+5, 
                  solver='ISTA',
                  SVD = 1)
    lines(res1$v,col=cnt+2)
    res1 <- sfpca("PCA",
                  X=X,
                  Y=matrix(runif(p*n),n,p),
                  Omega_u=O_u,Omega_v=O_v,
                  alpha_u=1,alpha_v=1,
                  lambda_u=l,lambda_v=l,
                  P_u="L1",P_v="L1",
                  non_neg=1,
                  scad_a=3.7,
                  EPS=1e-9,
                  MAX_ITER=1e+5, 
                  solver='ISTA',
                  SVD = 1)
    lines(res1$v,col=cnt+1)
    
  }
