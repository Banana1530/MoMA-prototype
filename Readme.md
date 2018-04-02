# MoMA Prototype demo

Currently my prototype supports the followings extensions:

#### Model: PCA, CCA

#### Optimization algorithm: ISTA, FISTA

#### Penalty: Lasso, MCP, Non-negativity Lasso, User-provide smoothing penalty

Right now I am trying to separate the implementation of functions from `.h` file to `.cpp` file but it would fail to compile for no reason. Still working on it.

The following gif shows how three different sparse penalties, L1(black), SCAD(green) and non-negativity L1 (red) effect the recovered signal with varying lambda. This example is from Allen, 2013, Sparse and Functional Principal Component Analysis, where a sinusoidal signal from a combination of three is extracted.

Primitive C++ code is posted now. R code for the demo is attached below.

![](Demo.gif)

```R
for(l in seq(1,10,0.5)){ # Changing lambda
       l1 <- sfpca("PCA",
                  X=X,
                  Y=NULL,
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
        scad <- sfpca("PCA",
                  X=X,
                  Y=NULL,
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
       l1nonneg <- sfpca("PCA",
                  X=X,
                  Y=NULL,
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
    plot(l1$v,type='l',xlab=paste("scad//l1//Non-nega l1","lambda=",l),xlim=c(0,70))
    lines(scad$v,col=cnt+2)
    lines(l1nonneg$v,col=cnt+1)
    save.image()
  }
}
```

