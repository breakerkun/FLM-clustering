library(MASS)

# some examples of coefficient functions
Beta89 <- function(t) {
  beta89.t <- rep(0,length(t))
  for(i in 1:length(t))
  {
    beta89.t[i] = 2-10*t[i]
  }
  return(beta89.t)
}
Betaquad <- function(t) {
  betaquad.t <- rep(0,length(t))
  for(i in 1:length(t))
  {
    betaquad.t[i] = -8-20*t[i]*t[i]+10*t[i];
  }
  return(betaquad.t)
}
Betacubic <- function(t) {
  betacubic.t <- rep(0,length(t))
  for(i in 1:length(t))
  {
    betacubic.t[i] = -6+10*t[i]*t[i]-10*t[i]+6*(t[i])^3;
  }
  return(betacubic.t)
}
Betaexp <- function(t) {
  betaexp.t <- rep(0,length(t))
  for(i in 1:length(t))
  {
    betaexp.t[i] = 2-exp(2*t[i])+3*sin(2*t[i]);
  }
  return(betaexp.t)
}

# L2 distance between two functions
L2_dist_func <- function(Beta1,Beta2){
  flog <- function(t) {(Beta1(t)-Beta2(t))^2}
  in.prod <- integrate(flog,0,1)
  return(sqrt(in.prod$value))
}

####################### Calculate InnerProd #######################
InnerProd <- function(Betaf,basisfd,j) {
  # compute the <beta_j, B_j>, integral of beta_j and B_j. 
  # Betaf: beta function, i.e. Beta1, Beta2, Beta3, Beta4, predefined.
  # basisfd: basis function, jth column of basismatrix (eval.basis) (t rows and nbasis columns)
  
  rng <- getbasisrange(basisfd)
  knots <- basisfd$param
  knots <- c(rng[1], knots, rng[2])
  nbasis <- basisfd$nbasis
  norder <- nbasis - length(knots) + 2
  
  a <-rng[1]
  if(j-norder > 0) {a <- knots[j-norder+1]}
  
  b <- rng[2]
  if(j <= (nbasis-norder)) {b <- knots[j+1]}
  
  BFun <- function(t) {
    basismatrix <- eval.basis(t,basisfd) # 71 by 74 matrix, t rows and nbasis column
    basismatrix.j <- t(basismatrix[,j]) #get jth column of basismatrix
    return(basismatrix.j)
  }
  
  flog <- function(t) {Betaf(t)*BFun(t)}
  in.prod <- integrate(flog,a,b) 
  return(in.prod$value)
}

# generate simulation data
DataYBeta <- function(n,n1,n2,n3,n4, rangeval, Beta.i, Beta.ii,Beta.iii,Beta.iv,aMat_Fun){
  #m and d of X(t)
  nnknots <- 15
  nnorder <- 5
  nnbasis <- nnknots + nnorder - 2
  basisfdData <- create.bspline.basis(rangeval=rangeval, nbasis=nnbasis, norder=nnorder)
  
  #inner product of <beta#_j, B_j>
  G1 <- matrix(0,nrow=nnbasis, ncol=1)
  for (j in 1:nnbasis){
    G1[j,1] <- InnerProd(Beta.i,basisfdData,j)
  } 
  G2 <- matrix(0,nrow=nnbasis, ncol=1)
  for (j in 1:nnbasis){
    G2[j,1] <- InnerProd(Beta.ii,basisfdData,j)
  } 
  G3 <- matrix(0,nrow=nnbasis, ncol=1)
  for (j in 1:nnbasis){
    G3[j,1] <- InnerProd(Beta.iii,basisfdData,j)
  } 
  G4 <- matrix(0,nrow=nnbasis, ncol=1)
  for (j in 1:nnbasis){
    G4[j,1] <- InnerProd(Beta.iv,basisfdData,j)
  } 
  # aij, coefficient of basis functionB_j
  if (aMat_Fun == "Norm") 
  {
    aMat <- matrix(rnorm(n1*nnbasis,mean=2,sd=1),n1,nnbasis)
    aMat2 <- matrix(rnorm(n2*nnbasis,mean=2,sd=1),n2,nnbasis)
    aMat3 <- matrix(rnorm(n3*nnbasis,mean=2,sd=1),n3,nnbasis)
    aMat4 <- matrix(rnorm(n4*nnbasis,mean=2,sd=1),n4,nnbasis)
  }else if (aMat_Fun=="Unif")
  {
    aMat <- matrix(runif(n1*nnbasis,min=-2,max=2),n1,nnbasis)
    aMat2 <- matrix(runif(n2*nnbasis,min=-2,max=2),n2,nnbasis)
    aMat3 <- matrix(runif(n3*nnbasis,min=-2,max=2),n3,nnbasis)
    aMat4 <- matrix(runif(n4*nnbasis,min=-2,max=2),n4,nnbasis)
  }
  
  y1Tru <- aMat%*%G1
  y2Tru <- aMat2%*%G2
  y3Tru <- aMat3%*%G3
  y4Tru <- aMat4%*%G4
  
  y1 <- y1Tru + rnorm(n1,mean=0,sd=1)
  y2 <- y2Tru + rnorm(n2,mean=0,sd=1)
  y3 <- y3Tru + rnorm(n3,mean=0,sd=1)
  y4 <- y4Tru + rnorm(n4,mean=0,sd=1)
  
  ycol = rbind(y1,y2,y3,y4)
  
  a<- rbind(aMat,aMat2,aMat3,aMat4)
  x<- fd(coef=t(a), basisobj=basisfdData)
  data <- list(x=x,ycol=ycol,y1=y1,y2=y2,y1Tru=y1Tru,y2Tru=y2Tru)
  return(data)
}

knots_eq3 <- function(x, k, m){
  #external knots are on boundary
  #return boundary with internal knots only
  #used in bs or bsplineS
  c(min(x), seq(from=min(x), to=max(x), length.out=m+2)[-c(1,m+2)], max(x))
}

create_adjacency <- function(V,n) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0);
  index = t(combn(n,2));
  i <- index[connected_ix,1]
  j <- index[connected_ix,2]
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  A[(j-1)*n + i] <- 1
  return(A)
}

# calculate GCV
GCV <- function(estimatey,K_est){
  fenzi = sum((estimatey[,1]-estimatey[,2])^2);
  Ulist_K = list()
  trace = 0
  for(i in 1:k_final){
    Ulist_K[[i]] = matrix(Uvec[groups(cls_final)[[i]],],nrow=length(groups(cls_final)[[i]]));
    trace = trace + sum(diag(Ulist_K[[i]]%*%ginv(t(Ulist_K[[i]])%*%Ulist_K[[i]]
                                                 +2*gamma1*D)%*%t(Ulist_K[[i]])))
  }
  fenmu = (1-trace/n)^2
  gcv = fenzi/fenmu
  return(gcv)
}

# BIC
cal_bic2_y <- function(estimatey,K_est){
  comp1 = log(sum((estimatey[,1]-estimatey[,2])^2)/n);
  dfsum = K_est*p*log(log(n+p));
  comp2 = (log(n)*dfsum)/n;
  bic2 = comp1 + comp2;
  return(bic2)
}

# calculate errors
error_inprod <- function(Beta,fdobj){
  fd_func <- function(t){
    eval.fd(fdobj,t)
  }
  flog <- function(t) {(Beta(t)-fd_func(t))^2}
  in.prod <- integrate(flog,0,1)
  return(in.prod$value)
}
error_list <- function(Beta,fdobj_all){
  fdlist = list();
  for(i in 1:n){
    fdlist[[i]] = fdobj_all[i];
  }
  return(as.numeric(Map(
    function(fn, value)
    {
      error_inprod(fn,value)
    },
    Beta,
    fdlist
  )))
}


