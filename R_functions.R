library(MASS)

# some examples of coefficient functions
Beta89 <- function(t) {
  beta89.t <- rep(0,length(t))
  for(i in 1:length(t))
  {
    beta89.t[i] = 2-6*t[i]
  }
  return(beta89.t)
}
Betacubic <- function(t) {
  betacubic.t <- rep(0,length(t))
  for(i in 1:length(t))
  {
    betacubic.t[i] = -3+10*t[i]*t[i]-10*t[i]+6*(t[i])^3;
  }
  return(betacubic.t)
}
Betaexp <- function(t) {
  betaexp.t <- rep(0,length(t))
  for(i in 1:length(t))
  {
    betaexp.t[i] = 1.5-exp(2*t[i])+3*sin(2*t[i]);
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


