library(Rcpp);
library(RcppArmadillo);
library(RcppEigen);
library(plyr);
library(fda);
library(fossil);
library(Matrix);
library(matrixcalc);
library(igraph);

set.seed(100)
sourceCpp("admm_related_functions.cpp");
source("R_functions.R")

###############Settings#################
# coefficient functions for 3 subgroups
Beta_set = c(Beta89,Betacubic,Betaexp)

rangeval <- c(0,1)
# m and d of B-spline
nknots = 8
norder = 4
nbasis <- nknots + norder - 2
p = nbasis
basisfd_beta <- create.bspline.basis(rangeval= rangeval, nbasis=nbasis, norder=norder)
# the sample size of each subgroup
n = 20
n1 = 20
n2 = 20
n3 = 20
Beta = c(rep(Beta_set[1],n1),rep(Beta_set[2],n2),rep(Beta_set[3],n3))
true_membership = c(rep(1,n1),rep(2,n2),rep(3,n3))

data = readRDS("sample_data.rds")
n = 3*n

# calculate some matrices
xcoef <- data$x$coefs
Jmat <- inprod(data$x$basis, basisfd_beta)
Uvec <- t(xcoef)%*%Jmat
Ulist = list()
for (i in 1:nrow(Uvec)){
  Ulist[[i]] = Uvec[i,]
}
U = bdiag(Ulist)
U = t(as.matrix(U))

nclass = n
AA <- matrix(0,nrow=nclass, ncol=nclass)
for(i in 1:nclass){
  for(j in 1:nclass){
    if(i==j){
      AA[i,j]=nclass-1
    }
    else{
      AA[i,j]=-1
    }
  }
}
K <- NULL
BB <- inprod(basisfd_beta, basisfd_beta) # M+d by M+d
K <- kronecker(AA,BB) # n(M+d) by n(M+d)

#compute 2-order diffrence matrix
C <- matrix(0, nrow=nknots+norder-2-2, ncol=nknots+norder-2)
for (j in 1:(nknots+norder-2-2)){
  d_j <- c(rep(0,j-1),1,-2,1,rep(0,(nknots+norder-2)-3-(j-1)))
  e_j <- c(rep(0,j-1), 1 ,rep(0,(nknots+norder-2)-3-(j-1)))
  C <- C + e_j%*%t(d_j)
}
D = t(C)%*%C;
diagD <- kronecker(diag(1, n), D);

####################### ADMM algorithm###########################
t1=proc.time();
####################### Choose gamma1,2 by GCV and BIC###########################
#Step1.Given gamma1=0.005, choose gamma2.
#Step2.Given optimal gamma2, choose gamma1.
gamma1_list = c(0.001,0.005,0.01);
gamma2_list = c(0.5,1,2);
index = t(combn(n,2));

#Step1:
gamma1 = 0.005;
best_bic = 10000;
for(i in gamma2_list){
  gamma2 = i;
  B_ini0 = matrix(rep(solve(t(Uvec)%*%Uvec+2*gamma1*D)%*%t(Uvec)%*%data$ycol,n),nrow=nbasis,ncol=n)
  sol_final = prclust_admm(U, y=as.vector(data$ycol), diagD, B_ini0,index,gamma1 = gamma1, gamma2 = gamma2,theta=0.005, tau = 2, n, p,  max_iter=1000,eps_abs=1e-4, eps_rel=1e-4)
  Ad_final <- create_adjacency(sol_final$V, n);
  G_final <- graph.adjacency(Ad_final, mode = 'upper')
  #clustering membership
  cls_final <- components(G_final);
  #number of clusters
  k_final <- cls_final$no;
  fdobj_all = fd(sol_final$B,basisfd_beta)
  estimatey = cbind(diag(inprod(data$x, fdobj_all)),data$ycol)
  bic <- cal_bic2_y(estimatey,k_final);
  if(bic<best_bic){
    best_gamma2 = gamma2;
    best_bic = bic;
  }
}

#Step2:
gamma2 = best_gamma2;
best_gcv = 10000;
for(i in gamma1_list){
  gamma1 = i;
  B_ini0 = matrix(rep(solve(t(Uvec)%*%Uvec+2*gamma1*D)%*%t(Uvec)%*%data$ycol,n),nrow=nbasis,ncol=n)
  sol_final = prclust_admm(U, y=as.vector(data$ycol), diagD, B_ini0,index,gamma1 = gamma1, gamma2 = gamma2,theta=0.005, tau = 2, n, p,  max_iter=1000,eps_abs=1e-4, eps_rel=1e-4)
  Ad_final <- create_adjacency(sol_final$V, n);
  G_final <- graph.adjacency(Ad_final, mode = 'upper')
  cls_final <- components(G_final);
  k_final <- cls_final$no;
  fdobj_all = fd(sol_final$B,basisfd_beta)
  estimatey = cbind(diag(inprod(data$x, fdobj_all)),data$ycol)
  gcv <- GCV(estimatey,k_final);
  if(gcv<best_gcv){
    best_gamma1 = i;
    best_gcv = gcv;
    best_sol = sol_final;
    best_cls_final = cls_final;
  }
}
sol_final = best_sol;
fdobj_all = fd(sol_final$B,basisfd_beta);

# mean estimated coefficient function of each subgroup 
fdobj = mean.fd(fdobj_all[groups(best_cls_final)$'1'])
fdobj2 = mean.fd(fdobj_all[groups(best_cls_final)$'2'])
fdobj3 = mean.fd(fdobj_all[groups(best_cls_final)$'3'])
tabl = cbind(inprod(data$x, fdobj),data$ycol)

# plot true and estimated coefficient functions
plot(fdobj,ylim = c(-6,4))
seqq=seq(from=0,to=1,length.out=100)
lines(x=seqq,y=Beta89(seqq),col="red")
lines(x=seqq,y=Betacubic(seqq),col="red")
lines(x=seqq,y=Betaexp(seqq),col="red")
lines(fdobj2)
lines(fdobj3)

# gamma2 chosen by BIC
best_gamma2
# gamma1 chosen by GCV
best_gamma1

t2=proc.time();
t2-t1