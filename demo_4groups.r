library(Rcpp);
library(RcppArmadillo);
library(RcppEigen);
library(plyr);
library(fda);
library(microbenchmark);
library(fossil);
library(Matrix);
library(matrixcalc);
library(igraph);

set.seed(100)
sourceCpp("admm_related_functions.cpp");
source("R_functions.R")

###############Settings#################
# coefficient functions for 4 subgroups
Beta_set = c(Betaquad,Beta89,Betacubic,Betaexp)

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
n4 = 20
Beta = c(rep(Beta_set[1],n1),rep(Beta_set[2],n2),rep(Beta_set[3],n3),rep(Beta_set[4],n4))
true_membership = c(rep(1,n1),rep(2,n2),rep(3,n3),rep(4,n4))

###################Data Generation##########################
data = DataYBeta(n,n1,n2,n3,n4,rangeval, Betaquad,Beta89,Betacubic,Betaexp,"Norm")
n = 4*n

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
result_1 <- data.frame(gamma1=numeric(), k_est=numeric(),ARI=numeric(),AMSE=numeric());
result_2 <- data.frame(gamma2=numeric(), k_est=numeric(),ARI=numeric(),AMSE=numeric());

gamma1_list = c(0.005,0.1);
gamma2_list = c(0.5,3,10);
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
  result_now = data.frame(gamma2=gamma2,K_est=k_final,
                          ARI=adj.rand.index(cls_final$membership, true_membership),
                          AMSE=sum(error_list(Beta,fdobj_all))/n);
  result_2 = rbind(result_2,result_now);
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
  result_now = data.frame(gamma1=gamma1,K_est=k_final,
                          ARI=adj.rand.index(cls_final$membership, true_membership),
                          AMSE=sum(error_list(Beta,fdobj_all))/n);
  result_1 = rbind(result_1,result_now);
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
fdobj4 = mean.fd(fdobj_all[groups(best_cls_final)$'4'])
tabl = cbind(inprod(data$x, fdobj),data$ycol)

# plot true and estimated coefficient functions
plot(fdobj,ylim = c(-20,6))
seqq=seq(from=0,to=1,length.out=100)
lines(x=seqq,y=Betaquad(seqq),col="red")
lines(x=seqq,y=Beta89(seqq),col="red")
lines(x=seqq,y=Betacubic(seqq),col="red")
lines(x=seqq,y=Betaexp(seqq),col="red")
lines(fdobj2,ylim = c(-20,6))
lines(fdobj3,ylim = c(-20,6))
lines(fdobj4,ylim = c(-20,6))
legend(-0.01,-15,c("estimated coefficient","true coefficient"),
       col=c("black","red"),text.col=c("black","red"),lty=c(1,1),cex=0.8)

result_2
# gamma2 chosen by BIC
best_gamma2
result_1
# gamma1 chosen by GCV
best_gamma1

t2=proc.time();
t2-t1