#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// calculate prox_L2
// [[Rcpp::export]]
arma::vec prox_L2(arma::vec x, double sigma){
    int n = x.n_elem;
    double lv;
    arma::vec px(n);
    
    lv = norm(x, 2);
    
    if (lv == 0.)
        px = x;
    else
        px = fmax(0., 1. - (sigma/lv))*x;
    
    return px;
}

// calculate the inverse 
// [[Rcpp::export]]
arma::mat update_invmat(arma::mat X, arma::mat diagD, arma::vec y,
                int n, double gamma1, arma::mat index, double theta,int p) {
    
    arma::mat eyemat_n = arma::eye<arma::mat>(n,n);
    arma::mat eyemat_p = arma::eye<arma::mat>(p,p);
    arma::mat onemat_n = arma::ones<arma::mat>(n,n);
    arma::mat AtA = arma::kron(n*eyemat_n - onemat_n,eyemat_p);
    
    arma::mat XtX = arma::trans(X)*X;
    arma::mat invmat = arma::inv(XtX + gamma1*diagD + theta*AtA) ;

    return invmat;
}

// Update B
// [[Rcpp::export]]
arma::mat update_B_ad(int m, int p, arma::mat eyemat_n, arma::mat eyemat_p, arma::mat onemat_n, arma::mat XtX, arma::mat invmat, arma::mat X, arma::mat diagD, arma::vec y, arma::mat V,
                   arma::mat Lambda, int n, double gamma1, arma::mat index, double theta) {
    
    arma::mat B_new(p, n);
    arma::mat V_tilde(p, m);
    int i, j;
    
    V_tilde = V + Lambda/theta;
    
    arma::vec lv = arma::zeros<arma::vec>(n*p);
    
    for (int l=0; l<m; l++) {
        i = index(l, 0);
        j = index(l, 1);
        arma::vec  e1 = arma::zeros<arma::vec>(n);
        arma::vec  e2 = arma::zeros<arma::vec>(n);
        e1(i-1) = 1;
        e2(j-1) = 1;
        
        arma::mat Alt = arma::kron(e1 - e2, eyemat_p);
        lv = lv + Alt*V_tilde.col(l);
    }
    
    arma::vec rhs = arma::trans(X)*y + theta*lv;
    arma::vec b_new = invmat*rhs;
    
    B_new = arma::reshape(b_new, p, n);
    return B_new;
}

// [[Rcpp::export]]
arma::mat update_B_ad_fast(int m, int p, arma::mat eyemat_n, arma::mat eyemat_p, arma::mat onemat_n, arma::mat XtX, arma::mat invmat, arma::mat X, arma::mat diagD, arma::vec y, arma::mat V,
                   arma::mat Lambda, int n, double gamma1, arma::mat index, double theta) {
    
    arma::mat B_new(p, n);
    arma::mat V_tilde(p, m);
    int i, j;
    
    V_tilde = V + Lambda/theta;
    
    arma::vec lv = arma::zeros<arma::vec>(n*p);
    
    for (int l=0; l<m; l++) {
        i = index(l, 0);
        j = index(l, 1);
        arma::mat zero1 = arma::zeros(p,n);
        zero1.col(i-1) = V_tilde.col(l);
        zero1.col(j-1) = -V_tilde.col(l);
		arma::vec each = arma::reshape(zero1, p*n, 1);
        lv = lv + each;
    }
    
    arma::vec rhs = arma::trans(X)*y + theta*lv;
    arma::vec b_new = invmat*rhs;
    
    B_new = arma::reshape(b_new, p, n);
    return B_new;
}

// Update V
// [[Rcpp::export]]
arma::mat update_V(int m, int p, arma::mat B, arma::mat Lambda,
                   double gamma2, double theta, double tau, arma::mat index) {
    
    arma::mat V(p, m);
    double sigma = gamma2/theta;
    
    int i, j;
    arma::vec x(p);
    double mcp = tau*theta/(tau*theta -1);
    
    
    for (int l=0; l<m; l++) {
        i = index(l, 0);
        j = index(l, 1);
        
        x = B.col(i-1) - B.col(j-1) - (1.0/(theta))*Lambda.col(l);
        
        if (arma::norm(x,2) >= tau*gamma2)
            V.col(l) = x;
        else
            V.col(l) = prox_L2(x, sigma)*mcp;
    }
    return V;
}


// Update Lambda
// [[Rcpp::export]]
arma::mat update_Lambda(int m, int p, arma::mat Lambda, arma::mat B, arma::mat V,
                        double theta, arma::mat index) {
    
    arma::mat Lambda_new(p, m);
    
    int i, j;
    arma::vec x(p);
    
    for (int l=0; l<m; l++) {
        i = index(l, 0);
        j = index(l, 1);
        
        x = V.col(l) - B.col(i-1) + B.col(j-1);
        
        Lambda_new.col(l) = Lambda.col(l) + theta*x;
    }
    return Lambda_new;
}

// Stopping criterion related
// [[Rcpp::export]]
double tolerance_primal(arma::mat B, arma::mat V, double eps_abs, double eps_rel, arma::mat index) {
    int m = index.n_rows;
    int p = V.n_rows;
    int i,j;
    arma::vec dbl(p);
    arma::vec dvl(p);
    double db = 0.0;
    double dv = 0.0;
    
    double output = sqrt(p*m)*(eps_abs);
    
    for (int l=0; l<m; l++) {
        i = index(l, 0);
        j = index(l, 1);
        dbl = pow(B.col(i-1) - B.col(j-1),2);
        db = db + sum(dbl);
        dvl = pow(V.col(l),2);
        dv = dv + sum(dvl);
    }
    
    db = sqrt(db);
    dv = sqrt(dv);
    output = output + (eps_rel)*fmax(db,dv);
    
    return output;
}

// [[Rcpp::export]]
double tolerance_dual(arma::mat eyemat_p, arma::mat Lambda, arma::mat index, int n, double eps_abs, double eps_rel) {
    
    int m = index.n_rows;
    int p = Lambda.n_rows;
    
    double output = sqrt(p*n)*(eps_abs);
    double s;
    arma::vec AtLambda = arma::zeros<arma::vec>(p*n);
    
    for (int l=0; l<m; l++) {
        int i = index(l, 0);
        int j = index(l, 1);
        arma::vec  e1 = arma::zeros<arma::vec>(n);
        arma::vec  e2 = arma::zeros<arma::vec>(n);
        e1(i-1) = 1;
        e2(j-1) = 1;
        
        arma::mat Alt = arma::kron(e1 - e2, eyemat_p);
        AtLambda = AtLambda + Alt*Lambda.col(l);
    }
    
    s = arma::norm(AtLambda,2);
    
    output = output + s*(eps_rel);
    
    return output;
}

// [[Rcpp::export]]
double tolerance_dual_fast(arma::mat eyemat_p, arma::mat Lambda, arma::mat index, int n, double eps_abs, double eps_rel) {
    
    int m = index.n_rows;
    int p = Lambda.n_rows;
    
    double output = sqrt(p*n)*(eps_abs);
    double s;
    arma::vec AtLambda = arma::zeros<arma::vec>(p*n);
    
    for (int l=0; l<m; l++) {
        int i = index(l, 0);
        int j = index(l, 1);
        
        arma::mat zero2 = arma::zeros(p,n);
        zero2.col(i-1) = Lambda.col(l);
        zero2.col(j-1) = -Lambda.col(l);
		arma::vec each = arma::reshape(zero2, p*n, 1);
        AtLambda = AtLambda + each;
    }
    
    s = arma::norm(AtLambda,2);
    
    output = output + s*(eps_rel);
    
    return output;
}

// This function computes the L2 norm of the primal residual.
// [[Rcpp::export]]
double residual_primal(arma::mat B, arma::mat V, arma::mat index) {
    int m = index.n_rows;
    int p = B.n_rows;
    
    arma::vec x(p);
    double r = 0.0;
    
    for (int l=0; l<m; l++) {
        int i = index(l, 0);
        int j = index(l, 1);
        
        x = pow(B.col(i-1) - B.col(j-1) - V.col(l),2);
        
        r = r + sum(x);
        
    }
    
    double residual = sqrt(r);
    
    return residual;
}

// This function computes the L2 norm of the dual residual.
// [[Rcpp::export]]
double residual_dual(arma::mat V, arma::mat V_old, arma::mat index, int n, double theta) {
    
    int m = index.n_rows;
    int p = V.n_rows;
    double s = 0.0;
    
    arma::mat diff = V - V_old;
    arma::vec v(p);
    
    for(int i=0;i<n;i++){
        arma::vec v1 = arma::zeros<arma::vec>(p);
        arma::vec v2 = arma::zeros<arma::vec>(p);
        for(int ix=0;ix<m;ix++){
            int l1 = index(ix,0);
            if(l1 == i + 1)
                v1 = v1 + diff.col(ix);
            else
                v1 = v1 + arma::zeros<arma::vec>(p);
            
            int l2 = index(ix,1);
            if(l2 == i + 1)
                v2 = v2 + diff.col(ix);
            else
                v2 = v2 + arma::zeros<arma::vec>(p);
        }
        v = pow(v1-v2,2);
        s = s + sum(v);
    }
    
    double residual = theta*sqrt(s);
    return residual;
}

// [[Rcpp::export]]
List prclust_admm(arma::mat X, arma::vec y, arma::mat diagD, arma::mat B_0, arma::mat index,
                  double gamma1, double gamma2, double theta, double tau, int n, int p, int max_iter,
                  double eps_abs, double eps_rel) {
    
    int m = index.n_rows;
    arma::mat V = arma::zeros<arma::mat>(p,m);
    arma::mat Lambda = arma::zeros<arma::mat>(p,m);
    double rp, rd;
    double tp, td;
    
    arma::mat B_new(p,n);
    arma::mat V_old(p,m);
    arma::mat Lambda_old(p,m);
    arma::mat invmat(n*p,n*p);
    arma::mat XtX = arma::trans(X)*X;
    arma::mat eyemat_n = arma::eye<arma::mat>(n,n);
    arma::mat eyemat_p = arma::eye<arma::mat>(p,p);
    arma::mat onemat_n = arma::ones<arma::mat>(n,n);
    
    arma::vec res_pri = arma::zeros<arma::vec>(max_iter);
    arma::vec res_dual = arma::zeros<arma::vec>(max_iter);

    V = update_V(m, p, B_0, Lambda, gamma2, theta, tau, index);
    invmat = update_invmat(X, diagD, y,  n, gamma1, index, theta,p);
    
    int iter = 1;
    for (int its=0; its< max_iter; its++) {
        V_old = V;
        Lambda_old = Lambda;
//        B_new = update_B_ad(m, p, eyemat_n, eyemat_p, onemat_n, XtX, invmat, X, diagD, y, V_old, Lambda_old, n, gamma1, index, theta);
		B_new = update_B_ad_fast(m, p, eyemat_n, eyemat_p, onemat_n, XtX, invmat, X, diagD, y, V_old, Lambda_old, n, gamma1, index, theta);
        V = update_V(m, p, B_new, Lambda_old, gamma2, theta, tau, index);
        Lambda = update_Lambda(m, p, Lambda_old, B_new, V, theta, index);
        
        rp = residual_primal(B_new, V, index);
        res_pri(its) = rp;
        
        tp = tolerance_primal(B_new, V, eps_abs, eps_rel, index);
        
        rd = residual_dual(V, V_old, index, n, theta);
        res_dual(its) = rd;
        
//        td = tolerance_dual(eyemat_p, Lambda, index, n, eps_abs, eps_rel);
		td = tolerance_dual_fast(eyemat_p, Lambda, index, n, eps_abs, eps_rel);
       
        if ((rp <= tp) & (rd <= td)) break;
        
        iter = iter + 1;
    }
    if (iter > max_iter)
        iter = max_iter;
    
    return  List::create(Named("B") = B_new,
                         Named("Lambda") = Lambda,
                         Named("V") = V,
                         Named("res_pri") = res_pri,
                         Named("res_dual") = res_dual);
}
