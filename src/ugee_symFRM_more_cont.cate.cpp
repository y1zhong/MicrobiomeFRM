#include <Rcpp.h>
#include <RcppEigen.h>


using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

SEXP ugeeSym_cont_cate_cov(NumericMatrix ymat,
             NumericMatrix xmat, // Main effect of Beta-div
            DataFrame cat,
            List cts,
          double tol = 1e-3,
          int maxiter = 25) {
  //Rcout << "Start "  << std::endl;
  int n(ymat.nrow());
  //int Kw(unique(wgroup1).size());
  //int p(1 + 2 + ( Kw*(Kw-1)/2));  // total number of parameters
  int l(cat.length());
  int s(cts.length());
  int p(1 + s);
  int iter(0);
  
  double h(0);
  double err(1);
  
  NumericVector Kp(l);
  NumericVector Kc(l);
  NumericVector group;
  for(int k = 0; k < l; k++){
    group = cat[k];
    Kp[k] = unique(group).size();
    Kc[k] = p;
    p = p + (Kp[k] - 1) * Kp[k] / 2;
  }
  
  VectorXd theta(VectorXd(p).setOnes());
  VectorXd theta_new(theta * 1.0);
  VectorXd delta(VectorXd(p).setZero());
  VectorXd U(VectorXd(p).setZero());
  MatrixXd U_diff(MatrixXd(p, p).setZero());
  
  // Using Newton's method to solve ugee
  while (err >= tol) {
    U.setZero();
    U_diff.setZero();
    
    for (int i = 0; i < n-1; ++i) {
      for (int j = i+1; j < n; ++j) {
        
        delta.setZero();
        //LogicalVector signn((i-j) > 0) ;
        int sig(i - j);
        delta[0] = sig / abs(sig) * xmat(i, j);
        //IntegerVector xx;
        //xx = sign(i - j);
        for(int k = 0; k < s; k++){
          NumericMatrix zmattemp = cts(k);
          delta[k+1] = zmattemp(i,j);
        }
        
        //delta[1] = zmat1(i,j);
        //delta[2] = zmat2(i,j);
        for(int k = 0; k < l; k++){
          group = cat[k];
          if (group[i] != group[j]) {
            delta[p - 1 - k] = sig / abs(sig);
          }
        }
        
        h = (delta.transpose() * theta)[0];
        U += delta * (ymat(i, j) - h);
        U_diff += delta * delta.transpose();
      }
    } // end of double for() loops
    
    MatrixXd U_diff_pinv = U_diff.completeOrthogonalDecomposition().pseudoInverse();
    theta_new = theta + U_diff_pinv * U;
    
    err = (theta_new - theta).lpNorm<Infinity>();
    theta = theta_new * 1.0;
    ++iter;
    
    if (iter > maxiter) {
      Rcerr << "Not converge!\n";
      break;
    }
  } // end of while()
  
  // Sandwich part:
  MatrixXd B(U_diff / (n * (n - 1) / 2));
           
  
  // Sigma_theta part:
  VectorXd v_i(VectorXd(p).setZero());
  MatrixXd Sigma_U(MatrixXd(p, p).setZero());
  
     
    for (int i = 0; i < n; ++i) {
       v_i.setZero();
       
       for (int j = 0; j < n; ++j) {
         if (j == i) continue;
         
      
      delta.setZero();
      int sig(i - j);
      delta[0] = sig / abs(sig) * xmat(i, j);
      
      for(int k = 0; k < s; k++){
        NumericMatrix zmattemp = cts(k);
        delta[k+1] = zmattemp(i,j);
      }
      
      //delta[1] = zmat1(i,j);
      //delta[2] = zmat2(i,j);
      for(int k = 0; k < l; k++){
        group = cat[k];
        if (group[i] != group[j]) {
          delta[p - 1 - k] = sig / abs(sig);
        }
      }
      
      h = (delta.transpose() * theta)[0];
      v_i += delta * (ymat(i, j) - h) / (n - 1);
    }
    
    Sigma_U += 4 * v_i * v_i.adjoint() / (n - 1);
  } // end of double for() loops
  
  MatrixXd B_pinv = B.completeOrthogonalDecomposition().pseudoInverse();
  MatrixXd Sigma_theta(B_pinv * Sigma_U * B_pinv  / n);
  VectorXd sigma_theta(Sigma_theta.diagonal());
  
  return List::create(Named("theta") = wrap(theta),
                      Named("Sigma_theta") = wrap(Sigma_theta),
                      Named("sigma_theta") = wrap(sigma_theta));
}

