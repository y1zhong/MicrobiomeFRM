#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <set>


using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

SEXP ugeeBinarySymReal(int n1,
                       int n2,
                       
                       int n1_n2,  //the # non-NA values OF 0 or 1
                       int effective_n,
                       
                       NumericMatrix ymat,
                       NumericMatrix xmat, // Main effect of Beta-div
                       
                       DataFrame cat,
                       List cts,
                       
                       double tol = 1e-3,
                       int maxiter = 50) {
  int n(ymat.nrow());
  //int Kw(unique(wgroup1).size());
  //int p(1 + 2 + ( Kw*(Kw-1)/2));  // total number of parameters
  int l(cat.length());
  int ctslen(cts.length());
  int p(1 + ctslen);
  int iter(0);
  
  double h(0);
  double Xb(0);
  double h_numerator(0);
  //double V(1);
  double err(1);
  double s(0);
  
  NumericVector Kp(l);
  NumericVector Kc(l);
  NumericVector group;
  for(int k = 0; k < l; k++){
    group = cat[k];
    Kp[k] = unique(group).size();
    Kc[k] = p;
    p = p + (Kp[k] - 1) * Kp[k] / 2;
  }
  
  VectorXd theta_00(VectorXd(p).setOnes());
  VectorXd theta(theta_00 * 0.2); // a smaller starting value
  VectorXd theta_new(theta * 1.0);
  VectorXd delta(VectorXd(p).setZero());
  //VectorXd V_inv(VectorXd(p).setOnes());
  VectorXd U(VectorXd(p).setZero());
  MatrixXd U_diff(MatrixXd(p, p).setZero());
  
  // Using Newton's method to solve ugee
  while (err >= tol) {
    U.setZero();
    U_diff.setZero();
    
    for (int i = 0; i < n-1; ++i) {
      for (int j = i+1; j < n; ++j) {
        //for (int j = 0; j < n; ++j) {
        
        s = ymat(i, j);
        if ( !Rcpp::internal::Rcpp_IsNA(s) ){
          delta.setZero();
          //LogicalVector signn((i-j) > 0) ;
          int sig(i - j);
          delta[0] = sig / abs(sig) * xmat(i, j);
          
          for(int k = 0; k < ctslen; k++){
            NumericMatrix zmattemp = cts(k);
            delta[k+1] = zmattemp(i,j);
          }
          
          for(int k = 0; k < l; k++){
            group = cat[k];
            if (group[i] != group[j]) {
              delta[p - 1 - k] = sig / abs(sig);
            }
          }
          
          // h = (delta.transpose() * theta)[0];
          //h_numerator = (delta.transpose() * theta.array().exp().matrix())[0];
          Xb = (delta.transpose() * theta)[0];
          h_numerator = exp(Xb);
          h = h_numerator/(1+h_numerator);
          
          // Vi
          //V = h*(1-h);
          
          U += delta * (ymat(i, j) - h);
          U_diff += h*(1-h) * delta * delta.transpose();
        }
      }
    } // end of double for() loops
    
    MatrixXd U_diff_pinv = U_diff.completeOrthogonalDecomposition().pseudoInverse();
    //theta_new = theta + V_inv * U_diff_pinv * U;
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
  MatrixXd B(U_diff / (n1_n2));
  
  ////////////////////////////////////////////////03/24 noon: estimates
  // Sigma_theta part:
  VectorXd v_i(VectorXd(p).setZero());
  MatrixXd Sigma_U(MatrixXd(p, p).setZero());
  
  //std::set<int> I1;
  //std::set<int> I2;
  //std::set<int> J1;
  //std::set<int> J2;
  
  // set some initial values:
  //for (int i = 0; i < n1; i++) I1.insert(i);    // setI1: 1, ..., n1
  //for (int i = n1; i < n; i++) I2.insert(i);    // setI2: n1+1, ..., n
  //J1 = I2;
  //J2 = I1;
  
  //print out the set element
  //std::copy(I1.begin(),
  //I1.end(),
  //std::ostream_iterator<int>(std::cout, " "));
  //std::copy(I2.begin(),
  //I2.end(),
  //std::ostream_iterator<int>(std::cout, " "));
  
  for (int i = 0; i < n; ++i) {
    v_i.setZero();
    
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      
      s = ymat(i, j);
      if ( !Rcpp::internal::Rcpp_IsNA(s) ){
        
        ///Main body
        delta.setZero();
        int sig(i - j);
        delta[0] = sig / abs(sig) * xmat(i, j);
        
        for(int k = 0; k < ctslen; k++){
          NumericMatrix zmattemp = cts(k);
          delta[k+1] = zmattemp(i,j);
        }
        
        for(int k = 0; k < l; k++){
          group = cat[k];
          if (group[i] != group[j]) {
            delta[p - 1 - k] = sig / abs(sig);
          }
        }
        
        
        Xb = (delta.transpose() * theta)[0];
        h_numerator = exp(Xb);
        h = h_numerator/(1+h_numerator);
        
        v_i += delta * (ymat(i, j) - h) / (effective_n - 1);
        
      }
      
      // if ( (I2.count(i) && J2.count(j)) ){
      //   
      //   ///Main body
      //   delta.setZero();
      //   int sig(i - j);
      //   delta[p-1] = sig / abs(sig) * xmat(i, j);
      //   
      //   Xb = (delta.transpose() * theta)[0];
      //   h_numerator = exp(Xb);
      //   h = h_numerator/(1+h_numerator);
      //   
      //   v_i += delta * (ymat(i, j) - h) / (effective_n - 1);
      //   
      // }
    }
    
    Sigma_U += 4 * v_i * v_i.adjoint() / (effective_n - 1);
  } // end of double for() loops
  
  MatrixXd B_pinv = B.completeOrthogonalDecomposition().pseudoInverse();
  MatrixXd Sigma_theta(B_pinv * Sigma_U * B_pinv  / effective_n);
  VectorXd sigma_theta(Sigma_theta.diagonal());
  
  return List::create(Named("theta") = wrap(theta),
                      Named("Sigma_theta") = wrap(Sigma_theta),
                      Named("sigma_theta") = wrap(sigma_theta));
}

