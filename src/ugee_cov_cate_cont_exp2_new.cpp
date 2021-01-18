#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP ugeecov_cate_cont_exp2_new(NumericMatrix dmatY,
                            DataFrame cat,
                            List cts,
                            double tol = 1e-3,
                            int maxiter = 25) {
  int n(dmatY.nrow());
  //int Kx(unique(group).size());
  //int Kw(unique(wgroup).size());
  int l(cat.length());
  int s(cts.length());
  int p(0);
  int iter(0);
  int ppos(0);
  
  double h(0);
  double Xb(0);
  double err(1);
  
  NumericVector Kp(l);
  NumericVector Kc(l);
  NumericVector group;
  for(int k = 0; k < l; k++){
    group = cat[k];
    Kp[k] = unique(group).size();
    Kc[k] = p;
    p = p + (Kp[k] + 1) * Kp[k] / 2;
  }
 // NumericVector Kc(cumsum(Kp)); // cumulative number of parameters
  //Rcout << "The unique size is " << Kp << std::endl;
  //Rcout << "The p is " << p << std::endl;
  //Rcout << "The Kc is " << Kc << std::endl;
  
  VectorXd theta_00(VectorXd(p).setOnes());
  VectorXd theta(theta_00 * 0.5); // a smaller starting value
  VectorXd theta_new(theta * 1.0);
  VectorXd delta(VectorXd(p).setZero());
  VectorXd D(VectorXd(p).setZero());
  VectorXd U(VectorXd(p).setZero());
  MatrixXd U_diff(MatrixXd(p, p).setZero());
  
  // Using Newton's method to solve ugee
  while (err >= tol) {
    U.setZero();
    U_diff.setZero();
    
    for (int i = 0; i < n-1; ++i) {
      for (int j = i+1; j < n; ++j) {
        delta.setZero();
        
        delta[0] = 1;
        
        // categorical
        for(int k = 0; k < l; k++){
          group = cat[k];
          ppos = Kc[k];
          if (group[i] < group[j]) {
            delta[ppos+Kp[k]+(2*Kp[k]-group[i])*(group[i]-1)/2+group[j]-group[i]-1-k] = 1;
          } else if (group[i] > group[j]) {
            delta[ppos+Kp[k]+(2*Kp[k]-group[j])*(group[j]-1)/2+group[i]-group[j]-1-k] = 1;
          } else if (group[i] != 1) {
            delta[ppos+group[i]-1-k] = 1;
            //Rcout << "The ppos+group[i]-1-k is " << ppos+group[i]-1-k << std::endl;
          }
        }
       
       
        
        // continous
        for(int k = 0; k < s; k++){
          NumericMatrix dmatX = cts(k);
          delta[p-1-k] = dmatX(i,j);
        }

        Xb = (delta.transpose() * theta)[0];
        //Rcout << "The Xb is " << Xb << std::endl;
        h = exp(Xb);
        
        U += delta * (dmatY(i, j) - h);
        U_diff += (h)* delta * delta.transpose();
      }
      //Rcout << "The p-1 is " << p-1 << std::endl;
      //Rcout << "The p-1 is " << dmatX(i,j) << std::endl;
      //Rcout << "The delta is " << delta << std::endl;
      //break;
    } // end of double for() loops
    //break;
    //Rcout << "The delta is " << delta << std::endl;
    MatrixXd U_diff_pinv = U_diff.completeOrthogonalDecomposition().pseudoInverse();
   // Rcout << "The U_diff_pinv is " << U_diff_pinv << std::endl;
    theta_new = theta + U_diff_pinv * U;
    err = (theta_new - theta).lpNorm<Infinity>();
    theta = theta_new * 1.0;
    ++iter;
    //Rcout << "The theta is " << theta << std::endl;
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
      
      delta[0] = 1;
      
      for(int k = 0; k < l; k++){
        group = cat[k];
        ppos = Kc[k];
        if (group[i] < group[j]) {
          delta[ppos+Kp[k]+(2*Kp[k]-group[i])*(group[i]-1)/2+group[j]-group[i]-1-k] = 1;
        } else if (group[i] > group[j]) {
          delta[ppos+Kp[k]+(2*Kp[k]-group[j])*(group[j]-1)/2+group[i]-group[j]-1-k] = 1;
        } else if (group[i] != 1) {
          delta[ppos+group[i]-1-k] = 1;
        }
      }
      
      
      for(int k = 0; k < s; k++){
        NumericMatrix dmatX = cts(k);
        delta[p-1-k] = dmatX(i,j);
      }
      
      Xb = (delta.transpose() * theta)[0];
      h = exp(Xb);
      
      v_i += delta * (dmatY(i, j) - h) / (n - 1);
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

