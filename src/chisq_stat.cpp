#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double chisq_stat(NumericMatrix X, NumericVector theta, NumericMatrix Sigma) {
  typedef Map<MatrixXd> MapMatd;
  typedef Map<VectorXd> MapVecd;
  const MapMatd A(as<MapMatd>(X));
  const MapMatd B(as<MapMatd>(Sigma));
  const MapVecd C(as<MapVecd>(theta));
  double Q(((A * C).transpose() * (A * B * A.transpose()).inverse() * A * C)[0]);
  
  return Q;
}