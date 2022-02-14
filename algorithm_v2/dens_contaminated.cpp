#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
//[[Rcpp::export]]
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
double dmvnorm_distance(arma::rowvec y, arma::mat Sigma);
double dmvnorm_rcpp(arma::rowvec y, arma::mat Sigma);
double  dens_contaminated(vec data, double beta_old, vec mu_P, double k_P, int nu_P, mat S_P, int p, int T);

// [[Rcpp::export]]
double  dens_contaminated(vec data, double beta_old, vec mu_P, double k_P, int nu_P, mat S_P, int p, int T) {
  double weight = 0;
  double log_weight = 0;
  double eval = 0;
  
  log_weight = log(1-beta_old);
  for (int i=0; i < T; i++){
    // Sample mu - TRONCATA 
    // 
    // Sample sigma
    // 
    // Eval K density of a normal - LOG
    eval += BOH ;
  }
  eval = eval/T; 
  
  log_weight += eval;
  weight = exp(log_weight);
  
  return weight;
}

// Evaluate a multivariate normal density multivariate
double dmvnorm_rcpp(arma::rowvec y, arma::mat Sigma) {
  
  int p = Sigma.n_rows;
  
  // inverse Sigma
  arma::mat Sigma1 = arma::inv(Sigma);
  // determinant Sigma
  double det_Sigma = arma::det(Sigma);
  // distance
  double dist = dmvnorm_distance( y, Sigma1);
  
  double pi1 = 3.14159265358979;
  double l1 = - p * std::log(2*pi1) - dist - std::log( det_Sigma );
  double ll = 0.5 * l1;
  
  return ll;
}

double dmvnorm_distance(arma::rowvec y, arma::mat Sigma)
{
  int n = Sigma.n_rows;
  double res=0;
  double fac=1;
  for (int ii=0; ii<n; ii++){
    for (int jj=ii; jj<n; jj++){
      if (ii==jj){ fac = 1; } else { fac = 2;}
      res += fac *y(0,ii) * Sigma(ii,jj) * y(0,jj);
    }
  }
  return res;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
// 
// HOW TO RUN MY FUNCTION
// library(Rcpp)
// sourceCpp('dens_contaminated.cpp')
//

