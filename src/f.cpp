#include <methsim.h>

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> f(NumericVector psi, NumericVector b1, NumericVector b2) {
  int n = psi.size();
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  arma::rowvec col_margins(2);
  arma::colvec row_margins(2);
  arma::mat joint_prob_matrix(2, 2);
  std::vector<int> Z(n);
  for (int i = 0; i < n; i++) {
    ipf_seed(0, 0) = psi[i];
    row_margins[0] = 1 - b1[i];
    row_margins[1] = b1[i];
    col_margins[0] = 1 - b2[i];
    col_margins[1] = b2[i];
    joint_prob_matrix = methsim::ipf(ipf_seed, row_margins, col_margins);
  }
  return Z;
}

// Test whether its methsim::ipf() that is the problem
// [[Rcpp::export]]
std::vector<int> f2(NumericVector psi, NumericVector b1, NumericVector b2) {
  int n = psi.size();
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  arma::rowvec col_margins(2);
  arma::colvec row_margins(2);
  arma::mat joint_prob_matrix(2, 2);
  std::vector<int> Z(n);
  for (int i = 0; i < n; i++) {
    ipf_seed(0, 0) = psi[i];
    row_margins[0] = 1 - b1[i];
    row_margins[1] = b1[i];
    col_margins[0] = 1 - b2[i];
    col_margins[1] = b2[i];
//    joint_prob_matrix = methsim::ipf(ipf_seed, row_margins, col_margins);
  }
  return Z;
}

// Test whether it is std::vector<int> Z(n) that is the problem
// [[Rcpp::export]]
int f3(NumericVector psi, NumericVector b1, NumericVector b2) {
  int n = psi.size();
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  arma::rowvec col_margins(2);
  arma::colvec row_margins(2);
  arma::mat joint_prob_matrix(2, 2);
//  std::vector<int> Z(n);
  for (int i = 0; i < n; i++) {
    ipf_seed(0, 0) = psi[i];
    row_margins[0] = 1 - b1[i];
    row_margins[1] = b1[i];
    col_margins[0] = 1 - b2[i];
    col_margins[1] = b2[i];
    joint_prob_matrix = methsim::ipf(ipf_seed, row_margins, col_margins);
  }
  return 0;
}

// Test whether it is the general creation of temporary matrices that is the problem
// [[Rcpp::export]]
std::vector<int> f4(NumericVector psi, NumericVector b1, NumericVector b2) {
  int n = psi.size();
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  arma::rowvec col_margins(2);
  arma::colvec row_margins(2);
  arma::mat joint_prob_matrix(2, 2);
  std::vector<int> Z(n);
  for (int i = 0; i < n; i++) {
    ipf_seed(0, 0) = psi[i];
    row_margins[0] = 1 - b1[i];
    row_margins[1] = b1[i];
    col_margins[0] = 1 - b2[i];
    col_margins[1] = b2[i];
    joint_prob_matrix << 1 - b1[i] << b1[i] << arma::endr
                      << 1 - b2[i] << b2[i] << arma::endr;
  }
  return Z;
}

// Initialise arma::mat joint_prob_matrix for each call of methsim::ipf()
// [[Rcpp::export]]
std::vector<int> f5(NumericVector psi, NumericVector b1, NumericVector b2) {
  int n = psi.size();
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  arma::rowvec col_margins(2);
  arma::colvec row_margins(2);
  std::vector<int> Z(n);
  for (int i = 0; i < n; i++) {
    ipf_seed(0, 0) = psi[i];
    row_margins[0] = 1 - b1[i];
    row_margins[1] = b1[i];
    col_margins[0] = 1 - b2[i];
    col_margins[1] = b2[i];
    arma::mat joint_prob_matrix = methsim::ipf(ipf_seed, row_margins, col_margins);
  }
  return Z;
}

// Don't initialise joint_prob_matrix with a particular size.
// [[Rcpp::export]]
std::vector<int> f6(NumericVector psi, NumericVector b1, NumericVector b2) {
  int n = psi.size();
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  arma::rowvec col_margins(2);
  arma::colvec row_margins(2);
  arma::mat joint_prob_matrix;
  std::vector<int> Z(n);
  for (int i = 0; i < n; i++) {
    ipf_seed(0, 0) = psi[i];
    row_margins[0] = 1 - b1[i];
    row_margins[1] = b1[i];
    col_margins[0] = 1 - b2[i];
    col_margins[1] = b2[i];
    joint_prob_matrix = methsim::ipf(ipf_seed, row_margins, col_margins);
  }
  return Z;
}