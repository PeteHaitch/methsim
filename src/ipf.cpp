// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

//' Two-dimensionsal Iterative Proportional Fitting.
//'
//' @description The function implements the iteratitive proportional fitting
//' (ipf) procedure for \strong{;two-dimensional tables only}. The code is
//' adapted from \code{mipfp::\link[Ipfp]{Ipfp}} for the special case of
//' two-dimensional tables without \code{NA} values.
//'
//' @param seed the initial two-dimensional array to be updated. Each cell must
//' be non-negative.
//' @param row_margins the contraints on the row margins.
//' @param col_margins the constraints on the column margins.
//' @param iter stopping criterion. The maximum number of iterations allowed;
//' must be greater than 0. Default is 1000.
//' @param tol stopping criterion. If the maximum absolute difference between
//' two iterations is lower than the value specified by \code{tol}, then
//' ipf has reached convergence; must be greater than 0. Default is 1e-10.
//'
//' @return \code{xi.hat} an \code{arma::mat} of the same dimension of
//' \code{seed} whose margins match the ones specified in \code{target.list}.
//'
//' @keywords internal
//'
//' @note Compared to \code{mipfp::\link[Ipfp]{Ipfp}}, this function only
//' returns "xi.hat" (as an arma::mat) and none of "stp.crit", "conv",
//' "check.margins", "evol.stp.crit". It is designed for internal use in the
//' methsim package by other C++ functions (although currently an R interface
//' is exposed for testing purposes).
//' This function has only been tested with 2x2 matrices.
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
arma::mat ipf(const arma::mat& seed, arma::vec row_margins,
           arma::rowvec col_margins, int iter = 1000,
           double tol = 1e-10) {

  // Variable initialisations.
  // Note that row_margins are a column vector and col_margins are a row vector.
  // Copy seed so that it is not modified.
  arma::mat result(seed.n_rows, seed.n_cols);
  result = seed;
  // result_temp stores the previous iteration of the algorithm.
  arma::mat result_temp(seed.n_rows, seed.n_cols);
  arma::vec row_sums(seed.n_rows);
  arma::rowvec col_sums(seed.n_cols);
  // update_factor needs to be arma::vec because is multiplied with rows/cols
  // of an arma::mat.
  arma::vec update_factor(2);
  double stp_crit;
  bool converged = false;
  // evol_stp_crit needs to be a std::vector and not arma::vec because
  // .push_back() is required.
  std::vector<double> evol_stp_crit;
  // diff_margins needs to be std::vector and not arma::vec (resp.
  // arma::rowvec) because the wrap method for arma::vec (resp. arma::rowvec)
  // is a matrix with 1 column (resp. row).
  std::vector<double> diff_margins(2);
  int nr = result.n_rows;
  int nc = result.n_cols;

  for (int i = 0; i < iter; i++ ) {
    result_temp = result;

    // Process rows
    row_sums = sum(result, 1);
    // TODO: Are these operations in-place?
    for (int r = 0; r < nr; r++) {
      if (row_margins[r] == 0 || row_sums[r] == 0) {
        update_factor[r] = 0;
      } else {
        update_factor[r] = row_margins[r] / row_sums[r];
      }
    }
    // This emulates R's sweep() function.
    result.each_col() %= update_factor;

    // Process columns
    col_sums = sum(result, 0);
    // TODO: Are these operations in-place?
    for (int c = 0; c < nc; c++) {
      if (col_margins[c] == 0 || col_sums[c] == 0) {
        update_factor[c] = 0;
      } else {
        update_factor[c] = col_margins[c] / col_sums[c];
      }
    }
    // This emulates R's sweep() function.
    result.each_row() %= update_factor.t();

    // Check for convergence.
    // Need max(max(X)) to find maximum element of a matrix X.
    evol_stp_crit.push_back(max(max(abs(result - result_temp))));
    if (evol_stp_crit[i] < tol) {
      converged = true;
      stp_crit = evol_stp_crit[i];
      break;
    }
  }

  if (!converged) {
    std::stringstream iter_string;
    iter_string << iter;
    std::string warning_msg = "IPFP did not converge after " +
      iter_string.str() + " iterations(s)!\n" +
      "This might be due to 0 cells in the seed, the maximum number of " +
      "iterations being too low, or the tolerance being too small.";
    Rcpp::warning(warning_msg);
  }

  diff_margins[0] = max(abs(row_margins - sum(result, 1)));
  diff_margins[1] = max(abs(col_margins - sum(result, 0)));

  return result;
}
