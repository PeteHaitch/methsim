// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// A C++ only function that creates a 2x2 matrix
arma::mat ipf2(const arma::mat& seed, 
               const arma::colvec& row_margins,
               const arma::rowvec& col_margins, 
               int iter = 1000,
               double tol = 1e-10) {
  // Variable initialisations.
  // Copy seed since it can't be modified.
  int nr = seed.n_rows;
  int nc = seed.n_cols;
  if (nr != nc) {
    Rcpp::stop("Currently only works with square matrices.");
  }
  arma::mat result(nr, nc);
  result = seed;
  // result_temp stores the previous iteration of the algorithm.
  arma::mat result_temp(nr, nc);
  arma::colvec row_sums(nr);
  arma::rowvec col_sums(nc);
  // update_factor needs to be arma::colvec, and not a std::vec, because is 
  // multiplied with rows/cols of an arma::mat.
  // NOTE: There would need to be two different `update_factor`s if nr != nc.
  arma::colvec update_factor(nr);
  double stp_crit;
  bool converged = false;
  // evol_stp_crit needs to be a std::vector and not arma::colvec because
  // .push_back() is required.
  std::vector<double> evol_stp_crit;
  // NOTE: diff_margins not used in this implementation. If wanting to use, 
  // then it needs to be std::vector and not arma::colvec (resp. arma::rowvec) 
  // because the wrap method for arma::colvec (resp. arma::rowvec) is a matrix 
  // with 1 column (resp. row).
//  std::vector<double> diff_margins(2);
  
  for (int i = 0; i < iter; i++) {
    // Store the result from the previous iteration of the algorithm
    result_temp = result;

    // Process rows
    row_sums = sum(result, 1);
    // TODO: Are these operations in-place?
    for (int r = 0; r < nr; r++) {
      if (row_margins(r) == 0 or row_sums(r) == 0) {
        update_factor(r) = 0;
      } else {
        update_factor(r) = row_margins(r) / row_sums(r);
      }
    }
    // This emulates R's sweep() function.
    result.each_col() %= update_factor;

    // Process columns
    col_sums = sum(result, 0);
    // TODO: Are these operations in-place?
    for (int c = 0; c < nc; c++) {
      if (col_margins(c) == 0 or col_sums(c) == 0) {
        update_factor(c) = 0;
      } else {
        update_factor(c) = col_margins(c) / col_sums(c);
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

  // NOTE: diff_margins isn't used in this implementation (see NOTE above)
//  diff_margins[0] = max(abs(row_margins - sum(result, 1)));
//  diff_margins[1] = max(abs(col_margins - sum(result, 0)));
  return result;
}


// [[Rcpp::export]]
std::vector<int> g(NumericVector psi, NumericVector b1, NumericVector b2) {
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
    joint_prob_matrix = ipf2(ipf_seed, row_margins, col_margins);
  }
  return Z;
}

// Test whether its ipf2() that is the problem
// [[Rcpp::export]]
std::vector<int> g2(NumericVector psi, NumericVector b1, NumericVector b2) {
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
//    joint_prob_matrix = ipf2(ipf_seed, row_margins, col_margins);
  }
  return Z;
}

// Test whether it is std::vector<int> Z(n) that is the problem
// [[Rcpp::export]]
int g3(NumericVector psi, NumericVector b1, NumericVector b2) {
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
    joint_prob_matrix = ipf2(ipf_seed, row_margins, col_margins);
  }
  return 0;
}

// Test whether it is the general creation of temporary matrices that is the problem
// [[Rcpp::export]]
std::vector<int> g4(NumericVector psi, NumericVector b1, NumericVector b2) {
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

// Initialise arma::mat joint_prob_matrix for each call of ipf2()
// [[Rcpp::export]]
std::vector<int> g5(NumericVector psi, NumericVector b1, NumericVector b2) {
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
    arma::mat joint_prob_matrix = ipf2(ipf_seed, row_margins, col_margins);
  }
  return Z;
}

// Don't initialise joint_prob_matrix with a particular size.
// [[Rcpp::export]]
std::vector<int> g6(NumericVector psi, NumericVector b1, NumericVector b2) {
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
    joint_prob_matrix = ipf2(ipf_seed, row_margins, col_margins);
  }
  return Z;
}

// Uses ipf2() rather than ipf()
// [[Rcpp::export(".simulateZ3")]]
std::vector<int> simulateZ3(NumericVector beta_by_region,
                        NumericVector lor_by_pair,
                        CharacterVector seqnames_one_tuples,
                        NumericVector u) {

  // Argument checks
  if (beta_by_region.size() != seqnames_one_tuples.size()) {
    stop("length(beta_by_region) != length(seqnames_one_tuples)");
  }
  if (beta_by_region.size() != u.size()) {
    stop("length(beta_by_region) != length(u)");
  }
  // There is only a value in lor_by_pair for pairs of methylation loci on the
  // same chromosome.
  if (lor_by_pair.size() !=
      (beta_by_region.size() - unique(seqnames_one_tuples).size())) {
    std::string stop_msg = "length(lor_by_pair) != ";
    stop_msg = stop_msg +
      "(length(beta_by_region) - length(unique(seqnames_one_tuples)))";
    stop(stop_msg);
  }

  // Initialise variables
  // TODO: n is a variable at runtime and this might be the cause of sefaults
  // (see http://stackoverflow.com/questions/17105555/rcpp-segfault-on-arrays-698152-if-integervector-is-declared)
  int n = beta_by_region.size();
  // Z stores the result.
  // TODO: This fills Z with 0; this might allow a speed-up where I can "skip" 
  // Z[i] that were simulated as 0.
  std::vector<int> Z(n);
  // ipf_seed is used to initialise ipf algorithm to get joint_prob_matrix.
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  // col_margins = (p_{0.}, p_{1.})
  arma::rowvec col_margins(2);
  // row_margins = (p_{.0}, p_{.1})
  arma::colvec row_margins(2);
  // The 2x2 matrix of joint probabilities (*not* the transition matrix).
  // joint_prob_matrix = | p_{0, 0}, p_{1, 0} |
  //                     | p_{1, 0}, p_{1, 1} |,
  // where p_{a, b} = Pr(Z_{i - 1} = a, Z_{i} = b)
  arma::mat joint_prob_matrix(2, 2);
  // p = Pr(Z_{i + 1} = 1 | Z_{i} = z_{i})
  double p;
  // j indexes the lor_by_pair vector.
  int j = 0;

  // Initialise the process (i = 0) by sampling from the marginal distribution.
  if (u[0] > beta_by_region[0]) {
    Z[0] = 0;
  } else {
    Z[0] = 1;
  }

  // Simulate the rest of the process.
  // TODO: Use iterators for i and j.
  for (int i = 1; i < n; i++) {

    // Check that the current methylation loci and the next are on the same
    // chromosome. If not, then simulate from the marginal distribution since
    // there is no lor_by_pair value
    if (seqnames_one_tuples[i] != seqnames_one_tuples[i - 1]) {
      if (u[i] > beta_by_region[i]) {
        Z[i] = 0;
      } else {
        Z[i] = 1;
      }
      // Don't increment j. There is only a value in lor_by_pair for pairs of
      // methylation loci on the same chromosome so when a pair is on different
      // chromosome we don't increment j.
      continue;
    } else {

      // Compute transition probability from beta_by_region and lor_by_pair.
      // Can get joint probability matrix by running iterative proportional
      // fitting on matrix(c(2 ^ (lor), 1, 1, 1), ncol = 2) with marginals
      // given by the average methylation level of the region for the first and
      // second methylation loci, respectively.
      // Then, compute transition probabilities using
      // Pr(Z_{i + 1} = z_{i + 1} | Z_i = z_i) =
        // Pr(Z_i = z_i, Z_{i + 1} = z_{i + 1}) /
        // Pr(Z_i = z_i).

      // NOTE: This assumes lor_by_pair uses base-2 logarithms.
      ipf_seed(0, 0) = pow(2.0, lor_by_pair[j]);
      // col_margins refer to the (i - 1)-th methylation loci
      col_margins[0] = 1 - beta_by_region[i - 1];
      col_margins[1] = beta_by_region[i - 1];
      // row_margins refer to the i-th methylation loci
      row_margins[0] = 1 - beta_by_region[i];
      row_margins[1] = beta_by_region[i];
      // ipf2, which is in the same source file as simulateZ3, works, whereas 
      // ipf, which is in a separate source file to simulateZ2 and is exported 
      // using Rcpp::interfaces(r, cpp), does not work.
      joint_prob_matrix = ipf2(ipf_seed, row_margins, col_margins,
                               1000, 1e-10);
      if (Z[i - 1] == 0) {
        p = joint_prob_matrix(0, 1) / (1 - beta_by_region[i - 1]);
      } else {
        p = joint_prob_matrix(1, 1) / beta_by_region[i - 1];
      }
      if (u[i] > p) {
        Z[i] = 0;
      } else {
        Z[i] = 1;
      }

      // Increment j.
      j += 1;
    }
  }
  return Z;
}
