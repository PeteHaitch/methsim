// Copyright (C) 2015 Peter Hickey
//
// This file is part of methsim.
//
// methsim is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// methsim is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with methsim  If not, see <http://www.gnu.org/licenses/>.

// The ipf function in this file is adapted from the Ipfp function available in
// the mipfp CRAN package version 2.0.
// (http://cran.r-project.org/web/packages/mipfp/index.html).
// mipfp is released under a GPL (>= 2) license.

// TODO (long term): Figure out why ipf() only works when included int he same
// source file as simulateZ(), i.e., why doesn't Rcpp::interfaces(r, cpp) +
// Makevars work as expected? Create a simple example that illustrates this
// problem. The bug seems to be similar to
// http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-May/005775.html

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// Two-dimensionsal Iterative Proportional Fitting.
//
// @description The function implements the iteratitive proportional fitting
// (ipf) procedure for \strong{two-dimensional square tables without
// \code{NA} values}. The code is adapted from \code{mipfp::\link[Ipfp]{Ipfp}}.
//
// @param seed the initial two-dimensional array to be updated. Each cell must
// be non-negative.
// @param row_margins the contraints on the row margins.
// @param col_margins the constraints on the column margins.
// @param iter stopping criterion. The maximum number of iterations allowed;
// must be greater than 0. Default is 1000.
// @param tol stopping criterion. If the maximum absolute difference between
// two iterations is lower than the value specified by \code{tol}, then
// ipf has reached convergence; must be greater than 0. Default is 1e-10.
//
// @return \code{xi.hat} an \code{arma::mat} of the same dimension of
// \code{seed} whose margins match the ones specified in \code{target.list}.
//
// @keywords internal
//
// @note Compared to \code{mipfp::\link[Ipfp]{Ipfp}}, this function only
// returns "xi.hat" (as an arma::mat) and none of "stp.crit", "conv",
// "check.margins", "evol.stp.crit". It is designed for internal use in the
// methsim package by other C++ functions and has no direct R interface.
// This function has only been tested with 2x2 matrices.
arma::mat ipf(const arma::mat& seed,
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

// TODO: Is it better to use std::vector<int> or Rcpp::IntegerVector
// TODO (long term): Allow beta_by_region and seqnames_one_tuple to be
// S4Vectors::Rleobjects to avoid the requirement to pre-expand these via
// as.vector().
//' Simulate a single "haplotype" of a methylome (Z).
//'
//' @param beta_by_region the beta-value (average methylation level) for each
//' methylation locus in the genome.
//' @param lor_by_pair the within-fragment co-methylation between each pair of
//' methylation loci in the genome. Should be log odds-ratios using base-2
//' logarithms. The length of this should be equal to the number of methylation
//' loci in the genome minus the number of chromosomes (seqnames).
//' @param seqnames_one_tuples the chromosome (seqname) of each methylation
//' locus in the genome, i.e., \code{seqnames(one_tuples)}.
//' @param u a vector of Uniform(0, 1) random variables used in choosing the
//' next state of the process. Must have same length as \code{beta_by_region}.
//'
//' @note Random number generation is performed outside of this function (at
//' the R level in a single thread) in order to simplify this function. If
//' the vector of Uniform(0, 1) random variables ('u') is to be simulated
//' within this function at the C++ level then great care must be taken,
//' especially if this function is subsequently called in parallel.
//'
//' @return an integer vector of simulated methylation states for each
//' methylation locus in the genome; 0 = unmethylated and 1 = methylated.
// [[Rcpp::export(".simulateZ")]]
std::vector<int> simulateZ(NumericVector beta_by_region,
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
  int n = beta_by_region.size();
  // Z stores the result.
  // TODO: This fills Z with 0; this might allow a speed-up where I can "skip"
  // Z[i] that were simulated as 0. Alternatively, since most loci will be 1,
  // might fill with ones and then flip to zeroes as appropriate.
  std::vector<int> Z(n);
  // ipf_seed is used to initialise ipf algorithm to get joint_prob_matrix.
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  // col_margins = (p_{0.}, p_{1.})
  arma::rowvec col_margins(2);
  // row_margins = (p_{.0}, p_{.1})
  arma::colvec row_margins(2);
  // The 2x2 matrix of joint probabilities (*not* the transition matrix).
  // joint_prob_matrix = | p_{0, 0}, p_{0, 1} |
  //                     | p_{1, 0}, p_{1, 1} |,
  // where p_{a, b} = Pr(Z_{i - 1} = a, Z_{i} = b), i.e., rows refer to the
  // (i - 1)-th methylation locus and columns refer to the i-th methylation
  // locus.
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
      //    Pr(Z_{i + 1} = z_{i + 1} | Z_i = z_i)
      //  = Pr(Z_i = z_i, Z_{i + 1} = z_{i + 1}) /  Pr(Z_i = z_i).

      // NOTE: This assumes lor_by_pair uses base-2 logarithms.
      ipf_seed(0, 0) = pow(2.0, lor_by_pair[j]);
      // row_margins refer to the (i - 1)-th methylation locus.
      row_margins[0] = 1 - beta_by_region[i - 1];
      row_margins[1] = beta_by_region[i - 1];
      // col_margins refer to the i-th methylation locus.
      col_margins[0] = 1 - beta_by_region[i];
      col_margins[1] = beta_by_region[i];
      joint_prob_matrix = ipf(ipf_seed,
                              row_margins,
                              col_margins,
                              1000,
                              1e-10);
      // Compute p = Pr(Z_{i + 1} = 1 | Z_{i} = z_{i})
      if (Z[i - 1] == 0) {
        p = joint_prob_matrix(0, 1) / (1 - beta_by_region[i - 1]);
      } else {
        p = joint_prob_matrix(1, 1) / beta_by_region[i - 1];
      }
      // Sanity check the resulting p; it should be a probability!
      if (p <= 0 or p >= 1) {
        Rcpp::Rcout << "i = " << i << std::endl;
        Rcpp::Rcout << "beta_by_region[i - 1] = " << beta_by_region[i - 1] <<
          std::endl;
        Rcpp::Rcout << "joint_prob_matrix = " << joint_prob_matrix << std::endl;
        Rcpp::Rcout << "p = " << p << std::endl;
        Rcpp::stop("Unexpected p; it should be a probability!");
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

// TODO: Switch beta_by_region = marginalProb and lor_by_pair = LOR.
// TODO (long term): Allow beta_by_region and seqnames_one_tuple to be
// S4Vectors::Rleobjects to avoid the requirement to pre-expand these via
// as.vector().
//' Compute the transition probabilities, P, of a first-order binary Markov
//' chain given a set of marginal probabilities and log odds ratios.
//'
//' @param marginalProb the marginal probability (i.e., average methylation
//' level) for each methylation locus in the genome.
//' @param LOR the within-fragment co-methylation between each pair of
//' methylation loci in the genome. Should be log odds-ratios using base-2
//' logarithms. The LOR is 0 for the first methylation locus on each seqlevel.
//' @param seqnames_one_tuples the chromosome (seqname) of each methylation
//' locus in the genome, i.e., \code{seqnames(one_tuples)}.
//'
//' @return A $n \times 2$ matrix of the transition probabilities. Let
//' $P_{i, .}$ be the $i$-th row of P ($i = 0, \ldots, n - 1$), then
//' $P_{i, .} = [Pr(Z_{i + 1} = 1 | Z_{i} = 0), Pr(Z_{i + 1} = 1 | Z_{i} = 1)]$;
//' NB: $P_{0, } = [Pr(Z_{0} = 1), Pr(Z_{0} = 1)], i.e., sampled from the
//' marginal distribution and similarly for all other $i$ that start a new
//' chromosome.
//' More generally, we might use a scheme where the $j$-th column of P
//' ($j = 1, \ldots, 2^{m - 1}$, with $m = 2$ under the default) stores
//' $Pr(Z_{i} = 1 | Z_{i - 1}, Z_{i - 2}, ...)$.
// [[Rcpp::export(".computeP")]]
Rcpp::NumericMatrix computeP(NumericVector marginalProb,
                             NumericVector LOR,
                             int mc_order = 1) {
  // Argument checks
  if (marginalProb.size() != LOR.size()) {
    stop("length(marginalProb) != length(LOR)");
  }

  if (mc_order != 1) {
    stop("Only mc_order = 1 is currently supported.");
  }

  // Initialise variables
  int n = marginalProb.size();
  // P stores the result.
  Rcpp::NumericMatrix P(n, pow(2.0, mc_order));
  // ipf_seed is used to initialise ipf algorithm to get joint_prob_matrix.
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  // col_margins = (p_{0.}, p_{1.})
  arma::rowvec col_margins(2);
  // row_margins = (p_{.0}, p_{.1})
  arma::colvec row_margins(2);
  // The 2x2 matrix of joint probabilities (*not* the transition matrix).
  // joint_prob_matrix = | p_{0, 0}, p_{0, 1} |
  //                     | p_{1, 0}, p_{1, 1} |,
  // where p_{a, b} = Pr(Z_{i - 1} = a, Z_{i} = b), i.e., rows refer to the
  // (i - 1)-th methylation locus and columns refer to the i-th methylation
  // locus.
  arma::mat joint_prob_matrix(2, 2);

  // Compute the rest of the transition probabilities
  for (int i = 0; i < n; i++) {

    // Compute transition probability from marginalProb and LOR.
    // Can get joint probability matrix by running iterative proportional
    // fitting on matrix(c(2 ^ (lor), 1, 1, 1), ncol = 2) with marginals
    // given by the average methylation level of the region for the first and
    // second methylation loci, respectively.
    // Then, compute transition probabilities using
    //    Pr(Z_{i + 1} = z_{i + 1} | Z_i = z_i)
    //  = Pr(Z_i = z_i, Z_{i + 1} = z_{i + 1}) /  Pr(Z_i = z_i).

    // NOTE: This assumes LOR uses base-2 logarithms.
    // ipf_seed(0, 0) = pow(2.0, LOR[j]);
    ipf_seed(0, 0) = pow(2.0, LOR[i]);
    // row_margins refer to the (i - 1)-th methylation locus.
    row_margins[0] = 1 - marginalProb[i - 1];
    row_margins[1] = marginalProb[i - 1];
    // col_margins refer to the i-th methylation locus.
    col_margins[0] = 1 - marginalProb[i];
    col_margins[1] = marginalProb[i];
    joint_prob_matrix = ipf(ipf_seed,
                            row_margins,
                            col_margins,
                            1000,
                            1e-10);
    P(i, 0) = joint_prob_matrix(0, 1) / (1 - marginalProb[i - 1]);
    P(i, 1) = joint_prob_matrix(1, 1) / marginalProb[i - 1];
  }

  // These column names are really only set so that P can be an assay in a
  // SummarizedExperiment-based object.
  // TODO (long-term): If mc_order > 1 then I want a strategy so that I can
  // quickly get a column based on the previous states of the chain,
  // e.g., column 5 = 1 * 2^2 + 0 * 2^1 + 1 * 2^0 for the pattern
  // (M, U, M) = (1, 0, 1).Any such scheme needs to have an order for the
  // history, i.e., left-to-right or right-to-left.
  P.attr("dimnames") = List::create(R_NilValue, seq_len(pow(2.0, mc_order)));
  return P;
}
