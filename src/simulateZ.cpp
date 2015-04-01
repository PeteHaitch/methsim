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

#include <methsim.h>

using namespace Rcpp;

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
      // TODO: May need to move ipf source code to this file and make it a
      // C++ only function.
      joint_prob_matrix = methsim::ipf(ipf_seed, row_margins, col_margins,
                                       1000, 1e-10);
      // Compute p = Pr(Z_{i + 1} = 1 | Z_{i} = z_{i})
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

