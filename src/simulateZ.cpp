#include <methsim.h>

using namespace Rcpp;

// TODO: I'm still getting occasional low-level errors, e.g., segfaults. These
// seem to occur when .simulateZ() is run in parallel, particularly via
// BiocParallel::bplapply(). It would be good to get to the bottom of these;
// see http://master.bioconductor.org/developers/how-to/c-debugging/,
// https://www.google.com.au/search?q=rcpp+gdb+site%3Alists.r-forge.r-project.org&oq=rcpp+gdb+site%3Alists.r-forge.r-project.org&aqs=chrome..69i57.388j0j9&sourceid=chrome&es_sm=91&ie=UTF-8,
//
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
IntegerVector simulateZ(NumericVector beta_by_region,
                        NumericVector lor_by_pair,
                        CharacterVector seqnames_one_tuples,
                        NumericVector u) {

  // Argument checks
  if (beta_by_region.length() != seqnames_one_tuples.length()) {
    stop("length(beta_by_region) != length(seqnames_one_tuples)");
  }
  if (beta_by_region.length() != u.length()) {
    stop("length(beta_by_region) != length(u)");
  }
  // There is only a value in lor_by_pair for pairs of methylation loci on the
  // same chromosome.
  if (lor_by_pair.length() !=
      (beta_by_region.length() - unique(seqnames_one_tuples).length())) {
    std::string stop_msg = "length(lor_by_pair) != ";
    stop_msg = stop_msg +
      "(length(beta_by_region) - length(unique(seqnames_one_tuples)))";
    stop(stop_msg);
  }

  // Initialise variables
  // TODO: n is a variable at runtime and this might be the cause of sefaults
  // (see http://stackoverflow.com/questions/17105555/rcpp-segfault-on-arrays-698152-if-integervector-is-declared)
  int n = beta_by_region.length();
  // Z stores the result.
  IntegerVector Z(n, NA_INTEGER);
  // ipf_seed is used to initialise ipf algorithm to get joint_prob_matrix.
  arma::mat ipf_seed(2, 2, arma::fill::ones);
  // col_margins = (p_{0.}, p_{1.e})
  arma::rowvec col_margins(2);
  // row_margins = (p_{.0}, p_{.1})
  arma::vec row_margins(2);
  // The 2x2 matrix of joint probabilities (*not* the transition matrix).
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
    }

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
    ipf_seed(0, 0) = pow(2, lor_by_pair[i - 1]);
    col_margins[0] = 1 - beta_by_region[i - 1];
    col_margins[1] = beta_by_region[i - 1];
    row_margins[0] = 1 - beta_by_region[i];
    row_margins[1] = beta_by_region[i];
    joint_prob_matrix = methsim::ipf(ipf_seed, row_margins, col_margins, 1000,
                                     1e-10);
    // Compute p = Pr(Z_{i + 1} = 1 | Z_{i} = z_{i})
    if (Z[i - 1] == 0) {
      p = joint_prob_matrix(0, 1) / (1 - beta_by_region[i - 1]);
    } else {
      p = joint_prob_matrix(1, 1) / beta_by_region[i];
    }
    if (u[i] > p) {
      Z[i] = 0;
    } else{
      Z[i] = 1;
    }

    // Increment j.
    j += 1;
  }
  return Z;
}
