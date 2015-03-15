#include <Rcpp.h>

using namespace Rcpp;

// TODO: Document
//' Sample Z.
//'
//' Simulate sequencing reads by sampling a subset of rows from a given column
//' of H for each read.
//'
//' @param Z an integer matrix where H[i, j] is the methylation state of the
//' i-th methylation locus on the j-th haplotype.
//' @param h an integer vector of length equal to the number of reads. Each
//' element is the haplotype from which each read is to be sampled, i.e.,
//' Z[, j].
//' @param fh an integer vector of length equal to the number of reads. Each
//' element is the position of the first methylation loci ("first hit") from
//' which read is to be sampled, i.e., Z[i, ].
//' @param cqh an integer vector of length equal to the number of reads. Each
//' element is the number of methylation loci that each read overlaps.
//' @param pos an integer vector of length equal to \code{sum(cqh)}. Each
//' element is the genomic position of a methylation locus in Z (ignoring the
//' chromosome, e.g., chr1:77 is simply 77).
//'
//' @keywords internal
//'
//' @return a list of length equal to the number of simulated reads. Each list
//' element is the sequenced/sampled methylation state for that read.
// [[Rcpp::export(".sampleZ")]]
List sampleZ(IntegerMatrix Z, IntegerVector h, IntegerVector fh,
             IntegerVector cqh, IntegerVector start_sm) {

  // Variable checks
  if (Z.nrow() < h.length()) {
    Rcpp::stop("nrow(Z) < length(h).");
  }
  if ((h.length() != fh.length()) ||
      (h.length())!= cqh.length()) {
    Rcpp::stop("length(h) != length(fh) != length(cqh).");
  }
  if (max(h) > Z.ncol()) {
    Rcpp::stop("max(h) < ncol(Z).");
  }
  if (min(cqh) < 1) {
    Rcpp::stop("min(cqh) < 1.");
  }

  // Initialise variables
  int n = h.length();
  List z(n);
  std::vector<int> zz;
  int col_idx;
  int row_idx;
  int pos_idx = 0;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < cqh[i]; j++) {
      row_idx = fh[i] + j - 1; // -1 for C++ indexing
      col_idx = h[i] - 1; // -1 for C++ indexing
      // TODO: Add sampling/copying with error to simulate sequencing error +
      // bisulfite-conversion error.
      zz.push_back(Z(row_idx, col_idx));
      // TODO: Add position information.
      pos_idx += 1;
    }
    // TODO: Add position information to read (e.g., chr1:34,56,99)?
    z[i] = zz;
    zz.clear();
  }

  return z;
}
