#include <Rcpp.h>
using namespace Rcpp;

// WARNING: This will produce zero-rows, i.e., unobserved m-tuples.

// [[Rcpp::export(".asMethPat")]]
std::map<std::string, std::vector<int> > asMethPat(IntegerVector readID,
                                           IntegerVector z,
                                           IntegerVector pos,
                                           int size) {

  // TODO: Argument checks
  if (size < 2) {
    Rcpp::stop("size < 2.");
  }

  // Variable initialisations
  // mtuples is a map where the values are the key is the co-ordinates of the
  // m-tuple (pos1,pos2,...,posm) and the values are the counts of each
  // methylation pattern.
  std::map<std::string, std::vector<int> > mtuples;
  std::stringstream mtuples_key;
  std::vector<int> mtuples_value(pow(2, size));

  // TODO: Do I need a corresponding map for the positions themselves?


  // Initialise the map of all possible adjacent m-tuples
  IntegerVector unique_pos = clone(unique(pos));
  std::sort(unique_pos.begin(), unique_pos.end());
  int n = unique_pos.length() - size + 1;
  for (int i = 0; i < n; i++) {
    mtuples_key << unique_pos[i] << ",";
    for (int j = (i + 1); j < (i + size - 1); j++) {
      mtuples_key << unique_pos[j] << ",";
    }
    mtuples_key << unique_pos[(i + size - 1)];
    mtuples[mtuples_key.str()] = mtuples_value;
    mtuples_key.str(std::string());
  }

  // Fill the map with the counts of each methylation pattern
  int k = 0;
  int N = readID.length() - size + 1;
//   Rcpp::Rcout << "N = " << N << std::endl;
  int idx;
  int idx0 = (pow(2, size) - 1);
  while (k < N) {
    // Rcpp::Rcout << "k = " << k << std::endl;
    // Check that these 'size' positions are from the same read.
    if (readID[k] == readID[k + size - 1]) {
      // idx converts the methylation pattern, e.g., (1, 1), to the index of
      // the corresponding element in mtuples_value.
      idx = idx0;
      for (int l = k; l < (k + size); l++) {
        idx -= z[l] * pow(2, size - (l - k) - 1);
      }
      // Rcpp::Rcout << "idx = " << idx << std::endl;
      mtuples_key << pos[k] << ",";
      for (int m = (k + 1); m < (k + size - 1); m++) {
        mtuples_key << pos[m] << ",";
      }
      mtuples_key << pos[(k + size - 1)];
      // Rcpp::Rcout << "mtuples_key = " << mtuples_key.str() << std::endl;
      mtuples[mtuples_key.str()][idx] += 1;
      mtuples_key.str(std::string());
      k += 1;
    } else {
      // TODO: Check this is the correct increment when size > 2
      k += (size - 1);
    }
  }

  // Rcpp::Rcout << "Final k = " << k << std::endl;
  return mtuples;
}
