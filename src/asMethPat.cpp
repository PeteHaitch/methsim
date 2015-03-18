#include <Rcpp.h>
using namespace Rcpp;

// WARNING: This will produce zero-rows, i.e., unobserved m-tuples.
// WARNIGN: This should only be used for size >= 2.

// TODO: Document
// TODO: Re-name to better reflect its function.
// TODO: What happens if zero-length vectors are passed to it.
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

  // Initialise the map of all adjacent m-tuples
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
  int idx;
  int idx0 = (pow(2, size) - 1);
  while (k < N) {
    // Check that these 'size' positions are from the same read.
    if (readID[k] == readID[k + size - 1]) {
      // idx converts the methylation pattern, e.g., (1, 1), to the index of
      // the corresponding element in mtuples_value.
      idx = idx0;
      for (int l = k; l < (k + size); l++) {
        idx -= z[l] * pow(2, size - (l - k) - 1);
      }
      mtuples_key << pos[k] << ",";
      for (int m = (k + 1); m < (k + size - 1); m++) {
        mtuples_key << pos[m] << ",";
      }
      mtuples_key << pos[(k + size - 1)];
      mtuples[mtuples_key.str()][idx] += 1;
      mtuples_key.str(std::string());
      k += 1;
    } else {
      // Can jump the k-index because we know if the read doesn't contain a
      // tuple then the next positions in the read can't either.
      k += (size - 1);
    }
  }

  return mtuples;
}
