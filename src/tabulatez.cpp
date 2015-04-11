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

#include <Rcpp.h>
using namespace Rcpp;

//' Tabulate methylation patterns.
//'
//' Tabulate methylation patterns of a given \code{size} from the elements of
//' the \code{z} slot of a \code{\link{SimulatedBS}} object.
//'
//' @param readID an integer vector of read IDs; the \code{readID} column of
//' an element of the \code{z} slot of a \code{\link{SimulatedBS}} object.
//' @param z an integer vector of methylation states; the \code{z} column of
//' an element of the \code{z} slot of a \code{\link{SimulatedBS}} object.
//' @param pos an integer vector of the positions of methylation loci sequenced
//' by each read; the \code{pos} column of an element of the \code{z} slots of
//' a \code{\link{SimulatedBS}} object.
//' @param size an integer greater than 1 specifying the size of the m-tuples
//' to be created.
//'
//' @return A list of tabulated methylation patterns. The name of each list
//' element is a comma-separated string of positions for that m-tuple and the
//' value of each element is a vector of associated counts (the order of this
//' vector is identical to that given by
//' \code{MethylationTuples:::\link[MethylationTuples]{.makeMethPatNames}}).
//'
//' @note \strong{WARNING}: Only adjacent m-tuples, with respect to the
//' sequenced methylation loci, are created. This means that (A) the list will
//' contain unobserved m-tuples (i.e., those with all counts set to zero) and
//' (B) that it will be possibly with paired-end reads to create m-tuples with
//' NIC > 0.
//' \strong{WARNING}: The special case where \code{size} = 1 is handled
//' separately by tabulating with the \code{data.table} package.
//'
//' @keywords internal
//'
// [[Rcpp::export(".tabulatez")]]
std::map<std::string, std::vector<int> > tabulatez(IntegerVector readID,
                                                   IntegerVector z,
                                                   IntegerVector pos,
                                                   int size) {

  // Argument checks.
  if (readID.size() != z.size()) {
    Rcpp::stop("length(readID) != length(z)");
  }
  if (readID.size() != pos.size()) {
    Rcpp::stop("length(readID) != length(pos)");
  }
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
  int n = unique_pos.size() - size + 1;
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
  int N = readID.size() - size + 1;
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
