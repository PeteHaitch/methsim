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

// WARNING: Does not support long vectors.
// TODO: Will to switch index variables from int to something else if
// supporting long vectors.
//' Simulate reads (z)
//'
//' @param fh An integer vector containing the "first hit" for each read.
//' @param nh An integer vector containing the "number of hits" for each read.
//' @param N An integervector giving the number of reads with the same
//' "first hit" and "number of hits".
//' @param marginalProb A matrix of marginal probabilities that each
//' methylation locus is methylated (rows) for each componentn (columns).
//' @param component An integer vector giving the mixture component from
//' which to simulate each read.
//' @param P A list of transition matrices. Each element of P should give the
//' transition probabilities that the methylation locus is methylated given the
//' state of the previous locus. The length of P is equal to the number of
//' mixture components.
//'
//' @return A list(h, readID, z), where h is the "hit", readID is a read ID,
//' and z is the observed methylation state. Each vector has length sum(nh * N).
//'
//' @note Assumes that sum(nh * N) < .Machine$integer.max, because Rcpp cannot
//' yet work with long vectors.
// [[Rcpp::export(".simulatez")]]
Rcpp::List simulatez(IntegerVector fh,
                     IntegerVector nh,
                     IntegerVector N,
                     IntegerVector component,
                     NumericMatrix MarginalProb,
                     ListOf<NumericMatrix> P) {

  // TODO: Do I need to manually "RNGScope scope" or does
  // Rcpp::compileAttributes() do this for me?
  RNGScope scope;

  // Argument checks
  if (fh.size() != nh.size() or fh.size() != N.size() or
        fh.size() != component.size()) {
    std::string error_msg = "length(fh) != length(nh) != length(N) "
      "!= length(component)";
    Rcpp::stop(error_msg);
  }
  if (MarginalProb.ncol() != P.size()) {
    Rcpp::stop("ncol(MarginalProb) != length(P)");
  }
  // TODO: Check that all elements of P have the identical and proper
  // dimensions.
  // The below only checks the first element has the proper nrow.
  if (MarginalProb.nrow() != P[1].nrow()) {
    Rcpp::stop("length(marginalProb) != nrow(P[[1]])");
  }
  // Make sure we don't try to access out of bounds elements of P.
  if (max(fh + nh - 1) > P[1].nrow()) {
    Rcpp::stop("max(fh + nh - 1) > nrow(P)");
  }

  // Variable initialisations
  // nr is the number of reads
  int nr = fh.size();
  // n is the length of the output vectors
  int n = sum(nh * N);
  IntegerVector h(n);
  IntegerVector readID(n);
  IntegerVector z(n);
  // l indexes h, readID and z.
  int l = 0;
  // rid is the current value of readID
  int rid = 0;

  // Loop over fh (i) and simulate N (j) paths of length nh (k).
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < N[i]; j++) {
      h[l] = fh[i];
      readID[l] = rid;
      // Indices of MarginalProb have -1 because C++ is 0-indexed.
      z[l] = R::rbinom(1, MarginalProb(fh[i] - 1, component[i] - 1));
      l += 1;
      for (int k = 1; k < nh[i]; k++) {
        h[l] = fh[i] + k;
        readID[l] = rid;
        // Indices of P have -1 because C++ is 0-indexed.
        z[l] = R::rbinom(1, P[component[i] - 1](fh[i] + k - 1, z[l - 1]));
        l += 1;
      }
      rid += 1;
    }
  }

  return List::create(_["h"] = h, _["readID"] = readID, _["z"] = z);
}
