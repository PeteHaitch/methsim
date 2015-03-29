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

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;


//' Sample W.
//'
//' Sample a column of W for each row of W with the probability of sampling
//' column j for row i given by W[i, j].
//'
//' @param W a numeric matrix where W[i, j] is the probability of sampling
//' column j for row i, i.e., rowMeans(W) == 1.
//'
//' @keywords internal
//'
//' @return an integer vector of length equal to nrow(W), where each element is
//' the column sampled for that row.
// [[Rcpp::export(".sampleW")]]
IntegerVector sampleW(NumericMatrix W) {

  RNGScope scope;

  // Variable initialisations.
  int nc = W.ncol();
  int nr = W.nrow();
  IntegerVector sampled_W(nr);
  IntegerVector W_column_idx = seq_len(nc);
  int size = 1;
  bool replace = false;
  NumericVector prob;

  for (int i = 0; i < nr; i++) {
    prob = W.row(i);
    sampled_W[i] = RcppArmadillo::sample(W_column_idx, size, replace, prob)[0];
  }

  return sampled_W;
}
