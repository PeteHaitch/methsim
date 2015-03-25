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
