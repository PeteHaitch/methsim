// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;


//' Sample H.
//'
//' Sample a column of H for each row of H with the probability of sampling
//' column j for row i given by H[i, j].
//'
//' @param H a numeric matrix where H[i, j] is the probability of sampling
//' column j for row i, i.e., rowMeans(H) == 1.
//'
//' @keywords internal
//'
//' @return an integer vector of length equal to nrow(H), where each element is
//' the column sampled for that row.
// [[Rcpp::export(".sampleH")]]
IntegerVector sampleH(NumericMatrix H) {

  RNGScope scope;

  // Variable initialisations.
  int nc = H.ncol();
  int nr = H.nrow();
  IntegerVector h(nr);
  IntegerVector hi = seq_len(nc);
  int size = 1;
  bool replace = false;
  NumericVector prob;

  for (int i = 0; i < nr; i++) {
    prob = H.row(i);
    h[i] = RcppArmadillo::sample(hi, size, replace, prob)[0];
  }

  return h;
}
