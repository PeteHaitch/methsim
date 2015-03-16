#include <Rcpp.h>
using namespace Rcpp;

//' Simulate sequencing errors.
//'
//' Simulate sequencing errors by modifying \code{z} in place.
//'
//' @param z an integer vector of sampled methylation states.
//' @param u an integer vector with the same length as \code{z} of random
//' variables distributed as Uniform(0, 1).
//' @param errorRate the (combined) probability of an "error", where an error
//' may be sequencing error or error in the bisulfite-conversion process.
//'
//' @note \strong{WARNING}: \code{z} is modified \strong{in place}.
//'
//' @return \code{NULL}. \code{z} is modified \strong{in place}.
//'
//' @keywords internal
// [[Rcpp::export(".simErrorInPlace")]]
void simErrorInPlace(IntegerVector z, NumericVector u, double errorRate) {
  // WARNING: This function modifies the input z, i.e., it acts on z
  // **in place**.

  // Argument checks
  if (z.length()!= u.length()) {
    Rcpp::stop("length(z) != length(u)");
  }

  // Variable initialisations
  int n = z.length();

  for (int i = 0; i < n; i++) {
    if (u[i] < errorRate) {
      z[i] = 1 - z[i];
    }
  }

  return;
}
