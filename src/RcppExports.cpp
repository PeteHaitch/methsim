// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sampleW
IntegerVector sampleW(NumericMatrix W);
RcppExport SEXP methsim_sampleW(SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type W(WSEXP);
    __result = Rcpp::wrap(sampleW(W));
    return __result;
END_RCPP
}
// sampleZ
List sampleZ(IntegerMatrix Z, IntegerVector sampled_W, IntegerVector fh, IntegerVector cqh);
RcppExport SEXP methsim_sampleZ(SEXP ZSEXP, SEXP sampled_WSEXP, SEXP fhSEXP, SEXP cqhSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sampled_W(sampled_WSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type fh(fhSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cqh(cqhSEXP);
    __result = Rcpp::wrap(sampleZ(Z, sampled_W, fh, cqh));
    return __result;
END_RCPP
}
// simErrorInPlace
void simErrorInPlace(IntegerVector z, NumericVector u, double errorRate);
RcppExport SEXP methsim_simErrorInPlace(SEXP zSEXP, SEXP uSEXP, SEXP errorRateSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type errorRate(errorRateSEXP);
    simErrorInPlace(z, u, errorRate);
    return R_NilValue;
END_RCPP
}
// simulateZ
std::vector<int> simulateZ(NumericVector beta_by_region, NumericVector lor_by_pair, CharacterVector seqnames_one_tuples, NumericVector u);
RcppExport SEXP methsim_simulateZ(SEXP beta_by_regionSEXP, SEXP lor_by_pairSEXP, SEXP seqnames_one_tuplesSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type beta_by_region(beta_by_regionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lor_by_pair(lor_by_pairSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type seqnames_one_tuples(seqnames_one_tuplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    __result = Rcpp::wrap(simulateZ(beta_by_region, lor_by_pair, seqnames_one_tuples, u));
    return __result;
END_RCPP
}
// computeP
arma::Mat<double> computeP(NumericVector beta_by_region, NumericVector lor_by_pair, CharacterVector seqnames_one_tuples, int mc_order);
RcppExport SEXP methsim_computeP(SEXP beta_by_regionSEXP, SEXP lor_by_pairSEXP, SEXP seqnames_one_tuplesSEXP, SEXP mc_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type beta_by_region(beta_by_regionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lor_by_pair(lor_by_pairSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type seqnames_one_tuples(seqnames_one_tuplesSEXP);
    Rcpp::traits::input_parameter< int >::type mc_order(mc_orderSEXP);
    __result = Rcpp::wrap(computeP(beta_by_region, lor_by_pair, seqnames_one_tuples, mc_order));
    return __result;
END_RCPP
}
// tabulatez
std::map<std::string, std::vector<int> > tabulatez(IntegerVector readID, IntegerVector z, IntegerVector pos, int size);
RcppExport SEXP methsim_tabulatez(SEXP readIDSEXP, SEXP zSEXP, SEXP posSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type readID(readIDSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pos(posSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    __result = Rcpp::wrap(tabulatez(readID, z, pos, size));
    return __result;
END_RCPP
}
