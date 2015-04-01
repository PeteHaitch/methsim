// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/methsim.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// ipf
arma::mat ipf(const arma::mat& seed, const arma::colvec& row_margins, const arma::rowvec& col_margins, int iter, double tol);
static SEXP methsim_ipf_try(SEXP seedSEXP, SEXP row_marginsSEXP, SEXP col_marginsSEXP, SEXP iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::traits::input_parameter< const arma::mat& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type row_margins(row_marginsSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type col_margins(col_marginsSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    __result = Rcpp::wrap(ipf(seed, row_margins, col_margins, iter, tol));
    return __result;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP methsim_ipf(SEXP seedSEXP, SEXP row_marginsSEXP, SEXP col_marginsSEXP, SEXP iterSEXP, SEXP tolSEXP) {
    SEXP __result;
    {
        Rcpp::RNGScope __rngScope;
        __result = PROTECT(methsim_ipf_try(seedSEXP, row_marginsSEXP, col_marginsSEXP, iterSEXP, tolSEXP));
    }
    Rboolean __isInterrupt = Rf_inherits(__result, "interrupted-error");
    if (__isInterrupt) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean __isError = Rf_inherits(__result, "try-error");
    if (__isError) {
        SEXP __msgSEXP = Rf_asChar(__result);
        UNPROTECT(1);
        Rf_error(CHAR(__msgSEXP));
    }
    UNPROTECT(1);
    return __result;
}
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
IntegerVector simulateZ(NumericVector beta_by_region, NumericVector lor_by_pair, CharacterVector seqnames_one_tuples, NumericVector u);
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
// simulateZ2
std::vector<int> simulateZ2(NumericVector beta_by_region, NumericVector lor_by_pair, CharacterVector seqnames_one_tuples, NumericVector u);
RcppExport SEXP methsim_simulateZ2(SEXP beta_by_regionSEXP, SEXP lor_by_pairSEXP, SEXP seqnames_one_tuplesSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type beta_by_region(beta_by_regionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lor_by_pair(lor_by_pairSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type seqnames_one_tuples(seqnames_one_tuplesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    __result = Rcpp::wrap(simulateZ2(beta_by_region, lor_by_pair, seqnames_one_tuples, u));
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

// validate (ensure exported C++ functions exist before calling them)
static int methsim_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::mat(*ipf)(const arma::mat&,const arma::colvec&,const arma::rowvec&,int,double)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP methsim_RcppExport_registerCCallable() { 
    R_RegisterCCallable("methsim", "methsim_ipf", (DL_FUNC)methsim_ipf_try);
    R_RegisterCCallable("methsim", "methsim_RcppExport_validate", (DL_FUNC)methsim_RcppExport_validate);
    return R_NilValue;
}
